#!/home/ssericksen/anaconda2/bin/python2.7

import numpy as np
import pandas as pd
from sklearn import metrics

############################################################################################
#####   informer selection procedures   ####################################################
############################################################################################
#
# cluster medoids	=	take as informers the medoids from agglomerative hierarchical 
#				clustering of Jaccard distance matrix of ECFP4 for compounds
#
# max promiscuous	=	for each left-out target take 16 most promiscuous cpds across 
#				the other 223 targets in matrix 
#
############################################################################################

def get_continuous( activity_matrix_file ):
    '''read in activity data for matrix (Huikun's space delimited pkis1.csv),
       return a dataframe with cpd molid indices and target columns

       note: the data file is transposed'''
    df = pd.read_csv( activity_matrix_file, delimiter=" ", index_col=0 ).T
    df.index = df.index.map(str)
    return df


def get_binary( df ):
    ''' input continous activity data frame, return binary activity dataframe '''
    df_binary = df[ df > (df.mean(axis=0) + 2*df.std(axis=0)) ].notnull()
    return df_binary


#### informer selection methods ####

def inf_sel_max_prom( n_informers, df_binary ):
    '''given n_informers and df_binary, return a dictionary where each key is kinase target
       and its value a list of maximally promiscuous cpd molids to be taken as informer set

       note: the maximally promiscuous set for a kinase is determined with that kinase held out'''
    prom_inf_dict = {}
    # iterate through targets, remove targ, determine 16 most promiscuous
    for targ in df_binary.columns:
        inf_list = df_binary.drop( columns=[ targ ]).sum(axis=1).sort_values(ascending=False, na_position='last' ).head( n_informers ).index
        prom_inf_dict[targ] = list(inf_list)
    return prom_inf_dict


def inf_sel_max_prom_global( n_informers, df_binary ):
    ''' given n_informers and df_binary, return a list of maximally promiscuous cpds to be
        used as the promiscuous informer set. Note, rather than LOTO, a one-time global 
        list is generated '''
    inf_list = df_binary.sum(axis=1).sort_values(ascending=False, na_position='last' ).head( n_informers ).index
    return inf_list


def inf_sel_clst_medoids( n_informers, fps_file ):
    '''given n_informers and a CSV fingerprints file (binary strings), return
       an informer list of cpd molids as the cluster medoids

       note: agglomerative hierarchical clustering is applied with average linkage
             and a distance matrix is computed using Jaccard distances between fingerprints'''
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.neighbors import DistanceMetric
    from itertools import combinations

    # read in CSV fingerprints to dataframe
    df = pd.read_csv( fps_file, index_col=0 )
    df.index = df.index.map(str)
    # split up bitstrings into integer bit vectors
    df['features'] = df['fingerprint'].apply( lambda x: np.array(list(map(int,list(x)))))

    tup = tuple( df['features'].values )
    # make feature matrix (rows:examples, columns:features)
    X = np.vstack( tup )

    # get mol identifiers "molids"
    molids = tuple(df.index)
    n_cpds = len(molids)

    # get the Jaccard distance matrix
    dist = DistanceMetric.get_metric('jaccard')
    dist_mat = dist.pairwise(X)

    # run Agglomerative Clustering with Average linkage using precomputed dist matrix
    model = AgglomerativeClustering(linkage='average', affinity='precomputed', connectivity=None, n_clusters=n_informers)
    model = model.fit(dist_mat)
    # note, sklearn will not allow k-means with precomputed dist--it uses Euclidean distances between points
    # not sure how this compares to Jaccard distances for our ECFP4 features (1024 element bit vectors)

    # get medoids (most central structure in each cluster)
    dist_dict = {}
    cids = list(set(model.labels_))

    for cid in cids:
        member_cpds = np.where( model.labels_ == cid )[0]
        cid_dist_matrix = np.zeros( (n_cpds, n_cpds) )
        for pair in combinations( member_cpds, 2 ):
            cid_dist_matrix[ pair[0], pair[1] ] = dist_mat[ pair[0], pair[1] ]

        # get dist sums
        cid_dist_sums = ( cid_dist_matrix.sum( axis=0 ) + cid_dist_matrix.sum( axis=1 ) ) / float(len(member_cpds))
        for cpd, dist_sum in zip( member_cpds, cid_dist_sums[member_cpds] ):
            # store clustering info for each frame
            if cid not in dist_dict:
                dist_dict[cid] = [ (cpd, dist_sum) ]
            else:
                dist_dict[cid].append( (cpd, dist_sum) )

        dist_dict[cid] = sorted( dist_dict[cid], key=lambda x:x[1] )

    # find representative frame for each cid
    rep_frames = [ dist_dict[cid][0][0] for cid in dist_dict ]

    inf_list = []
    for cid in dist_dict:
        inf_list.append( molids[dist_dict[cid][0][0]] )
    return inf_list


############################################################################################
#####   noninformer ranking procedures   ###################################################
############################################################################################

# simple expansion      =       "active expansion" -- rank noninformers according to minimal
#                               distance to any active informer

# loop expansion    =       "loop active expansion" -- rank noninformers by looping through 
#                               active informers for each informer, append nearest non-ranked
#                               non-informer to rank list 

# weighted expansion  =       rank noninformers by sum of activity-weighted similarity 
#                               vectors to informer set

############################################################################################

def gen_rand_inf_act_vec( lb, ub, k ):
    '''generate a random activity vector of percentage inhibition values
       with lb=lowerbound(float), ub=upperbound, and k (integer) the number 
       of values (# informers)'''
    # generate vector random activites (% inhib) with flat dist
    A = np.random.uniform( low=lb, high=ub, size=( k, ) )
    return A

def process_informer_list( inf_molid_list, inf_activity_list, thresh=None ):
    '''read in full informer list and return list of [active] informers in
       order of descending activity, threshold may be input as an option'''
    # get activity sorted informer list for only the "active" informers
    inf_tups = zip( inf_molid_list, inf_activity_list )
    inf_tups = sorted( inf_tups, reverse=True, key=lambda x: x[1] )
    if thresh != None:
        inf_tups = [ i for i in inf_tups if i[1] >= float(thresh) ]
    new_inf_list = [ x[0] for x in inf_tups ]
    return new_inf_list 

def get_distance_submatrix_df( distmat_csv, raw_inf_list, active_inf_list ):
    '''read in full distance matrix (csv), full informer list, and active informer
       list, then return dataframe with index=noninformers and columns=informers'''
    df = pd.read_csv( distmat_csv, index_col=0 )
    df.index = df.index.map(str)
    # to make rows noninformers and columns informers:
    df = df[ active_inf_list ]
    df.drop( index=raw_inf_list, inplace=True )
    # return noninf-by-inf distance submatrix as dataframe
    return df

#### ranking procedures ####

def rank_by_loop_active_expansion( df_dist_submat ):
    ''' procedure for NN_loop_active_informers (balanced expansion), input distance
        submatrix, return ranked list of molids (descending in priority) '''
    ranked_noninf_molids = []
    act_inf_molids = list( df_dist_submat.columns )
    # while there are still unranked noninformers out there, keep looping through informers
    while len(df_dist_submat.index) > 0:
        for c in act_inf_molids:
            # check to see if noninformers are still neighbors out there for informer c, find and append to ranked list
            if len( df_dist_submat[c] > 0 ):
                # find row (noninformer molid) with minimal distance to informer column
                NN_molid = df_dist_submat[c].idxmin(axis=0)
                ranked_noninf_molids.append(NN_molid)
                # make sure to delete this noninformer since it has been added to ranked list
                df_dist_submat.drop( index=NN_molid, inplace=True )
    
    d = { i:ranked_noninf_molids.index(i) for i in ranked_noninf_molids }
    s = pd.Series(d)
    return s

def rank_by_active_expansion( df_dist_submat ):
    ''' procedure for NN_active_informers (hit expansion), input distance
        submatrix, return series of noninformers scored by min distance to any active informer '''
    inf_molids = list(df_dist_submat.columns)
    s = df_dist_submat.min(axis=1)
    return s

def rank_by_weighted_expansion( inf_activities, df_dist_submat ):
    ''' rank noninformers by sum of activity-weighted similarities to all informers,
        input list of informer activities and distance submatrix, return ranked list of
        molids (descending in priority) ''' 
    ranked_noninf_molids = []
    A_norm = np.array(inf_activities)
    df_sim_submat = 1.00 - df_dist_submat
    s = -1.00 * df_sim_submat.apply( lambda x: np.dot( A_norm, np.array(x) ), axis=1 )
    return s

###########################################
# metrics #
###########################################

def compute_roc_auc( labels_arr, scores_arr ):
    '''use an sklearn function to compute ROC AUC
        probably should add some other metrics to this'''
    if len(np.unique(labels_arr)) == 2:
        auc = metrics.roc_auc_score( labels_arr, scores_arr )
    else:
        auc = 'ND'
    return auc

def enrichment_factor_single(labels_arr, scores_arr, percentile):
    '''
    calculate the enrichment factor based on some upper fraction
    of library ordered by docking scores. upper fraction is determined
    by percentile (actually a fraction of value 0.0-1.0)
    -1 represents missing value
    and remove them when in evaluation
    '''
    non_missing_indices = np.argwhere(labels_arr != -1)[:, 0]
    labels_arr = labels_arr[non_missing_indices]
    scores_arr = scores_arr[non_missing_indices]

    sample_size = int(labels_arr.shape[0] * percentile)           # determine number mols in subset
    pred = np.sort(scores_arr, axis=0)[::-1][:sample_size]        # sort the scores list, take top subset from library
    indices = np.argsort(scores_arr, axis=0)[::-1][:sample_size]  # get the index positions for these in library
    n_actives = np.nansum(labels_arr)                             # count number of positive labels in library
    total_actives = np.nansum(labels_arr)
    total_count = len(labels_arr)
    n_experimental = np.nansum(labels_arr[indices])               # count number of positive labels in subset
    temp = scores_arr[indices]

    if n_actives > 0.0:
        ef = float(n_experimental) / n_actives / percentile       # calc EF at percentile
        ef_max = min(n_actives, sample_size) / (n_actives * percentile)
    else:
        ef = 'ND'
        ef_max = 'ND'
    return n_actives, ef, ef_max

def normalized_enrichment_factor_single_hz(labels_arr, scores_arr, percentile):
    '''this is an adjusted NEF that enables better comparison across targets with wide
       variation in class imbalance--this was the metric implement in our informer paper'''
    n_actives, ef, ef_max = enrichment_factor_single(labels_arr, scores_arr, percentile)
    return ( 1.00 + (ef - 1.00) / (ef_max - 1.00) ) / 2.00

