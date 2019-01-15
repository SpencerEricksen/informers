#!/home/ssericksen/anaconda2/bin/python2.7

import numpy as np
import pandas as pd

############################################################################################
#####   informer selection procedures   ####################################################
############################################################################################

# cluster medoids	=	take as informers the medoids from agglomerative hierarchical 
#				clustering of Jaccard distance matrix of ECFP4 for compounds

# max promiscuous	=	for each left-out target take 16 most promiscuous cpds across 
#				the other 223 targets in matrix 

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
            #print cid, cpd, dist_sum

        dist_dict[cid] = sorted( dist_dict[cid], key=lambda x:x[1] )

    # find representative frame for each cid
    rep_frames = [ dist_dict[cid][0][0] for cid in dist_dict ]

    inf_list = []
    for cid in dist_dict:
        #tup = tuple( [ cid, molids[dist_dict[cid][0][0]], len(dist_dict[cid]) ] )
        #inf_list.append(tup)
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
    #A = np.array( [ 98.5, 10.4, 33.3, 6.2, 0.1, 0.3 ] )
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
                #dist = df_dist_submat[c].loc[NN_molid]
                #tup = (NN_molid, c, dist)
                #ranked_noninf_molids.append( tup )
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
    #df_dist_submat['min_dist'] = df_dist_submat.min(axis=1)
    #df_dist_submat['min_inf'] = df_dist_submat[inf_molids].idxmin(axis=1)
    #temp = df_dist_submat[ [ 'min_inf', 'min_dist'] ].sort_values( by='min_dist', ascending=True )
    #temp['molid'] = temp.index
    #ranked_noninf_molids = [ tuple(c) for c in temp[ ['molid', 'min_inf', 'min_dist'] ].values ]
    #return ranked_noninf_molids
    return s

def rank_by_weighted_expansion( inf_activities, df_dist_submat ):
    ''' rank noninformers by sum of activity-weighted similarities to all informers,
        input list of informer activities and distance submatrix, return ranked list of
        molids (descending in priority) ''' 
    ranked_noninf_molids = []
    #A_norm = np.array(inf_activities) / 100.0
    A_norm = np.array(inf_activities)
    df_sim_submat = 1.00 - df_dist_submat
    s = -1.00 * df_sim_submat.apply( lambda x: np.dot( A_norm, np.array(x) ), axis=1 )
    #df_sim_submat['WE_score'] = df_sim_submat.apply( lambda x: np.dot( A_norm, np.array(x) ), axis=1 )
    #df_sim_submat.sort_values( by='WE_score', ascending=False, inplace=True )
    #df_sim_submat['molid'] = df_sim_submat.index
    #ranked_noninf_molids = [ tuple(c) for c in df_sim_submat[ ['molid', 'WE_score'] ].values ]
    #return ranked_noninf_molids
    return s

###################################################################################

def main():

    import sys
    
    activity_data_csv = '../data/pkis1.csv'
    distance_data_csv = '../data/compounds/pkis1_jaccard_distmat_ecfp4.csv'
    fps_file = '../data/compounds/pkis1_uniq_ecfp4_1024.fps'
 
    try:
        n_informers = int(sys.argv[1])
        inf_selection = str(sys.argv[2])
        ranking = str(sys.argv[3])
    except:
        print('')
        print("usage: ./script.py   number_informers: 16  inf_sel_methods:'clst','prom','prom_glob'  ranking:'se','le','we' ") 
        print('       arg0          arg1                  arg2                                       arg3')  
        print('')
        exit(1) 

    # read in data matrix
    df_act_mat = get_continuous( activity_data_csv )
    
    # get binary data matrix
    df_binary = get_binary( df_act_mat )

    # read in thresholds data (provided by Huikun--predicted 2sigma thresholds based on informer activities)
    if inf_selection == 'prom':
        df_thresh = pd.read_csv( '../data/thresholds_2sigma/prom_le_thresh.txt', index_col=0, header=None, names=['Activity'], delimiter=" " )
    elif inf_selection == 'clst':
        df_thresh = pd.read_csv( '../data/thresholds_2sigma/clst_le_thresh.txt', index_col=0, header=None, names=['Activity'], delimiter=" " )    
    else:
        print('issue with 2sigma thresholds!')
        exit(1)

    # obtain list of informer cpds via different informer selection methods
    if inf_selection == 'clst':
        inf_molids = inf_sel_clst_medoids( n_informers, fps_file )
    elif inf_selection == 'prom':
        # this is a little tricky--LOTO causes slight variation in maximally promiscuous set
        # therefore, a dictionary of targets with associated max prom informers is provided
        inf_molids_dict = inf_sel_max_prom( n_informers, df_binary )
        # just use the informer set provided for the target of interest
        #inf_molids = inf_molids_dict[targ]
    elif inf_selection == 'prom_glob':
        inf_molids = inf_sel_max_prom_global( n_informers, df_binary )
    else:
        print("no valid informer selection method provided!")
        exit(1)
    
    ranked_sets = []
    
    print("target,n_active_infs")
    for targ in df_act_mat.columns:
        if inf_selection == 'prom':
            inf_molids = inf_molids_dict[targ]
    
        # get the activities for the informer set on target of interest
        inf_activities = list( df_act_mat[targ].loc[ inf_molids ] )
     
        # determine activity threshold for target based on its activity data
        #thresh = df_act_mat[targ].mean(axis=0) + 2*df_act_mat[targ].std(axis=0)
        thresh = float( df_thresh.loc[targ] )

        # get list of the active informers--remove informers that don't satify threshold
        act_inf_molids = process_informer_list( inf_molids, inf_activities, thresh=thresh )
        print("{},{}").format( targ, len(act_inf_molids) )
        if len(act_inf_molids) == 0 and ranking in ('se','le'):
            noninf_molids = list( set(df_act_mat.index) ^ set(inf_molids) ) 
            ranked_noninf_molids = pd.Series( [np.nan]*len(noninf_molids), index=noninf_molids )
    
        else:
            # get noninf-by-inf distance submatrix
            df_dist_submat = get_distance_submatrix_df( distance_data_csv, inf_molids, act_inf_molids )
        
            df_dist_submat_temp = df_dist_submat.copy()
        
            # choose ranking method, score noninformers, build rank-ordered list
            if ranking == 'se':
                # simple expansion
                ranked_noninf_molids = rank_by_active_expansion( df_dist_submat_temp )
            elif ranking == 'le':
                # loop expansion
                ranked_noninf_molids = rank_by_loop_active_expansion( df_dist_submat_temp )
            elif ranking == 'we':
                # weighted expansion
                df_dist_submat = get_distance_submatrix_df( distance_data_csv, inf_molids, inf_molids )
                ranked_noninf_molids = rank_by_weighted_expansion( inf_activities, df_dist_submat )
            else:
                print('missing or invalid selection method indicated')
                exit(1)
       
        inf_series = pd.Series( ['informer']*len(inf_molids), index=inf_molids )
        ranked_molids = pd.concat( [ranked_noninf_molids, inf_series ], axis=0 )
        ranked_molids.rename(targ, inplace=True)
        ranked_sets.append( ranked_molids )
    
    
    ranked_df = pd.concat( ranked_sets, axis=1, sort=False )
    ranked_df.to_csv( 'ranked_df_'+str(n_informers)+'_'+inf_selection+'_'+ranking+'_thresh_v24.csv', index_label='molid' )
    #print ranked_df

if __name__ == "__main__":
    main()

