#!/home/ssericksen/anaconda2/bin/python2.7

import numpy as np
import pandas as pd
import informer_functions as inf
import sys
    
activity_data_csv = '../data/pkis1.csv'
distance_data_csv = '../data/compounds/pkis1_jaccard_distmat_ecfp6.csv'
fps_file = '../data/compounds/pkis1_uniq_ecfp6_1024.fps'
 
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
df_act_mat = inf.get_continuous( activity_data_csv )

# get binary data matrix
df_binary = inf.get_binary( df_act_mat )

# read in thresholds data (provided by Huikun--predicted \
# 2sigma thresholds based on informer activities)
if inf_selection == 'prom':
    df_thresh = pd.read_csv( '../data/thresholds_2sigma/prom_le_thresh.txt', \
            index_col=0, header=None, names=['Activity'], delimiter=" " )
elif inf_selection == 'clst':
    df_thresh = pd.read_csv( '../data/thresholds_2sigma/clst_le_thresh.txt', \
            index_col=0, header=None, names=['Activity'], delimiter=" " )    
else:
    print('issue with 2sigma thresholds!')
    exit(1)

# obtain list of informer cpds via different informer selection methods
if inf_selection == 'clst':
    inf_molids = inf.inf_sel_clst_medoids( n_informers, fps_file )
elif inf_selection == 'prom':
    # this is a little tricky--LOTO causes slight variation in maximally promiscuous set
    # therefore, a dictionary of targets with associated max prom informers is provided
    inf_molids_dict = inf.inf_sel_max_prom( n_informers, df_binary )
    # just use the informer set provided for the target of interest
    #inf_molids = inf_molids_dict[targ]
elif inf_selection == 'prom_glob':
    inf_molids = inf.inf_sel_max_prom_global( n_informers, df_binary )
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

    '''Use predicted activity thresholds for target based on informer activity data.
       Active/Inactive threshold used in our work is z-score >= +2, but without knowledge
       of cpd activity distribution on target, activity at z-score=2 is estimated using
       informer activities and inference to activity distributions of nearest off-targets 
       in PKIS1 or PKIS2 matrices.'''
    
    # ground truth threshold
    #thresh = df_act_mat[targ].mean(axis=0) + 2*df_act_mat[targ].std(axis=0)

    # predicted threshold
    thresh = float( df_thresh.loc[targ] )

    # get list of the active informers--remove informers that don't satify threshold
    act_inf_molids = inf.process_informer_list( inf_molids, inf_activities, thresh=thresh )
    print("{},{}").format( targ, len(act_inf_molids) )
    if len(act_inf_molids) == 0 and ranking in ('se','le'):
        default_thresh = max(inf_activities) # if no actives in informer set use the most active cpd
        act_inf_molids = inf.process_informer_list( inf_molids, inf_activities, thresh=default_thresh )

    # get noninf-by-inf distance submatrix
    df_dist_submat = inf.get_distance_submatrix_df( distance_data_csv, inf_molids, act_inf_molids )
    
    df_dist_submat_temp = df_dist_submat.copy()
   
    # choose ranking method, score noninformers, build rank-ordered list
    if ranking == 'se':
        # simple expansion
        ranked_noninf_molids = inf.rank_by_active_expansion( df_dist_submat_temp )
    elif ranking == 'le':
        # loop expansion
        ranked_noninf_molids = inf.rank_by_loop_active_expansion( df_dist_submat_temp )
    elif ranking == 'we':
        # weighted expansion
        df_dist_submat = inf.get_distance_submatrix_df( distance_data_csv, inf_molids, inf_molids )
        ranked_noninf_molids = inf.rank_by_weighted_expansion( inf_activities, df_dist_submat )
    else:
        print('missing or invalid selection method indicated')
        exit(1)
   
    inf_series = pd.Series( ['informer']*len(inf_molids), index=inf_molids )
    ranked_molids = pd.concat( [ranked_noninf_molids, inf_series ], axis=0 )
    ranked_molids.rename(targ, inplace=True)
    ranked_sets.append( ranked_molids )

ranked_df = pd.concat( ranked_sets, axis=1, sort=False )
ranked_df.to_csv( '../output_pkis1loto/rankings/bl/ranked_df_'+str(n_informers)+'_'+inf_selection+'_'+ranking+'_thresh_ecfp6_v1.2.csv', index_label='molid' )

