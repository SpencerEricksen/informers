#!/home/ssericksen/anaconda2/bin/python2.7

import baseline_inf_sel_rank_methods_v23 as blm
import sys
import pandas as pd
import numpy as np


try:
    n_informers = int(sys.argv[1])   # 16
    matrix_dataset = str(sys.argv[2])       # 1 or 2
    newtarg_dataset = str(sys.argv[3])      # datafile
except:
    print('')
    print("usage: ./script.py num_informers matrix_dataset newtarg_dataset")
    print('       arg0        arg1          arg2           arg3           ')
    print('')
    print(" running with defaults:")
    print("		n_informers	=16")
    print("		matrix_dataset	='1'")
    print("		newtarg_dataset	='./data/data_pkis1.csv' ")
    print('')
    #n_informers     = 16           		   # any integer
    #matrix_dataset  = '1'  			   # also: 'pkis2'
    #newtarg_dataset = './data/data_pkis1.csv'  # also: 'BGLF4'
    exit(1)

# load data
if matrix_dataset == '1':
    activity_data_csv = '../data/pkis1.csv'
    fps_file = '../data/pkis1_uniq_ecfp4_1024.fps'
    distance_data_csv = '../data/pkis1_jaccard_distmat_ecfp4.csv'
elif matrix_dataset == '2':
    activity_data_csv = '../data/pkis2.csv'
    fps_file = "../data/pkis2_uniq_415_ecfp4_1024.fps"
    distance_data_csv = '../data/pkis2_415_jaccard_distmat_ecfp4.csv'
else:
    print('matrix dataset (arg4) was invalid')
    exit(1)


ranked_sets = []

print('targ, n_inf, inf_sel, ranking, matrix, est_thresh, num_act_inf')
for targ in [ 'bglf4', 'pknb' ]:
    for inf_selection in [ 'c', 'p' ]:
        for ranking in [ 's', 'l', 'w' ]:
       
            # read in data matrix (PKIS1 or PKIS2)
            df_act_mat = blm.get_continuous( activity_data_csv )

            # get binary data matrix (labels real, not inferred)
            df_binary = blm.get_binary( df_act_mat )
            
            # get informers
            if inf_selection == 'c':
                inf_molids = blm.inf_sel_clst_medoids( n_informers, fps_file )
            elif inf_selection == 'p':
                inf_molids = blm.inf_sel_max_prom_global( n_informers, df_binary )
            else:
                print('informer selection method, {}, is invalid.').format( inf_selection )
                exit(1)
            
            # get the thresholds
            df_thresh = pd.read_csv('../data/thresholds_2sigma/newtarget_thresholds.csv', index_col='dataset')
            thresh = df_thresh[targ+'_'+inf_selection].loc['PKIS'+matrix_dataset]
            
            # read in new target data
            df_targ_act = pd.read_csv( newtarg_dataset, index_col='molid' )
            df_targ_act.index = df_targ_act.index.map(str)            

            # get informer activities
            inf_activities = list( df_targ_act[targ].loc[ inf_molids ] )
            
            # get list of just the active informers--those that pass activity threshold
            act_inf_molids = blm.process_informer_list( inf_molids, inf_activities, thresh=thresh )
            
            # if no active informers and method was active-informer dependent, assign np.nans for scores
            if len(act_inf_molids) == 0 and ranking in ('s','l'):
                noninf_molids = list( set(df_act_mat.index) ^ set(inf_molids) )
                ranked_noninf_molids = pd.Series( [np.nan]*len(noninf_molids), index=noninf_molids )
            else:
                # get noninf-by-inf distance submatrix
                df_dist_submat = blm.get_distance_submatrix_df( distance_data_csv, inf_molids, act_inf_molids )

                # choose ranking method, score noninformers, build rank-ordered list
                if ranking == 's':
                            # simple expansion
                            ranked_noninf_molids = blm.rank_by_active_expansion( df_dist_submat )
                elif ranking == 'l':
                            # loop expansion
                            ranked_noninf_molids = blm.rank_by_loop_active_expansion( df_dist_submat )
                elif ranking == 'w':
                            # weighted expansion
                            df_dist_submat = blm.get_distance_submatrix_df( distance_data_csv, inf_molids, inf_molids )
                            ranked_noninf_molids = blm.rank_by_weighted_expansion( inf_activities, df_dist_submat )
                else:
                            print('missing or invalid selection method indicated')
                            exit(1)
            
            inf_series = pd.Series( ['informer']*len(inf_molids), index=inf_molids )
            ranked_molids = pd.concat( [ranked_noninf_molids, inf_series ], axis=0 )
            ranked_molids.rename( '{}_{}_{}_{}'.format(  targ, inf_selection, ranking, matrix_dataset ), inplace=True )
            ranked_sets.append( ranked_molids )
            print('{}, {}, {}, {}, {}, t{:.3f}, {}').format( targ, str(n_informers), inf_selection, ranking, matrix_dataset, thresh, len(act_inf_molids) )

ranked_df = pd.concat( ranked_sets, axis=1, sort=False )
ranked_df.to_csv( 'rankings_baseline_pkis'+matrix_dataset+'.csv', na_rep=np.nan, index_label='molid' )

