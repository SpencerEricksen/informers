import informer_functions as inf
import sys
import pandas as pd
import numpy as np

try:
    n_informers = int(sys.argv[1])   # 16
    matrix_dataset = str(sys.argv[2])       # 1 or 2
except:
    print("")
    print("usage: python script.py num_informers matrix_dataset ")
    print("                arg0        arg1          arg2           ")
    print("")
    print("		num_informers: 16")
    print("		matrix_dataset	= '1' for pkis1 or '2' for pkis2")
    print("")
    exit(1)

# load data
if matrix_dataset == '1':
    newtarg_dataset = '../data/data_newtargs_pkis1cpds.csv'
    activity_data_csv = '../data/pkis1.csv'
    fps_file = '../data/compounds/pkis1_uniq_ecfp4_1024.fps'
    distance_data_csv = '../data/compounds/pkis1_jaccard_distmat_ecfp4.csv'
elif matrix_dataset == '2':
    newtarg_dataset = '../data/data_newtargs_pkis2cpds.csv'
    activity_data_csv = '../data/pkis2.csv'
    fps_file = "../data/compounds/pkis2_uniq_415_ecfp4_1024.fps"
    distance_data_csv = '../data/compounds/pkis2_415_jaccard_distmat_ecfp4.csv'
else:
    print('matrix dataset (arg4) was invalid')
    exit(1)


ranked_sets = []

print('targ, n_inf, inf_sel, ranking, matrix, est_thresh, num_act_inf')
for targ in [ 'bglf4', 'pknb', 'rop18' ]:
    for inf_selection in [ 'c', 'p' ]:
        for ranking in [ 's', 'l', 'w' ]:
       
            # read in data matrix (PKIS1 or PKIS2)
            df_act_mat = inf.get_continuous( activity_data_csv )

            # get binary data matrix (labels real, not inferred)
            df_binary = inf.get_binary( df_act_mat )
            
            # get informers
            if inf_selection == 'c':
                inf_molids = inf.inf_sel_clst_medoids( n_informers, fps_file )
            elif inf_selection == 'p':
                inf_molids = inf.inf_sel_max_prom_global( n_informers, df_binary )
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
            act_inf_molids = inf.process_informer_list( inf_molids, inf_activities, thresh=thresh )
            
            # if no active informers and method was active-informer dependent, assign np.nans for scores
            if len(act_inf_molids) == 0 and ranking in ('s','l'):
                thresh = max(inf_activities)
                act_inf_molids = inf.process_informer_list( inf_molids, inf_activities, thresh=thresh )

            # get noninf-by-inf distance submatrix
            df_dist_submat = inf.get_distance_submatrix_df( distance_data_csv, inf_molids, act_inf_molids )

            # choose ranking method, score noninformers, build rank-ordered list
            if ranking == 's':
                # simple expansion
                ranked_noninf_molids = inf.rank_by_active_expansion( df_dist_submat )
            elif ranking == 'l':
                # loop expansion
                ranked_noninf_molids = inf.rank_by_loop_active_expansion( df_dist_submat )
            elif ranking == 'w':
                # weighted expansion
                df_dist_submat = inf.get_distance_submatrix_df( distance_data_csv, inf_molids, inf_molids )
                ranked_noninf_molids = inf.rank_by_weighted_expansion( inf_activities, df_dist_submat )
            else:
                print('missing or invalid selection method indicated')
                exit(1)
            
            inf_series = pd.Series( ['informer']*len(inf_molids), index=inf_molids )
            ranked_molids = pd.concat( [ranked_noninf_molids, inf_series ], axis=0 )
            ranked_molids.rename( '{}_{}_{}_{}'.format(  targ, inf_selection, ranking, matrix_dataset ), inplace=True )
            ranked_sets.append( ranked_molids )
            print('{}, {}, {}, {}, {}, t{:.3f}, {}').format( targ, str(n_informers), inf_selection, ranking, matrix_dataset, thresh, len(act_inf_molids) )

ranked_df = pd.concat( ranked_sets, axis=1, sort=False )
ranked_df.to_csv( '../output_newtargs/bl/rankings_baseline_pkis'+matrix_dataset+'.csv', na_rep=np.nan, index_label='molid' )

