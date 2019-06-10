#!/home/ssericksen/anaconda2/bin/python2.7

import informer_functions as inf
import sys
import pandas as pd
import numpy as np

try:
    n_informers = int(sys.argv[1])   # 16
    matrix_dataset = str(sys.argv[2])       # 1 or 2
except:
    print("")
    print("usage: ./script.py num_informers matrix_dataset ")
    print("       arg0        arg1          arg2           ")
    print("")
    print("		num_informers: 16")
    print("		matrix_dataset	= '1' for pkis1 or '2' for pkis2")
    print("")
    exit(1)

# load data
if matrix_dataset == '1':
    #newtarg_dataset = '../data/data_newtargs_pkis1cpds.csv'
    activity_data_csv = '../data/pkis1.csv'
    fps_file = '../data/compounds/pkis1_uniq_ecfp4_1024.fps'
    distance_data_csv = '../data/compounds/pkis1_jaccard_distmat_ecfp4.csv'
elif matrix_dataset == '2':
    #newtarg_dataset = '../data/data_newtargs_pkis2cpds.csv'
    activity_data_csv = '../data/pkis2.csv'
    fps_file = "../data/compounds/pkis2_uniq_415_ecfp4_1024.fps"
    distance_data_csv = '../data/compounds/pkis2_415_jaccard_distmat_ecfp4.csv'
else:
    print('need matrix dataset (arg4) was invalid')
    exit(1)


ranked_sets = []

print('targ, n_inf, inf_sel, ranking, matrix, est_thresh, num_act_inf')
#for targ in [ 'bglf4', 'pknb' ]:
for targ in [ 'rop18' ]:
    for inf_selection in [ 'c', 'p' ]:
        for ranking in [ 's', 'l', 'w' ]:
       
            # read in data matrix (PKIS1 or PKIS2)
            df_act_mat = inf.get_continuous( activity_data_csv )

            # get informers
            df_inf = pd.read_excel('../data/rop18/rop18_baseline_informer_activities_2019-05-15.xlsx', sheet_name='pkis'+matrix_dataset+'_'+inf_selection, index_col='molid')
            inf_molids = df_inf.index.map(str)
            inf_activities = df_inf['activity'].to_list()

            thresh = 0.50            

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
ranked_df.to_csv( '../output_newtargs/bl/rankings_3_baseline_pkis'+matrix_dataset+'test.csv', na_rep=np.nan, index_label='molid' )

