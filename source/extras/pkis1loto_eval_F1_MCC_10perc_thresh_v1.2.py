#!/home/ssericksen/anaconda2/bin/python2.7

import numpy as np
import pandas as pd
import informer_functions as inf
import glob
import sys
import sklearn as sk

rankings_files = glob.glob('../output_pkis1loto/rankings/bl/ranked_df_16_*.csv')
rankings_files = rankings_files + glob.glob('../output_pkis1loto/rankings/hz/*.csv')
rankings_files = rankings_files + glob.glob('../output_pkis1loto/rankings/cp/*.csv')
# alphabetize the file list so we can re-arrange later
rankings_files.sort()

activity_matrix_file = '../data/pkis1.csv'
df_continuous = inf.get_continuous( activity_matrix_file )
df_binary = inf.get_binary( df_continuous )

s_list_f1 = []
s_list_mcc = []

#load rankings
for f in rankings_files:
    print('{}').format(f)
    df_rankings = pd.read_csv( f, index_col='molid' )
    df_rankings.index = df_rankings.index.map(str)
    # set up dictionaries to store F1 and MCC for each target
    temp_dict_f1 = {}
    temp_dict_mcc = {}

    for targ in df_rankings.columns:
        if df_rankings[targ].count() < 300:
            print("target {} missing significant portion of scored cpds, skipping metric eval".format(targ))
            temp_dict_f1[targ] = np.nan 
            temp_dict_mcc[targ] = np.nan
        else:
            df_rankings[targ].replace( 'informer', -1000.0, inplace=True )
            s_labels = df_binary[targ].rename('labels')
            s_rankings = df_rankings[targ].astype(float).rename('scores')
            df_temp = pd.concat( [s_labels, s_rankings], axis=1, sort=False )
            # for F1 and MCC evaluations do not consider negative informers as false positives
            df_temp = df_temp[ ~((df_temp['scores'] == -1000.0) & (df_temp['labels'] == False)) ]

            # predict compounds scoring in top 0.10 quantile predicted as active; all others inactive
            df_temp['binary_predictions'] = df_temp['scores'] < df_temp['scores'].quantile(0.10)
            labels_arr = df_temp['labels'].values
            predictions_arr = df_temp['binary_predictions'].values

            f1 = sk.metrics.f1_score( labels_arr, predictions_arr )
            mcc = sk.metrics.matthews_corrcoef( labels_arr, predictions_arr )
            temp_dict_f1[targ] = f1
            temp_dict_mcc[targ] = mcc

    s_f1 = pd.Series( temp_dict_f1 ).rename(f)
    s_list_f1.append(s_f1)
    s_mcc = pd.Series( temp_dict_mcc ).rename(f)
    s_list_mcc.append(s_mcc)

df_f1 = pd.concat( s_list_f1, axis=1 )
df_mcc = pd.concat( s_list_mcc, axis=1 )

df_f1.columns = [ 'BC_l', 'BC_s', 'BC_w', 'BF_l', 'BF_s', 'BF_w', 'RS', 'AS', 'CS' ]
df_mcc.columns = [ 'BC_l', 'BC_s', 'BC_w', 'BF_l', 'BF_s', 'BF_w', 'RS', 'AS', 'CS' ]
new_order = [ 'BC_s', 'BC_l', 'BC_w', 'BF_s', 'BF_l', 'BF_w', 'RS', 'CS', 'AS' ]

df_f1 = df_f1[ new_order ]
df_mcc = df_mcc[ new_order ]

df_f1.to_csv('../output_pkis1loto/metrics/pkis1loto_eval_F1_10perc_thresh_v1.2.csv', index_label='target')
df_mcc.to_csv('../output_pkis1loto/metrics/pkis1loto_eval_MCC_10perc_thresh_v1.2.csv', index_label='target')

