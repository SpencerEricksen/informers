#!/home/ssericksen/anaconda2/bin/python2.7

import numpy as np
import pandas as pd
import informer_functions as inf
import glob
import sys
import sklearn as sk

inf_selection = 'prom'

rankings_files = glob.glob('../output_pkis1loto/rankings/bl/ranked_df_16_*.csv')
rankings_files = rankings_files + glob.glob('../output_pkis1loto/rankings/hz/*.csv')
rankings_files = rankings_files + glob.glob('../output_pkis1loto/rankings/cp/*.csv')
# alphabetize the file list so we can re-arrange later
rankings_files.sort()

activity_matrix_file = '../data/pkis1.csv'

df_continuous = inf.get_continuous( activity_matrix_file )

df_binary = inf.get_binary( df_continuous )

# read in thresholds data (provided by Huikun--predicted 2sigma thresholds based on informer activities)
'''
if inf_selection == 'prom':
    df_thresh = pd.read_csv( '../data/thresholds_2sigma/prom_le_thresh.txt', \
            index_col=0, header=None, names=['Activity'], delimiter=" " )
elif inf_selection == 'clst':
    df_thresh = pd.read_csv( '../data/thresholds_2sigma/clst_le_thresh.txt', \
            index_col=0, header=None, names=['Activity'], delimiter=" " )
else:
    print('issue with 2sigma thresholds!')
    exit(1)
'''

# possible to back-calculate the expected active fraction?

s_list_f1 = []
s_list_mcc = []

#load rankings
for f in rankings_files:
  print('{}').format(f)
  df_rankings = pd.read_csv( f, index_col='molid' )
  df_rankings.index = df_rankings.index.map(str)
  # pandas doesn't like the types in the problematic 'LRRK2' target column so
  # will handle this target column by target column later in loop
  #df_rankings[df_rankings == 'informer' ] = -10.0
  temp_dict_f1 = {}
  temp_dict_mcc = {}

  for targ in df_rankings.columns:
    if df_rankings[targ].count() < 300:
        temp_dict_f1[targ] = 0.000 
        temp_dict_mcc[targ] = 0.000
    else:
        df_rankings[targ][ df_rankings[targ] == 'informer' ] = -100.0
        s_labels = df_binary[targ].rename('labels')
        s_rankings = df_rankings[targ].astype(float).rename('scores')
        df_temp = pd.concat( [s_labels, s_rankings], axis=1, sort=False )

        #filter df_temp to keep only informers that were active for evaluation
        df_temp = df_temp[ ~((df_temp['scores'] == -100.000) & (df_temp['labels'] == False)) ]

        # compounds scored in top nth quantile predicted as active; all others inactive
        #df_temp['predictions'] = df_temp['scores'] < df_temp['scores'].quantile(0.10)

        # binary prediction: assign 'active' label to N top ranking \ 
        # cpds where N is number of true actives on target
        N = df_temp['labels'].sum()
        df_temp['binary_predictions'] = df_temp['scores'] <= df_temp['scores'].sort_values()[N]

        # or probably better, binary predictions: assign 'active' label to N top ranking \
        # cpds where N is predicted number of true actives on target
        #N = df_continuous.loc[ df_continuous[targ] >= df_thresh.loc[targ][0], targ ].count()
        #df_temp['binary_predictions'] = df_temp['scores'] <= df_temp['scores'].sort_values()[N]

        labels_arr = df_temp['labels'].values
        predictions_arr = df_temp['binary_predictions'].values
        #scores_arr = -1.00 * df_temp['scores'].astype(float).values

        #nef10 = inf.normalized_enrichment_factor_single_hz( labels_arr, scores_arr, 0.10 )
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

df_f1.to_csv('../output_pkis1loto/metrics/pkis1loto_eval_F1_real_thresh.csv', index_label='target')
df_mcc.to_csv('../output_pkis1loto/metrics/pkis1loto_eval_MCC_real_thresh.csv', index_label='target')

