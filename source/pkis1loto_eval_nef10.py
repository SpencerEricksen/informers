#!/home/ssericksen/anaconda2/bin/python2.7

import numpy as np
import pandas as pd
import informer_functions as inf
import glob

rankings_files = glob.glob('../output_pkis1loto/rankings/bl/ranked_df_16_*.csv')
rankings_files = rankings_files + glob.glob('../output_pkis1loto/rankings/hz/*.csv')
rankings_files = rankings_files + glob.glob('../output_pkis1loto/rankings/cp/*.csv')
# alphabetize the file list so we can re-arrange later
rankings_files.sort()

activity_matrix_file = '../data/pkis1.csv'

df_continuous = inf.get_continuous( activity_matrix_file )

df_binary = inf.get_binary( df_continuous )

s_list = []
#load rankings
for f in rankings_files:
  print('{}').format(f)
  df_rankings = pd.read_csv( f, index_col='molid' )
  df_rankings.index = df_rankings.index.map(str)
  # pandas doesn't like the types in the problematic 'LRRK2' target column so
  # will handle this target column by target column later in loop
  #df_rankings[df_rankings == 'informer' ] = -10.0
  temp_dict = {}
  for targ in df_rankings.columns:
    if df_rankings[targ].count() < 300:
        nef10 = 0.500
        #nef10 = np.nan
        temp_dict[targ] = nef10
    else:
        df_rankings[targ][ df_rankings[targ] == 'informer' ] = -100.0
        s_labels = df_binary[targ].rename('labels')
        s_rankings = df_rankings[targ].astype(float).rename('scores')
        df_temp = pd.concat( [s_labels, s_rankings], axis=1, sort=False )
        #filter df_temp to include only informers that were active for ROCAUC
        df_temp = df_temp[ ~((df_temp['scores'] == -100.000) & (df_temp['labels'] == False)) ]
        labels_arr = df_temp['labels'].values
        scores_arr = -1.00 * df_temp['scores'].astype(float).values

        nef10 = inf.normalized_enrichment_factor_single_hz( labels_arr, scores_arr, 0.10 )
        temp_dict[targ] = nef10 

  s = pd.Series( temp_dict ).rename(f)
  s_list.append(s)

df = pd.concat( s_list, axis=1 )

df.columns = [ 'BC_l', 'BC_s', 'BC_w', 'BF_l', 'BF_s', 'BF_w', 'RS', 'AS', 'CS' ]
new_order = [ 'BC_s', 'BC_l', 'BC_w', 'BF_s', 'BF_l', 'BF_w', 'RS', 'CS', 'AS' ]

df = df[ new_order ]

df.to_csv('../output_pkis1loto/metrics/pkis1loto_eval_NEF10.csv', index_label='target')
