#!/home/ssericksen/anaconda2/bin/python2.7

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
import glob


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

def compute_roc_auc( labels_arr, scores_arr ):
    '''use an sklearn function to compute ROC AUC
        probably should add some other metrics to this'''
    if len(np.unique(labels_arr)) == 2:
        auc = roc_auc_score( labels_arr, scores_arr )
    else:
        auc = 'ND'
    return auc



rankings_files = glob.glob('../rankings/bl/ranked_df_16_*.csv')
rankings_files = rankings_files + glob.glob('../rankings/hz/*.csv')
rankings_files = rankings_files + glob.glob('../rankings/cp/*.csv')
# alphabetize the file list so we can re-arrange later
rankings_files.sort()

activity_matrix_file = '../../../data/pkis1.csv'

df_continuous = get_continuous( activity_matrix_file )

df_binary = get_binary( df_continuous )

s_list = []
#load rankings
for f in rankings_files:
  print('{}').format(f)
  df_rankings = pd.read_csv( f, index_col='molid' )
  df_rankings.index = df_rankings.index.map(str)
  # pandas doesn't like the types in the problematic 'LRRK2' target column so
  # let's handle this target column by target column later in loop
  #df_rankings[df_rankings == 'informer' ] = -10.0
  temp_dict = {}
  for targ in df_rankings.columns:
    if df_rankings[targ].count() < 300:
        rocauc = 0.500
        #rocauc = np.nan
        temp_dict[targ] = rocauc
    else:
        df_rankings[targ][ df_rankings[targ] == 'informer' ] = -100.0
        s_labels = df_binary[targ].rename('labels')
        s_rankings = df_rankings[targ].astype(float).rename('scores')
        df_temp = pd.concat( [s_labels, s_rankings], axis=1, sort=False )
        #filter df_temp to include only informers that were active for ROCAUC
        df_temp = df_temp[ ~((df_temp['scores'] == -100.000) & (df_temp['labels'] == False)) ]
        labels_arr = df_temp['labels'].values
        scores_arr = -1.00 * df_temp['scores'].astype(float).values

        rocauc = compute_roc_auc( labels_arr, scores_arr )
        temp_dict[targ] = rocauc

  s = pd.Series( temp_dict ).rename(f)
  s_list.append(s)

df = pd.concat( s_list, axis=1 )

# , ./bl/ranked_df_16_clst_se_thresh.csv, ./bl/ranked_df_16_prom_le_thresh.csv, ./bl/ranked_df_16_prom_we_thresh.csv,
#   ./bl/ranked_df_16_clst_we_thresh.csv, ./bl/ranked_df_16_clst_le_thresh.csv, ./bl/ranked_df_16_prom_se_thresh.csv,
#   ./hz/coding_cv_ranking.csv, ./hz/adaptive_cv_ranking.csv, ./cp/df_CP_rocauc_pkis1loto.csv

df.columns = [ 'BC_l', 'BC_s', 'BC_w', 'BF_l', 'BF_s', 'BF_w', 'RS', 'AS', 'CS' ]

new_order = [ 'BC_s', 'BC_l', 'BC_w', 'BF_s', 'BF_l', 'BF_w', 'RS', 'CS', 'AS' ]

df = df[ new_order ]

df.to_csv('pkis1loto_eval_ROCAUC.csv', index_label='target')
