#!/home/ssericksen/anaconda2/bin/python2.7

import numpy as np
import pandas as pd
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



rankings_files = glob.glob('../rankings/bl/ranked_df_16_*.csv')
rankings_files = rankings_files + glob.glob('../rankings/hz/*.csv')
rankings_files = rankings_files + glob.glob('../rankings/cp/*.csv')
rankings_files.sort()

activity_matrix_file = '../../data/pkis1.csv'

df_continuous = get_continuous( activity_matrix_file )

df_binary = get_binary( df_continuous )

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

        nef10 = normalized_enrichment_factor_single_hz( labels_arr, scores_arr, 0.10 )
        temp_dict[targ] = nef10 

  s = pd.Series( temp_dict ).rename(f)
  s_list.append(s)

df = pd.concat( s_list, axis=1 )

df.columns = [ 'BC_l', 'BC_s', 'BC_w', 'BF_l', 'BF_s', 'BF_w', 'RS', 'AS', 'CS' ]
new_order = [ 'BC_s', 'BC_l', 'BC_w', 'BF_s', 'BF_l', 'BF_w', 'RS', 'CS', 'AS' ]

df = df[ new_order ]

df.to_csv('pkis1loto_eval_NEF10.csv', index_label='target')
