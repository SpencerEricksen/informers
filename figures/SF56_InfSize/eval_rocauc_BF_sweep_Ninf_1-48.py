#!/home/ssericksen/anaconda2/bin/python2.7

import sys
import numpy as np
import pandas as pd
import vsmetrics
import baseline_inf_sel_rank_methods_v23 as baseline
import glob

n_inf_set = range(1,49)

activity_matrix_file = './data/pkis1.csv'

df_continuous = baseline.get_continuous( activity_matrix_file )

df_binary = baseline.get_binary( df_continuous )

s_list = []

#load output from each pkis1 loto (each informer set size)
for n in n_inf_set: #rankings_files:
  f = 'ranked_df_'+str(n)+'_prom_we_thresh_v23.csv'
  print('{}').format(f)
  df_rankings = pd.read_csv( f, index_col='molid' )
  df_rankings.index = df_rankings.index.map(str)
  df_rankings[df_rankings == 'informer' ] = -100.0
  temp_dict = {}
  for targ in df_rankings.columns:
    if df_rankings[targ].count() < 300:
        rocauc = 0.500
        temp_dict[targ] = rocauc
    else:
        s_labels = df_binary[targ].rename('labels')
        s_rankings = df_rankings[targ].astype(float).rename('scores')
        df_temp = pd.concat( [s_labels, s_rankings], axis=1, sort=False )
        # remove informers that were inactive so they are not counted as false positives
        df_temp = df_temp[ ~( (df_temp['labels'] == False) & (df_temp['scores'] == -100.0) ) ]
        labels_arr = df_temp['labels'].values
        scores_arr = -1.00 * df_temp['scores'].astype(float).values
        rocauc = vsmetrics.compute_roc_auc( labels_arr, scores_arr )
        temp_dict[targ] = rocauc

  s = pd.Series( temp_dict ).rename(f)
  s_list.append(s)

df = pd.concat( s_list, axis=1 )
df.columns = n_inf_set
df.to_csv('./data/baseline_pkis1_rocauc_evaluations_BF_sweep_Ninf_1-48.csv')
