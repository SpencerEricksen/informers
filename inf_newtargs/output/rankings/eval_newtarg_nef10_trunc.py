#!/home/ssericksen/anaconda2/bin/python2.7

import numpy as np
import pandas as pd
import vsmetrics
import baseline_inf_sel_rank_methods_v22 as baseline
import sys

try:
    matrix = sys.argv[1] # 1 or 2
    targ = sys.argv[2]   # pknb or bglf4
except:
    print('')
    print(' eval_rocauc_newtarg.py  matrix   targ')
    print('')
    print('                         1 or 2   pknb or bglf4')
    print('')
    exit

rankings_file = 'pkis'+matrix+'_'+targ+'_model_rankings.csv'
activity_matrix_file = '/home/ssericksen/matrix_informer/models/data/newtargs/data_newtargs_pkis'+matrix+'cpds.csv'
#molid,pknb,bglf4
#CID006539592,0.0000,0.0000

df_continuous = pd.read_csv( activity_matrix_file, index_col='molid')

df_binary = baseline.get_binary( df_continuous )
df_binary.index = df_binary.index.map(str)


df_rankings = pd.read_csv( rankings_file, index_col='molid' )
df_rankings.index = df_rankings.index.map(str)

print('model,inf_hits,hits_recovered,tot_hits,NEF10')
for model in df_rankings.columns:
    if df_rankings[model].count() < 300:
        df_rankings[model][ df_rankings[model] == 'informer' ] = -100.0
        s_labels = df_binary[targ].rename('labels')
        s_rankings = df_rankings[model].astype(float).rename('scores')
        df_temp = pd.concat( [s_labels, s_rankings], axis=1, sort=False )
        inf_hits = df_temp[ df_temp['scores'] == -100.000 ]['labels'].sum()
        tot_hits = df_temp['labels'].sum()
        #rocauc = np.nan
        print('{},{},{},{},fail').format( model, inf_hits, inf_hits, tot_hits )
    else:
        df_rankings[model][ df_rankings[model] == 'informer' ] = -100.0
        s_labels = df_binary[targ].rename('labels')
        s_rankings = df_rankings[model].astype(float).rename('scores')
        df_temp = pd.concat( [s_labels, s_rankings], axis=1, sort=False )
        #filter df_temp to include only informers that were active for ROCAUC
        df_temp = df_temp[ ~((df_temp['scores'] == -100.000) & (df_temp['labels'] == False)) ]
        df_temp = df_temp.dropna( how='any')
        labels_arr = df_temp['labels'].values
        scores_arr = -1.00 * df_temp['scores'].astype(float).values
        tot_hits = labels_arr.sum()
        inf_hits = df_temp[ df_temp['scores'] == -100.000 ]['labels'].sum()
        nef10 = vsmetrics.normalized_enrichment_factor_single_hz( labels_arr, scores_arr, 0.10 )
        # so with truncated dataset (with negative informers removed), how many cpds in 10% of dataset?
        n_10perc = int(round(len(df_temp)/10.0))
        hits_recovered = df_temp.sort_values('scores').head( n_10perc )['labels'].sum()
        #rocauc = vsmetrics.compute_roc_auc( labels_arr, scores_arr )
        #print model, rocauc, nef10
        print('{},{},{},{},{}').format( model, inf_hits, hits_recovered, tot_hits, nef10 )
        #print df_temp.sort_values('scores', ascending=True)
        #print ""

