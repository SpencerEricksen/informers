import numpy as np
import pandas as pd
import informer_functions as inf
import sys

try:
    matrix = sys.argv[1] # 1 or 2
    targ = sys.argv[2]   # pknb, bglf4, or rop18
except:
    print('')
    print(' python script.py   matrix   targ')
    print('')
    print('                    1 or 2   pknb, bglf4, or rop18')
    print('')
    exit

rankings_file = '../output_newtargs/pkis'+matrix+'_'+targ+'_model_rankings.csv'
activity_matrix_file = '../data/data_newtargs_pkis'+matrix+'cpds.csv'

df_continuous = pd.read_csv( activity_matrix_file, index_col='molid')

df_binary = inf.get_binary( df_continuous )
df_binary.index = df_binary.index.map(str)

df_rankings = pd.read_csv( rankings_file, index_col='molid' )
df_rankings.index = df_rankings.index.map(str)
df_rankings.replace('informer', -1000.0, inplace=True)

print('model,inf_hits,hits_recovered,tot_hits,NEF10')
for model in df_rankings.columns:
    if df_rankings[model].count() < 300:
        s_labels = df_binary[targ].rename('labels')
        s_rankings = df_rankings[model].astype(float).rename('scores')
        df_temp = pd.concat( [s_labels, s_rankings], axis=1, sort=False )
        inf_hits = df_temp[ df_temp['scores'] == -1000.0 ]['labels'].sum()
        tot_hits = df_temp['labels'].sum()
        hits_recovered = np.nan
        nef10 = np.nan
    else:
        s_labels = df_binary[targ].rename('labels')
        s_rankings = df_rankings[model].astype(float).rename('scores')
        df_temp = pd.concat( [s_labels, s_rankings], axis=1, sort=False )
        # do not count inactive informers as false positives
        df_temp = df_temp[ ~((df_temp['scores'] == -1000.0) & (df_temp['labels'] == False)) ]
        df_temp = df_temp.dropna( how='any')
        labels_arr = df_temp['labels'].values
        scores_arr = -1.00 * df_temp['scores'].astype(float).values
        tot_hits = labels_arr.sum()
        inf_hits = df_temp[ df_temp['scores'] == -1000.0 ]['labels'].sum()
        nef10 = inf.normalized_enrichment_factor_single_hz( labels_arr, scores_arr, 0.10 )
        # so with truncated dataset (with negative informers removed), how many cpds in 10% of dataset?
        N = int( round( len(df_temp) * 0.10 ) )
        hits_recovered = df_temp.sort_values('scores').head( N )['labels'].sum()
    print('{},{},{},{},{}').format( model, inf_hits, hits_recovered, tot_hits, nef10 )
