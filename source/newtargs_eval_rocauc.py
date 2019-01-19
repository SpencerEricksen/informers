#!/home/ssericksen/anaconda2/bin/python2.7

import numpy as np
import pandas as pd
import informer_functions as inf
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
activity_matrix_file = '../../../data/data_newtargs_pkis'+matrix+'cpds.csv'
#molid,pknb,bglf4
#CID006539592,0.0000,0.0000

df_continuous = pd.read_csv( activity_matrix_file, index_col='molid')

df_binary = inf.get_binary( df_continuous )
df_binary.index = df_binary.index.map(str)

df_rankings = pd.read_csv( rankings_file, index_col='molid' )
df_rankings.index = df_rankings.index.map(str)

print('model,ROCAUC')
for model in df_rankings.columns:
    if df_rankings[model].count() < 300:
        #rocauc = np.nan
        print model, 'fail'
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
        rocauc = inf.compute_roc_auc( labels_arr, scores_arr )
        print('{},{}').format( model, str(rocauc) )

