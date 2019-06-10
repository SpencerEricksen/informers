#!/home/ssericksen/anaconda2/bin/python2.7

# evaluate F1 and MCC metrics on new targets. Assume 10% hit fractions,
# and predict top 10% of cpds by score as the actives

import numpy as np
import pandas as pd
import informer_functions as inf
import sklearn as sk
import sys

try:
    matrix = sys.argv[1] # 1 or 2
    targ = sys.argv[2]   # pknb, bglf4, or rop18
except:
    print('')
    print(' eval_rocauc_newtarg.py  matrix   targ')
    print('')
    print('                         1 or 2   pknb, bglf4, or rop18')
    print('')
    exit

rankings_file = '../output_newtargs/pkis'+matrix+'_'+targ+'_model_rankings_v1.2.csv'
activity_matrix_file = '../data/data_newtargs_pkis'+matrix+'cpds.csv'

df_continuous = pd.read_csv( activity_matrix_file, index_col='molid')
df_binary = inf.get_binary( df_continuous )
df_binary.index = df_binary.index.map(str)

df_rankings = pd.read_csv( rankings_file, index_col='molid' )
df_rankings.index = df_rankings.index.map(str)

df_rankings.replace('informer', -1000.0, inplace=True)

print('model,inf_hits,hits_recovered,tot_hits,F1,MCC')
for model in df_rankings.columns:
    if df_rankings[model].count() < 300:
        print("model:{} and target:{} missing significant portion of scored cpds, skipping metric eval".format(model,targ))
        s_labels = df_binary[targ].rename('labels')
        s_rankings = df_rankings[model].astype(float).rename('scores')
        df_temp = pd.concat( [s_labels, s_rankings], axis=1, sort=False )
        inf_hits = df_temp[ df_temp['scores'] == -1000.0 ]['labels'].sum()
        tot_hits = df_temp['labels'].sum()
        hits_recovered = np.nan 
        f1 = np.nan
        mcc = np.nan
    else:
        s_labels = df_binary[targ].rename('labels')
        s_rankings = df_rankings[model].astype(float).rename('scores')
        df_temp = pd.concat( [s_labels, s_rankings], axis=1, sort=False )
        df_temp = df_temp[ ~((df_temp['scores'] == -1000.0) & (df_temp['labels'] == False)) ]
        df_temp = df_temp.dropna( how='any')
        # to convert scores to binary predictions use thresholds based on median or mean active fractions across pkis(1/2) targets
        #act_fract_mean = df_binary.sum().mean() / len(df_binary)
        act_fract_median = df_binary.sum().median() / len(df_binary)
        df_temp['binary_predictions'] = df_temp['scores'] <= df_temp['scores'].quantile( act_fract_median )
        predictions_arr = df_temp['binary_predictions'].values
        labels_arr = df_temp['labels'].values
        tot_hits = labels_arr.sum()
        inf_hits = df_temp[ df_temp['scores'] == -1000.0 ]['labels'].sum()
        f1 = sk.metrics.f1_score( labels_arr, predictions_arr )
        mcc = sk.metrics.matthews_corrcoef( labels_arr, predictions_arr )

        # so with truncated dataset (with negative informers removed), how many cpds in predicted active fraction of dataset?
        N        = int( round( len(df_temp) * act_fract_median ) )
        hits_recovered = df_temp.sort_values('scores').head( N )['labels'].sum()

    print('{},{},{},{},{},{}').format( model, inf_hits, hits_recovered, tot_hits, f1, mcc )
