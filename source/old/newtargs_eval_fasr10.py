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

# get the data files for scores, labels, and scaffolds
rankings_file = 'pkis'+matrix+'_'+targ+'_model_rankings.csv'
activity_matrix_file = '../../../data/data_newtargs_pkis'+matrix+'cpds.csv'
df_scaffolds = pd.read_csv('../../../data/compounds/scaffolds/pkis'+matrix+'_gen_scaffids.csv', index_col='molid')
df_scaffolds.index = df_scaffolds.index.map(str)

# get the cpd labels
df_continuous = pd.read_csv( activity_matrix_file, index_col='molid')
df_binary = inf.get_binary( df_continuous )
df_binary.index = df_binary.index.map(str)

# get the cpd scores
df_rankings = pd.read_csv( rankings_file, index_col='molid' )
df_rankings.index = df_rankings.index.map(str)

print('model,recognized_act_scaffold_IDs,all_active_scaffold_IDs,fasr10')
for model in df_rankings.columns:
    df_rankings[model][ df_rankings[model] == 'informer' ] = -100.0
    s_labels = df_binary[targ].rename('labels')
    s_rankings = df_rankings[model].astype(float).rename('scores')
    df_temp = pd.concat( [s_labels, s_rankings, df_scaffolds], axis=1, sort=False )

    # get scaffids in informer hits--we might get credit for these even if method fails
    informer_act_scaffids = list( df_temp.sort_values('scores').head(16)['scaffid'][ df_temp['labels'] == True ].unique() )

    # get scaffid ids for all actives
    all_act_scaffids = list( df_temp.sort_values('scores')['scaffid'][ df_temp['labels'] == True ].unique() )

    #filter df_temp to include only informers that were active for ROCAUC
    df_temp = df_temp[ ~((df_temp['scores'] == -100.000) & (df_temp['labels'] == False)) ]
    df_temp = df_temp.dropna( how='any')

    # so with truncated dataset (with negative informers removed), how many cpds in 10% of dataset?
    n_10perc = int(round(len(df_temp)/10.0))

    # get scaffids for actives in top 10%
    if df_temp['scores'].count() >= 300:
        top_act_scaffids = list( df_temp.sort_values('scores').head( n_10perc)['scaffid'][ df_temp['labels'] == True ].unique() )
    else:
        # can only take credit for "active" informers since ranking is random if evaluation failed
        top_act_scaffids = informer_act_scaffids

    # fraction active scaffolds recognized in top 10% (includes informers)
    fasr10 = len(top_act_scaffids) * 1.000 / len(all_act_scaffids)

    print('{},{},{},{}').format( model, " ".join([ str(i) for i in top_act_scaffids]), " ".join([str(i) for i in all_act_scaffids]), str(fasr10) )


