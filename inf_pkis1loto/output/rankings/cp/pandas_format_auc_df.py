#!/home/ssericksen/anaconda2/bin/python2.7

import pandas as pd
import numpy as np

df1 = pd.read_csv('./original/pkis1_auc_loo.csv', header=None)
df1.index = df1.index + 1
df1 = df1 * -1.000

df2 = pd.read_csv( './original/pkis1_auc_informer.csv', header=None)

df3 = pd.read_csv( '../hz/adaptive_cv_ranking.csv', index_col='molid')
idx_list = df3.index
col_list = df3.columns

#df3 = pd.read_csv( 'hz_molid_index_pkis1.csv', index_col='molid')
#idx_list = df3.index

for c in df1.columns:
    df1[c].loc[ df2[c] ] = 'informer'

# re-label the indices and columns using Huikun's order (from his dataframe, df3)
df1.index = idx_list
df1.columns = col_list


df1.to_csv('df_CP_rocauc_pkis1loto.csv', index_label='molid')


