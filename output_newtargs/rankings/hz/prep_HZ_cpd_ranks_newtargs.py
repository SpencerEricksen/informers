#!/home/ssericksen/anaconda2/bin/python2.7

import pandas as pd
import numpy as np
import sys

# read in the pkis matrix number (as string)
matrix = str(sys.argv[1])  # 1 or 2
model = str(sys.argv[2])   # AS or CS
targ = str(sys.argv[3])    # pknb or bglf4

df2 = pd.read_csv( './missing_inf/'+model+'_'+targ+'_pkis'+matrix+'_noinf.csv', index_col='molid')
df2.index = df2.index.map(str)

df1 = pd.read_csv('../../../data/pkis'+matrix+'_continuous_labels.csv', index_col='molid' )
df1.index = df1.index.map(str)

df3 = df2.reindex( df1.index )

df3[ df3.isnull() ] = 'informer'

df3.to_csv( model+'_'+targ+'_pkis'+matrix+'.csv', index_label='molid' )

'''
df_merge = pd.concat( [df_coding, df_as], axis=1, sort=False )
df_merge.index.name = 'molid'

# chop off that CID prefix if matrix is '2'
if matrix == '2':
    df_as.index = [ str(int(c[3:])) for c in df_as.index ]

df_merge.to_csv('df_'+targ+'_pkis'+matrix+'.csv', index_label='molid' )
'''

