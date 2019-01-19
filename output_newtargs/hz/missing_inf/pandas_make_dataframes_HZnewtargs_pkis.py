#!/home/ssericksen/anaconda2/bin/python2.7

import pandas as pd
import numpy as np
import sys

# read in the pkis matrix number (as string)
matrix = str(sys.argv[1])
targ = str(sys.argv[2])

df = pd.read_csv( '../original/'+targ+'_pkis'+matrix+'.csv', sep=" " )

df_coding = df[[ c for c in df.columns if 'coding_' in c ]]
df_coding.set_index( 'coding_CID', inplace=True )
df_coding.index = df_coding.index.map(str)
df_coding.index.name = 'molid'
df_coding.columns = ['CS']
df_coding.to_csv( 'CS_'+targ+'_pkis'+matrix+'_noinf.csv', index_label='molid' )

df_as = df[[ c for c in df.columns if 'as_' in c ]]
df_as.set_index('as_CID', inplace=True )
df_as.index = df_as.index.map(str)
df_as.index.name = 'molid'
df_as.columns = ['AS']
df_as.to_csv( 'AS_'+targ+'_pkis'+matrix+'_noinf.csv', index_label='molid' )

'''
df_merge = pd.concat( [df_coding, df_as], axis=1, sort=False )
df_merge.index.name = 'molid'

# chop off that CID prefix if matrix is '2'
if matrix == '2':
    df_as.index = [ str(int(c[3:])) for c in df_as.index ]

df_merge.to_csv('df_'+targ+'_pkis'+matrix+'.csv', index_label='molid' )
'''

