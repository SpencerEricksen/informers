#!/home/ssericksen/anaconda2/bin/python2.7

import pandas as pd
import numpy as np


df1 = pd.read_csv('bglf4_pkis1', sep=" ")
df1.set_index('fid', inplace=True)
df1.columns = ['BGLF4']
df1.index.rename('molid', inplace=True)
df1.index = df1.index.map(str)

df2 = pd.read_csv('new_pkis1_informers_CP.csv', header=None)
df2.set_index(0, inplace=True)
df2.index.rename('molid', inplace=True)
df2.columns = ['BGLF4']
df2.index = df2.index.map(str)

df3 = pd.concat( [df1, df2], axis=0 )

df3['BGLF4'][ df3['BGLF4'] == 'Informer' ] = -10.0
df4 = df3[~df3.index.duplicated(keep='last')]

df4['BGLF4'][ df4['BGLF4'] == -10 ] = 'informer'
df4.to_csv( 'ching_pei_bglf4_pkis1.csv', index_label='molid' )

