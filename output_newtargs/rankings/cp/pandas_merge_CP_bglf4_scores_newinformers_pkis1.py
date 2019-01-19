#!/home/ssericksen/anaconda2/bin/python2.7

import pandas as pd
import numpy as np

# load Ching-Pei's compound scores for BGLF4 with PKIS1
df1 = pd.read_csv('bglf4_pkis1', sep=" ")
df1.set_index('fid', inplace=True)
df1.columns = ['BGLF4']
df1.index.rename('molid', inplace=True)
df1.index = df1.index.map(str)

# load informer list as dataframe
df2 = pd.read_csv('new_pkis1_informers_CP.csv', header=None)
df2.set_index(0, inplace=True)
df2.index.rename('molid', inplace=True)
df2.columns = ['BGLF4']
df2.index = df2.index.map(str)

# merge dataframes
df3 = pd.concat( [df1, df2], axis=0 )
print("duplicated indices: {}").format( df3.duplicated().sum() )

# check duplicates for PKIS1 molid '11959682'
print( df3.loc['11959682'] )

