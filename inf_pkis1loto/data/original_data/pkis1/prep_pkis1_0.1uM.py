#!/home/ssericksen/anaconda2/bin/python2.7

import pandas as pd
import numpy as np


df = pd.read_csv('PKIS_screening_data.csv')

# extract lines using 1 uM threshold
df = df[ df['ASSAY'].str.contains(" 0.1 uM") ]

# make new column with unique target ID--concatenate first word in 'ASSAY' and the 'ASSAY_CHEMBL_ID'
df['targid'] = df['ASSAY'].apply(lambda x: x.split()[0]) + "_" + df['ASSAY_CHEMBL_ID']

# now build a list of series--one series for each targid
s_list = []

for t in df['targid'].unique():
    s = df[['CHEMBL_ID','VALUE']][df['targid'] == t ]
    s.set_index('CHEMBL_ID', inplace=True)
    s.rename( columns={"VALUE":t}, inplace=True )
    s_list.append( s )

'''
for t in df['targid'].unique():
    s = df[['SMILES','VALUE']][df['targid'] == t ]
    s.set_index('SMILES', inplace=True)
    s.rename( columns={"VALUE":t}, inplace=True )
    s_list.append( s )
'''

# merge the series into a dataframe (columns are targids, index is CHEMBL_ID for mols)
df2 = pd.concat( s_list, axis=1 )

# 3 mols were duplicates (both SMILES and CHEMBL_ID)
# keep only the first occurrence:
#df2 = df2[ ~df2.index.duplicated() ]

# better yet, use the average for the two occurrences:
df3 = df2.groupby('CHEMBL_ID').mean()

df3.to_csv('pkis1_0.1uM_raw.csv', index_label='CHEMBL_ID' )

