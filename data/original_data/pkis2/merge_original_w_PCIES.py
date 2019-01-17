#!/home/ssericksen/anaconda2/bin/python2.7

import pandas as pd
import numpy as np

df1 = pd.read_csv('520888970260023795.txt', sep='\t', header=None)
df1.set_index(0, inplace=True)

df2 = pd.read_csv('pkis2_smiles_rdkit_desalt.can', sep=" ", header=None)
df2.set_index(0, inplace=True)

d1  = df1[1].to_dict()
d2 = df2[1].to_dict()

# make new dictionary with smiles as keys and [ Compound, CID ] as the values
d3 = {}
for k in d1:
    print d2[k], d1[k]
    d3[k] = [ d2[k], d1[k] ]


df3 = pd.DataFrame.from_dict( d3, orient='index' )

# make unique Compound names the index

# make new column for smiles
df3[2] = df3.index

# make the 'Compound' column the new index
df3.set_index(0, inplace=True)

# name  the columns and index
df3.columns = ['CID', 'smiles_nosalt']
df3.index.name = 'Compound'

# write it out
df3.to_csv('pkis2_Compound_CID_DesaltCanSmi.csv', index_label='Compound' )

'''
will have to add in some duplicates by hand:


In [113]: for c in list(df1[ df1.index.duplicated() ].index):
     ...:     print df1[1].loc[c], df2[1].loc[c]

0
NCCNS(=O)(=O)c1ccc(-c2ccnc3[nH]ccc23)cc1    23648890.0
NCCNS(=O)(=O)c1ccc(-c2ccnc3[nH]ccc23)cc1    23648890.0
Name: 1, dtype: float64 0
NCCNS(=O)(=O)c1ccc(-c2ccnc3[nH]ccc23)cc1    UNC10225023A
NCCNS(=O)(=O)c1ccc(-c2ccnc3[nH]ccc23)cc1    UNC10225023B
Name: 1, dtype: object
0
CC(C)(CN)c1nc(-c2ccc(Cl)c(O)c2)c(-c2ccncc2)[nH]1    22185475.0
CC(C)(CN)c1nc(-c2ccc(Cl)c(O)c2)c(-c2ccncc2)[nH]1    22185475.0
Name: 1, dtype: float64 0
CC(C)(CN)c1nc(-c2ccc(Cl)c(O)c2)c(-c2ccncc2)[nH]1    UNC10225206A
CC(C)(CN)c1nc(-c2ccc(Cl)c(O)c2)c(-c2ccncc2)[nH]1    UNC10225206B
Name: 1, dtype: object
0
COc1ccc(-c2nc3n(c2-c2ccncc2)CCC3)cc1    469038.0
COc1ccc(-c2nc3n(c2-c2ccncc2)CCC3)cc1    469038.0
Name: 1, dtype: float64 0
COc1ccc(-c2nc3n(c2-c2ccncc2)CCC3)cc1    UNC10225334A
COc1ccc(-c2nc3n(c2-c2ccncc2)CCC3)cc1    UNC10225334B
Name: 1, dtype: object
0
Fc1ccc(-c2ncn(C3CCNCC3)c2-c2ccnc(Nc3ccccc3)n2)cc1    9888013.0
Fc1ccc(-c2ncn(C3CCNCC3)c2-c2ccnc(Nc3ccccc3)n2)cc1    9888013.0
Name: 1, dtype: float64 0
Fc1ccc(-c2ncn(C3CCNCC3)c2-c2ccnc(Nc3ccccc3)n2)cc1    UNC10225420A
Fc1ccc(-c2ncn(C3CCNCC3)c2-c2ccnc(Nc3ccccc3)n2)cc1    UNC10225420B
Name: 1, dtype: object
0
CSc1nn2c(=O)cc(N3CCOCC3)nc2n1Cc1cccc(C(F)(F)F)c1C    70685442.0
CSc1nn2c(=O)cc(N3CCOCC3)nc2n1Cc1cccc(C(F)(F)F)c1C    70685442.0
Name: 1, dtype: float64 0
CSc1nn2c(=O)cc(N3CCOCC3)nc2n1Cc1cccc(C(F)(F)F)c1C    UNC10225417A
CSc1nn2c(=O)cc(N3CCOCC3)nc2n1Cc1cccc(C(F)(F)F)c1C    UNC10243881A

'''

