#!/home/ssericksen/anaconda2/bin/python2.7

import pandas as pd
import numpy as np

# load in the predictions plus ground truth assay data for rop18
df = pd.read_excel( "rop18_predictions_20190524.xlsx", worksheet='ROP18_Predictions')

# create molid column with PKIS1 cpds having no CID prefix
df.loc[ df['Group'] == 1, 'molid'] = df.loc[ df['Group'] == 1, 'CID'].apply( lambda x: str(int(x[3:])))
df.loc[ df['Group'] == 2, 'molid'] = df.loc[ df['Group'] == 2, 'CID']

# make 'molid' the index
df.set_index( 'molid', inplace=True )

# normalize inhibition data
df['rop18'] = np.nan
df['rop18'] = df["ROP18"] / 100.0

# dump to files
df.loc[ df['Group'] == 1, ['rop18']].to_csv('data_rop18_pkis1cpds.csv', index_label='molid')
df.loc[ df['Group'] == 2, ['rop18']].to_csv('data_rop18_pkis2cpds.csv', index_label='molid')

