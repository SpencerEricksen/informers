#!/home/ssericksen/anaconda2/bin/python2.7

import pandas as pd
import numpy as np

# load my CSV that matched the "Compound" identifiers from original PKIS2 matrix to desalted rdkit cansmi and CIDs
# from PubChem Identifier Exchange Service
df1 = pd.read_csv('pkis2_download_journal.pone.0181585_s004_smiles_pandas_addedCIDs_shuffled.csv', index_col='Compound')

# the original spreadsheet with PFE-PKIS series of compounds edited with added underscores so 
# "PFE-PKIS 10" becomes PFE-PKIS_10 to avoid any possible issues with space delimitation
df2 = pd.read_csv('pkis2_download_journal.pone.0181585.s004.csv', index_col='Compound')

# merge together using the unique "Compound" identifier as the index
df3 = pd.concat( [df1, df2], axis=1 )
df3.to_csv('pkis2_download_journal.pone.0181585.s004_wCIDs_NoSaltSmis.csv', index_label='Compound')

