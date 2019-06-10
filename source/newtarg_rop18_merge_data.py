import pandas as pd
import numpy as np

df = pd.read_excel( "../data/rop18/rop18_predictions_20190524.xlsx", worksheet='ROP18_Predictions')

# make new 'molid' column with no CID prefix for PKIS1 cpds
df.loc[ df['Group'] == 1, 'molid'] = df.loc[ df['Group'] == 1, 'CID'].apply( lambda x: str(int(x[3:])))
df.loc[ df['Group'] == 2, 'molid'] = df.loc[ df['Group'] == 2, 'CID']

# make 'molid' the index
df.set_index( 'molid', inplace=True )

# rename columns for HZ's and CP's methods
df.rename( columns={ 'HK_sd':'CS', 'HK_km':'AS', 'CP':'RS' }, inplace=True )

# spotted case issue with "informer"/"Informer"
df.replace( "Informer", "informer", inplace=True )

# flip sign on values for CP's "RS" method
df.loc[ df['RS']=='informer', 'RS' ] = 1000.0
df['RS'] = df['RS'] * -1.000
df.loc[ df['RS'] == -1000.0, 'RS'] = "informer"

# dump rop18 rankings for PKIS1 and PKIS2 cpd sets
df[['b_cs','b_cl','b_cw','b_bs','b_bl','b_bw','RS','CS','AS']].loc[ df['Group']==1].to_csv( '../output_newtargs/pkis1_rop18_model_rankings.csv', index_label='molid' )
df[['b_cs','b_cl','b_cw','b_bs','b_bl','b_bw','RS','CS','AS']].loc[ df['Group']==2].to_csv( '../output_newtargs/pkis2_rop18_model_rankings.csv', index_label='molid' )

