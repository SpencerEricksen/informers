#!/home/ssericksen/anaconda2/bin/python2.7

import pandas as pd
import numpy as np


df1 = pd.read_csv( '../../../data/pkis2_continuous_labels.csv', index_col='molid' )
df1.index = df1.index.map(str)

df2 = pd.read_csv( 'chingpei_predictions_pknb.csv' )
df2['molid'] = df2['CID']
#df2['molid'] = df2['CID'].apply( lambda x: str( int( x[3:] )))
df2.set_index('molid', inplace=True)

df3 = df2['ChingPei']
df3.index = df3.index.map(str)

df4 = df3.loc[ df1.index ]
# in future, passing list-likes to .loc or [] with any missing label will raise
# a KeyError, so use .reindex() as alternative
#df4 = df3.reindex( df1.index )

df4[ df4 == "Informer" ] = "informer"

df4.to_csv('ching_pei_pknb_pkis1.csv', index_label='molid' )

