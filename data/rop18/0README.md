import pandas as pd
import numpy as np

df1 = pd.read_csv('../data_newtargs_pkis1cpds.csv', index_col='molid')

df2 = pd.read_csv('./data_rop18_pkis1cpds.csv', index_col='molid')

df3 = pd.concat( [df1, df2], axis=1, sort=False )

df3.to_csv('../data_newtargs3_pkis1cpds.csv', index_label='molid')
