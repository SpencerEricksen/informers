#!/home/ssericksen/anaconda2/bin/python2.7

import numpy as np
import pandas as pd
import sys

try:
    matrix = sys.argv[1]
except:
    print('')
    print('specify a matrix (pkis1 or pkis2)')
    print('')
    exit


df_rank_pknb = pd.read_csv('../../output_newtargs/'+matrix+'_pknb_model_rankings_v1.2.csv', index_col='molid' )
df_rank_pknb.index = df_rank_pknb.index.map(str)

df_rank_bglf4 = pd.read_csv('../../output_newtargs/'+matrix+'_bglf4_model_rankings_v1.2.csv', index_col='molid' )
df_rank_bglf4.index = df_rank_bglf4.index.map(str)

df_rank_rop18 = pd.read_csv('../../output_newtargs/'+matrix+'_rop18_model_rankings_v1.2.csv', index_col='molid' )
df_rank_rop18.index = df_rank_rop18.index.map(str)

df_act = pd.read_csv('../../data/data_newtargs_'+matrix+'cpds.csv', index_col='molid' )
df_act.index = df_act.index.map(str)

w_list = []

for col in ['b_cw', 'b_bw', 'RS', 'CS', 'AS']:
    w = df_rank_pknb[col][ df_rank_pknb[col] == 'informer' ]
    w.rename( col, inplace=True)
    w_list.append(w)

df = pd.concat( w_list, axis=1 )
df['pknb'] = df_act['pknb'].loc[df.index]
df['bglf4'] = df_act['bglf4'].loc[df.index]
df['rop18'] = df_act['rop18'].loc[df.index]

new_order = ['pknb', 'bglf4', 'rop18', 'b_cw', 'b_bw', 'RS', 'CS', 'AS' ]
df = df[new_order]

df.index.name = 'PubChem CID'
#df.round(3)
df.to_csv('informer_activities_'+matrix+'_S1table.csv', float_format='%.2f' )

