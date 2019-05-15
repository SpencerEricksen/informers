#!/home/ssericksen/anaconda2/bin/python2.7

import pandas as pd
import numpy as np
import sys

matrix = str(sys.argv[1])

df = pd.read_csv('../output_newtargs/bl/rankings_baseline_pkis'+matrix+'.csv', index_col='molid' )

df_pknb = df[[ c for c in df.columns if 'pknb' in c ]]
df_pknb.columns = ['b_cs', 'b_cl', 'b_cw', 'b_bs', 'b_bl', 'b_bw' ]
df_pknb.to_csv('../output_newtargs/bl/BL_pknb_pkis'+matrix+'.csv', index_label='molid' )

df_bglf4 = df[[ c for c in df.columns if 'bglf4' in c ]]
df_bglf4.columns = ['b_cs', 'b_cl', 'b_cw', 'b_bs', 'b_bl', 'b_bw' ]
df_bglf4.to_csv('../output_newtargs/bl/BL_bglf4_pkis'+matrix+'.csv', index_label='molid' )
