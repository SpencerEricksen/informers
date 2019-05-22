#!/home/ssericksen/anaconda2/bin/python2.7

import pandas as pd
import numpy as np
import sys

matrix = str(sys.argv[1])

df = pd.read_csv('../output_newtargs/bl/rankings_3_baseline_pkis'+matrix+'_v1.2.csv', index_col='molid' )

df_rop18 = df[[ c for c in df.columns if 'rop18' in c ]]
df_rop18.columns = ['b_cs', 'b_cl', 'b_cw', 'b_bs', 'b_bl', 'b_bw' ]
df_rop18.to_csv('../output_newtargs/bl/BL_rop18_pkis'+matrix+'_v1.2.csv', index_label='molid' )

