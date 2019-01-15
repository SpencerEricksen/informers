#!/home/ssericksen/anaconda2/bin/python2.7

import sys
import pandas as pd
import numpy as np
import baseline_inf_sel_rank_methods_v22 as bl

matrix = str(sys.argv[1])    # 1 or 2

df_act_mat = bl.get_continuous( './data/pkis'+matrix+'.csv' )

df_binary = bl.get_binary( df_act_mat )
fract_act = df_binary.sum() / len(df_binary)
fract_act.index.name = 'target'
fract_act.rename('active_fraction', inplace=True)
fract_act.to_csv( './data/active_fraction_pkis'+matrix+'.csv', header=True)

