#!/usr/bin/env python2

import pandas as pd
import numpy as np
import sys

matrix = sys.argv[1] # 1 or 2
targ = sys.argv[2]   # bglf4, pknb, or rop18

df1 = pd.read_csv('../output_newtargs/bl/BL_'+targ+'_pkis'+matrix+'.csv', index_col='molid' )
df2 = pd.read_csv('../output_newtargs/cp/RS_'+targ+'_pkis'+matrix+'.csv', index_col='molid' ) 
df3 = pd.read_csv('../output_newtargs/hz/CS_'+targ+'_pkis'+matrix+'.csv', index_col='molid' ) 
df4 = pd.read_csv('../output_newtargs/hz/AS_'+targ+'_pkis'+matrix+'.csv', index_col='molid' )

df5 = pd.concat( [df1, df2, df3, df4], axis=1, sort=False )


df5[ ['b_cs', 'b_cl', 'b_cw', 'b_bs', 'b_bl', 'b_bw', 'RS', 'CS', 'AS' ] ].to_csv('../output_newtargs/pkis'+matrix+'_'+targ+'_model_rankings.csv', index_label='molid')

