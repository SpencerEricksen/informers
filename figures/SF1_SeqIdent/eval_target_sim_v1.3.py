#!/home/ssericksen/anaconda2/bin/python2.7

import pandas as pd
import numpy as np
import sys

try:
    sheet = sys.argv[1] # 'pkis1_set' or pkis2_set'
    outcsv = sys.argv[2]
except:
    print ""
    print "script.py   sheetname   outputcsv            "
    print "            'pkis1_set' ./data/newtarg_sim_pkis1.csv"
    print ""
    exit

excel_file = "./data/set1_set2_identities_matrix.xlsx"

df0 = pd.read_excel( excel_file, sheet_name=sheet )
df0.index = df0.columns

# remove sequence fragments
if 'sp|Q9UI2' in df0.index:
    df0.drop( columns=['sp|Q9UI2'], index=['sp|Q9UI2'], inplace=True )

if 'sp|P5496' in df0.index:
    df0.drop( columns=['sp|P5496'], index=['sp|P5496'], inplace=True )

# bglf4 =  'AIE89081'
# pknb = 'sp|P9WI8'

dataset = sheet.split('_')[0]
df_12 = df0.max(axis=1).sort_values()
df_12.rename( 'sim12', inplace=True )
df_21 = df0.max(axis=0).sort_values()
df_21.rename( 'sim21', inplace=True )

df1 = pd.concat( [df_12, df_21], axis=1 )
df1.round(2).to_csv( outcsv, index_label='target' )

