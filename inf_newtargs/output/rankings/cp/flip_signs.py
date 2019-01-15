#!/home/ssericksen/anaconda2/bin/python2.7

import pandas as pd
import numpy as np
import sys

rs_scores = sys.argv[1] # 'RS_bglf4_pkis1.csv'

# read in old data file
df = pd.read_csv( rs_scores, index_col='molid' )

# cpds with 'informer' as score, convert to large float
df['RS'][ df['RS'] == 'informer' ] = 1000.0
df = df.astype(float)

#flip the sign for RS model scores
df = df * -1.000

# convert cpds with -1000.0 score back to informers
df['RS'][ df['RS'] == -1000.0 ] = "informer"

# dump
df.to_csv( rs_scores+'2', index_label='molid' )

