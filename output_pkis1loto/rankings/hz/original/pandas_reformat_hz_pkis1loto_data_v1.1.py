#!/home/ssericksen/anaconda2/bin/python2.7

import pandas as pd
import numpy as np
import sys

input_data = sys.argv[1]

df = pd.read_csv( input_data, index_col=0, delimiter=" ")

basename = input_data.split('.')[0]

df.to_csv( basename+'.csv', index_label='molid' )

