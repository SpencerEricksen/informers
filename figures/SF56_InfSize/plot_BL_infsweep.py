#!/home/ssericksen/anaconda2/bin/python2.7

import sys
import glob
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
mpl.style.use('ggplot')
import seaborn as sb

try:
  datafile1 = sys.argv[1]
  datafile2 = sys.argv[2]
  plotfile = sys.argv[3]
except:
  print('')
  print('usage: ./plot_violin_vXX.py  data_csv_rocauc data_csv_nef10  plotfilename')
  print('')
  exit(1)

df1 = pd.read_csv( datafile1, index_col=0 )
df2 = pd.read_csv( datafile2, index_col=0 )

cols = range(1,49)
df1.columns = cols
df2.columns = cols

fig, (ax1, ax2) = plt.subplots( ncols=2, sharey=True )

orientation = 'h'
sb.boxplot( data=df1, palette="Set3", notch=False, linewidth=1.0, fliersize=0.4, orient=orientation, ax=ax1)
sb.boxplot( data=df2, palette="Set3", notch=False, linewidth=1.0, fliersize=0.4, orient=orientation, ax=ax2)

ax1.set_xlabel('ROCAUC')
ax1.set_ylabel('N informers')
ax1.tick_params(axis='y', labelsize=8)
ax1.set_xlim( 0.0, 1.0 )

ax2.set_xlabel('NEF10')
ax2.set_xlim( 0.4, 1.0 )

fig.set_size_inches(8,6)
plt.savefig( plotfile, bbox_inches='tight', dpi=600 )

