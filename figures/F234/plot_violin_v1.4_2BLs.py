#!/home/ssericksen/anaconda2/bin/python2.7

import sys
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
mpl.style.use('ggplot')
import seaborn as sb

try:
  datafile = sys.argv[1]
  metric = sys.argv[2]
  plotfile = sys.argv[3]
  orientation = sys.argv[4]
except:
  print('')
  print('usage: ./plot_violin_vXX.py  datafilename  metric  plotfilename  orientation')
  print('')
  exit(1)


df = pd.read_csv( datafile, index_col='target' )
df.columns = [ 'b_cs', 'b_cl', 'b_cw', 'b_bs', 'b_bl', 'b_bw', 'RS', 'CS', 'AS' ]
df.columns = [ 'BC_s', 'BC_l', 'BC_w', 'BF_s', 'BF_l', 'BF_w', 'RS', 'CS', 'AS' ]

df = df[ ['BC_w', 'BF_w', 'RS', 'CS', 'AS'] ]

fig, ax = plt.subplots()

#ax.set_title('PKIS1 LOTO (N=224 targets)') #, fontsize=10)

ax = sb.violinplot( data=df, palette="Set3", inner='box', scale="count", bw=0.1, alpha=1.0,  cut=0, linewidth=1.5, orient=orientation, zorder=0 )
#sb.violinplot( data=df, palette="Set3",    inner='stick', scale="count", bw=0.1, alpha=0.5,  cut=0, linewidth=0.5, orient=orientation, ax=ax )
sb.swarmplot( data=df, color='k', size=2, alpha=0.25, ax=ax, orient=orientation )

if orientation == 'v':
    ax.set_xlabel('IBR model') #, fontsize=8)
    ax.set_ylabel(metric)
    ax.tick_params(axis='x', labelsize=8)

elif orientation == 'h':
    ax.set_xlabel(metric)
    ax.set_ylabel('IBR model')
    ax.tick_params(axis='y', labelsize=8)

#ax.grid(False)

# width * height
fig.set_size_inches( 8,5 )

#plt.show()
plt.savefig( plotfile, bbox_inches='tight', dpi=600 )

