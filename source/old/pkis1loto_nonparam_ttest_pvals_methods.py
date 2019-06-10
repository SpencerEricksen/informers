#!/home/ssericksen/anaconda2/bin/python2.7

import sys
import pandas as pd
import numpy as np
from scipy import stats 

incsv = sys.argv[1]
#df = pd.read_csv('./rocauc/pkis1loto_rocauc_metrics_20181106.csv', index_col='target')
df = pd.read_csv( incsv, index_col='target' )

cols = list(df.columns)

p_array = np.empty( (len(cols), len(cols)) )
p_array[:] = np.nan

# run 2-sided pairwise t-test to obtain p-values
for i in range(len(cols)):
    for j in range(i):
        # remove rows with missing values for pairwise comparison (conservative)
        temp = df[ [ cols[i], cols[j] ] ].dropna()
        a = temp[ cols[i] ].values
        b = temp[ cols[j] ].values
        #s, p = stats.ttest_rel(a,b)
        #s, p = stats.ranksums(a,b)
        s, p = stats.wilcoxon(a,b)
        print('{} vs {}, tstat:{:.6f}, pval:{:.9f}, N:{}').format( cols[i], cols[j], s, p, len(temp) )
        p_array[i,j] = p
print('')

# print stats
print('model, mean, median, stdevp')
for c in cols:
    print('{}, {:.3f}, {:.3f}, {:.3f}').format( c, df[c].mean(), df[c].median(), df[c].std() )
print('*'*40)

print(',{}').format( ",".join(cols) )
dim = np.shape(p_array)
for i in range(dim[0]):
    print('{},{}').format( cols[i], ",".join( [ str(x) for x in list(p_array[i]) ] ) )

