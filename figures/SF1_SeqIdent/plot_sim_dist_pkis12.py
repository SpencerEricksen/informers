#!/home/ssericksen/anaconda2/bin/python2.7

import numpy as np
## agg backend is used to create plot as a .png file
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys
#from matplotlib.ticker import MultipleLocator, FormatStrFormatter


try:
    numbers_file1  = sys.argv[1]
    numbers_file2  = sys.argv[2]
    data_column = 3
    bin_width = 2
    xmin = 0.0
    xmax = 100.0
    ymax = 24.0
    # list of minimal sequence identities to nearest neighbor (obtained from matrix)
    pknb_1 = 16.1
    bglf4_1 = 13.8
    pknb_2 = 16.1
    bglf4_2 = 14.2
except:
    print " "
    print " plot_hist.py   numbers1.csv   numbers2.csv "
    print "                (datafile1)    (datafile2)  "
    print ""

def getNumData(dfile, col_num):
    '''read in data file, get values from specified column'''
    with open(dfile,'r') as fh:
	data = fh.readlines()
    data_values_list = []
    for line in data[1:]:
        try:
            data_values_list.append( float(line.split(',')[ col_num - 1 ]) )
        except:
            continue
    return data_values_list

def getStats(num_list):
    '''calc stats based on the data values'''
    num_arr  = np.array(num_list)
    mean     = np.mean(num_arr)
    median   = np.median(num_arr)
    stdev    = np.std(num_arr)
    N        = len(num_arr)
    return mean, median, stdev, N

# main

# load in data file and specify column
# returns a list of data values
num_list1 = getNumData( numbers_file1, data_column )
num_list2 = getNumData( numbers_file2, data_column )

# calc stats for distribution data
mu1, median1, sigma1, N1 = getStats( num_list1 )
mu2, median2, sigma2, N2 = getStats( num_list2 )

x_limits_data = ( xmin, xmax )
y_limits_data = ( 0, ymax )

fig = plt.figure()
ax1, ax2 = fig.subplots(1, 2, sharey=True)
# build histogram
bin_positions = np.arange( x_limits_data[0], x_limits_data[1] + bin_width, bin_width)
n, bins, patches = ax1.hist( num_list1, bins=bin_positions, facecolor='skyblue', edgecolor='k', linewidth=0.5, alpha=0.75 )
n, bins, patches = ax2.hist( num_list2, bins=bin_positions, facecolor='skyblue', edgecolor='k', linewidth=0.5, alpha=0.75 )

fig.suptitle('Sequence Identities of Nearest Neighbor Targets', fontsize=14)
ax1.set_title( "PKIS1, median={:.1f}, sigma={:.1f}".format( median1, sigma1 ), fontsize=8 )
ax2.set_title( "PKIS2, median={:.1f}, sigma={:.1f}".format( median2, sigma2 ), fontsize=8 )
ax1.set_ylabel('number of targets')
ax1.set_xlabel('sequence identity (%)')
ax2.set_xlabel('sequence identity (%)')
ax1.plot( [pknb_1], [1], marker="D", linestyle="", color='tomato' ) 
ax1.plot( [bglf4_1], [1], marker="D", linestyle="", color='slateblue' )
ax2.plot( [pknb_2], [1], marker="D", linestyle="", color='tomato' )
ax2.plot( [bglf4_2], [1], marker="D", linestyle="", color='slateblue' )
ax1.set_xlim( x_limits_data[0], x_limits_data[1]) #20, 50)
ax2.set_xlim( x_limits_data[0], x_limits_data[1])
ax1.set_ylim( y_limits_data[0], y_limits_data[1])

fig.set_size_inches(7,4)
plt.savefig('hist_pkis1_pkis2_21_max_seq_ident.png', bbox_inches='tight')
