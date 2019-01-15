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
    data_column = 2 
    bin_width = 0.005
    xmin = 0.0
    xmax = 0.14
    ymax = 55.0 
    # active fractions for each target (pknb or bglf4) with each matrix (1 or 2)
    pknb_1 = 0.0219     # 8/366
    bglf4_1 = 0.0301    # 11/366
    pknb_2 = 0.0169     # 7/415
    bglf4_2 = 0.0241    # 10/415
except:
    print " "
    print " plot_hist.py   numbers1.csv   numbers2.csv   title"
    print "                (datafile1)    (datafile2)    title"
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

# build histograms
bin_positions = np.arange( x_limits_data[0], x_limits_data[1] + bin_width, bin_width)
n, bins, patches = ax1.hist( num_list1, bins=bin_positions, facecolor='skyblue', edgecolor='k', linewidth=0.5, alpha=0.75 )
n, bins, patches = ax2.hist( num_list2, bins=bin_positions, facecolor='skyblue', edgecolor='k', linewidth=0.5, alpha=0.75 )

fig.suptitle('Fraction of Active Compounds on Targets', fontsize=14)
ax1.set_title( "PKIS1, m={}, median={:.3f}, sigma={:.3f}".format( N1, median1, sigma1 ), fontsize=8 )
ax2.set_title( "PKIS2, m={}, median={:.3f}, sigma={:.3f}".format( N2, median2, sigma2 ), fontsize=8 )
ax1.set_ylabel('number of targets')
ax1.set_xlabel('fraction active')
ax2.set_xlabel('fraction active')
ax1.plot( [pknb_1], [1], marker="D", linestyle="", color='tomato' ) 
ax1.plot( [bglf4_1], [1], marker="D", linestyle="", color='slateblue' )
ax2.plot( [pknb_2], [1], marker="D", linestyle="", color='tomato' )
ax2.plot( [bglf4_2], [1], marker="D", linestyle="", color='slateblue' )
ax1.set_xlim( x_limits_data[0], x_limits_data[1]) #20, 50)
ax2.set_xlim( x_limits_data[0], x_limits_data[1])
ax1.set_ylim( y_limits_data[0], y_limits_data[1])

fig.set_size_inches(7,4)
plt.savefig('hist_pkis12_actfract.png', bbox_inches='tight')
