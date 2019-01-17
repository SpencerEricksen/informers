#!/home/ssericksen/anaconda2/bin/python2.7

import numpy as np
import pandas as pd

def get_continuous( activity_matrix_file ):
    '''read in activity data for matrix (Huikun's space delimited pkis1.csv),
       return a dataframe with cpd molid indices and target columns

       note: the data file is transposed'''
    df = pd.read_csv( activity_matrix_file, delimiter=" ", index_col=0 ).T
    df.index = df.index.map(str)
    return df

def get_binary( df ):
    ''' input continuous activity data frame, return binary activity dataframe '''
    df_binary = df[ df > (df.mean(axis=0) + 2*df.std(axis=0)) ].notnull()
    return df_binary

def process_informer_list( inf_molid_list, inf_activity_list, thresh=None ):
    '''read in full informer list and return list of [active] informers in
       order of descending activity, threshold may be input as an option'''
    # get activity sorted informer list for only the "active" informers
    inf_tups = zip( inf_molid_list, inf_activity_list )
    inf_tups = sorted( inf_tups, reverse=True, key=lambda x: x[1] )
    if thresh != None:
        inf_tups = [ i for i in inf_tups if i[1] >= float(thresh) ]
    new_inf_list = [ x[0] for x in inf_tups ]
    return new_inf_list 

###################################################################################

import sys
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
mpl.style.use('ggplot')

try:
    colorby = sys.argv[1]
    colormap = sys.argv[2]
    rankings_data_csv = sys.argv[3]
    rankthresh = float(sys.argv[4])
except:
    print('')
    print('')
    print('usage: ./plot...py   activity  cividis  ../ranked_clst_we_thresh.csv  clst_we_civ_act     1.00          ')
    print('')
    print('                     colorby   colormap   rankingdata                 graph_image_rootname  rankthresh    ')
    print('')
    print('')
    print('       colorby: activity or ranking')
    print('')
    print('       colormap: cividis, magma, viridis, plasma, inferno, ...')
    print('')
    print('')
    exit


activity_data_csv = './data/pkis1.csv'
fps_file = './data/pkis1_uniq_ecfp4_1024.fps'

# read in compound feature data
df_fps = pd.read_csv( fps_file, index_col=0 )
df_fps.index = df_fps.index.map(str)
# split up bitstrings into integer bit vectors (feature vectors)
df_fps['features'] = df_fps['fingerprint'].apply( lambda x: np.array(list(map(int,list(x)))))
tup = tuple( df_fps['features'].values )
# stack feature vectors into matrix (rows:examples, columns:features)
X = np.vstack( tup )
# get mol identifiers "molids"
molids = list(df_fps.index)


# read in data matrix
df_act_mat = get_continuous( activity_data_csv )
# re-order the activity dataframe indices to match with feature data matrix rows
df_act_mat = df_act_mat.reindex( index=molids )

# get binary data matrix
df_binary = get_binary( df_act_mat )
#re-order the binary dataframe indices to match with feature data matrix rows
df_binary = df_binary.reindex( index=molids )

# read in rankings matrix
df_rankings = pd.read_csv( rankings_data_csv, index_col='molid' )
df_rankings.index = df_rankings.index.map(str)
df_rankings = df_rankings.reindex( index=molids )

# boolean mask for informers
df_informer_bool = df_rankings == 'informer'

# normalize scores
df_rankings_perc = df_rankings.copy()
df_rankings_perc[ df_rankings_perc == 'informer' ] = np.nan
df_rankings_perc = pd.DataFrame( df_rankings_perc, dtype='float')
df_rankings_perc = df_rankings_perc.rank( pct=True )
df_rankings_perc.index = df_rankings_perc.index.map(str)

#run PCA
pca = PCA(n_components=10)
X_r = pca.fit(X).transform(X)
# this helps spacing better in molecular viewer
X_r = 10.0 * X_r

targlist = ['EGFR', 'LOK']

fig = plt.figure()
for i in range(len(targlist)):
        # get ordinal index values for informers and noninformers
        targ = targlist[i]
        idx_inf  = np.flatnonzero(  df_informer_bool[targ] )
        idx_ninf = np.flatnonzero( ~df_informer_bool[targ] )
        idx_ninf_rankthresh = np.flatnonzero( df_rankings_perc[targ] <= rankthresh )
        # obtain list of informer cpds via different informer selection methods
        inf_molids = list( df_rankings[targ].iloc[ idx_inf ].index )
        # plot projection of fps onto 2 PC dimensions
        ax = fig.add_subplot(1,2,i+1, projection='3d')
        #axes[i].grid( b=True )
        ax.set_title('PCA:PKIS1 cpds, target:{}'.format( targ),  fontsize=8 )
        # get informer and noninformer PCA coords
        X_r_i = X_r[ idx_inf ]
        X_r_n = X_r[ idx_ninf_rankthresh ]
        
        # get informer/noninformer edge colors (red=active, black=inactive)
        df_binary[targ].iloc[ idx_inf ]
        inf_colors = [ 'red' if i == True else "None" for i in df_binary[targ].iloc[ idx_inf ] ]
        ninf_colors = [ 'red' if i == True else "None" for i in df_binary[targ] ]
        # normalize plot point colors
        normalize = mpl.colors.Normalize(vmin=0.0, vmax=1)
        if colorby == 'activity':
            scat_i = ax.scatter( X_r_i[:,0], X_r_i[:,1], X_r_i[:,2], c=df_act_mat[targ].iloc[ idx_inf ].values, cmap=plt.get_cmap(colormap), norm=normalize, alpha=0.7, lw=0.1, s=10, edgecolor='none' )
            scat_n = ax.scatter( X_r_n[:,0], X_r_n[:,1], X_r_n[:,2], c=df_act_mat[targ].iloc[ idx_ninf_rankthresh ].values,   cmap=plt.get_cmap( colormap ), norm=normalize, alpha=0.7,  lw=0.5, s=10, edgecolor=list( np.array(ninf_colors)[ idx_ninf_rankthresh ]) ) #edgecolor='black' )

        elif colorby == 'ranking':
            scat_i = ax.scatter( X_r_i[:,0], X_r_i[:,1], X_r_i[:,2], c=inf_colors, alpha=1.0, lw=0, s=10, edgecolor=inf_colors )
            scat_n = ax.scatter( X_r_n[:,0], X_r_n[:,1], X_r_n[:,2], c=df_rankings_perc[targ].iloc[ idx_ninf_rankthresh ].values, cmap=plt.get_cmap( colormap ), norm=normalize, alpha=1.0,  lw=0, s=10, edgecolor=list( np.array(ninf_colors)[ idx_ninf_rankthresh ]) )
        ax.set_xlim(-40.0, 40.0)
        ax.set_ylim(-40.0, 40.0)
        ax.set_zlim(-40.0, 40.0)
        ax.set_xlabel('PCA-0', fontsize=6)
        ax.set_ylabel('PCA-1', fontsize=6)
        ax.set_zlabel('PCA-2', fontsize=6)
        ax.tick_params(axis = 'both', which = 'major', labelsize = 4)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes( [0.85, 0.15, 0.05, 0.7] )
fig.colorbar( scat_n, cax=cbar_ax)
fig.set_size_inches(10,4)
plt.savefig( 'sfig3_pca_morgan_pkis1_egfr_lok.png', dpi=600 )
plt.close()

