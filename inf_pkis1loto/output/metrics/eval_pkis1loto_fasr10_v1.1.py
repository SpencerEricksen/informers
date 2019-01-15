#!/home/ssericksen/anaconda2/bin/python2.7

import numpy as np
import pandas as pd
import glob


def get_continuous( activity_matrix_file ):
    '''read in activity data for matrix (Huikun's space delimited pkis1.csv),
       return a dataframe with cpd molid indices and target columns

       note: the data file is transposed'''
    df = pd.read_csv( activity_matrix_file, delimiter=" ", index_col=0 ).T
    df.index = df.index.map(str)
    return df

def get_binary( df ):
    ''' input continous activity data frame, return binary activity dataframe '''
    df_binary = df[ df > (df.mean(axis=0) + 2*df.std(axis=0)) ].notnull()
    return df_binary



# load the rankings matrix (PKIS1) for each method
rankings_files = glob.glob('../rankings/bl/ranked_df_16_*.csv')
rankings_files = rankings_files + glob.glob('../rankings/hz/*.csv')
rankings_files = rankings_files + glob.glob('../rankings/cp/*.csv')
rankings_files.sort()

activity_matrix_file = '../../data/pkis1.csv'

df_continuous = get_continuous( activity_matrix_file )

df_binary = get_binary( df_continuous )

# for fasr need scaffold data on PKIS1
df_scaffolds = pd.read_csv('../../data/compounds/scaffolds/pkis1_gen_scaffids.csv', index_col='molid')
df_scaffolds.index = df_scaffolds.index.map(str)

s_list = []
#load rankings
for f in rankings_files:
  print('{}').format(f)
  df_rankings = pd.read_csv( f, index_col='molid' )
  df_rankings.index = df_rankings.index.map(str)
  # pandas doesn't like the types in the problematic 'LRRK2' target column so
  # will handle this target column by target column later in loop
  #df_rankings[df_rankings == 'informer' ] = -10.0
  temp_dict = {}
  for targ in df_rankings.columns:
    try:
        df_rankings[targ][ df_rankings[targ] == 'informer' ] = -100.0
        s_labels = df_binary[targ].rename('labels')
        s_rankings = df_rankings[targ].astype(float).rename('scores')
        df_temp = pd.concat( [s_labels, s_rankings, df_scaffolds], axis=1, sort=False )

        # get scaffids in informer hits--we might get credit for these even if method fails
        informer_act_scaffids = list( df_temp.sort_values('scores').head(16)['scaffid'][ df_temp['labels'] == True ].unique() )

        #filter df_temp to include only informers that were active in ranked list
        df_temp = df_temp[ ~((df_temp['scores'] == -100.000) & (df_temp['labels'] == False)) ]

        # so with truncated dataset (with inactive informers removed), how many cpds in 10% of dataset?
        n_10perc = int(round(len(df_temp)/10.0))

        # get scaffid ids for all actives 
        all_act_scaffids = list( df_temp.sort_values('scores')['scaffid'][ df_temp['labels'] == True ].unique() )

        # get scaffids for actives in top 10%
        if df_temp['scores'].count() >= 300:
            top_act_scaffids = list( df_temp.sort_values('scores').head( n_10perc)['scaffid'][ df_temp['labels'] == True ].unique() )
        else:
            # can only take credit for "active" informers since ranking is random if evaluation failed
            top_act_scaffids = informer_act_scaffids

        # fraction active scaffolds recovered in top 10% (includes informers)
        fasr10 = len(top_act_scaffids) * 1.000 / len(all_act_scaffids)

        temp_dict[targ] = fasr10
        
    except:
        #temp_dict[targ] = np.nan
        temp_dict[targ] = 0.000
        continue

  s = pd.Series( temp_dict ).rename(f)
  s_list.append(s)

df = pd.concat( s_list, axis=1 )

df.columns = [ 'BC_l', 'BC_s', 'BC_w', 'BF_l', 'BF_s', 'BF_w', 'RS', 'AS', 'CS' ]
new_order = [ 'BC_s', 'BC_l', 'BC_w', 'BF_s', 'BF_l', 'BF_w', 'RS', 'CS', 'AS' ]

df = df[ new_order ]

df.to_csv('pkis1loto_eval_FASR10.csv', index_label='target')
