import numpy as np
import pandas as pd
import informer_functions as inf
import glob

# load the rankings matrix (PKIS1) for each method
rankings_files = glob.glob('../output_pkis1loto/rankings/bl/ecfp6/ranked_df_16_*.csv')
rankings_files = rankings_files + glob.glob('../output_pkis1loto/rankings/hz/*.csv')
rankings_files = rankings_files + glob.glob('../output_pkis1loto/rankings/cp/*.csv')
# alphabetize the file list so we can re-arrange later
rankings_files.sort()

activity_matrix_file = '../data/pkis1.csv'
df_continuous = inf.get_continuous( activity_matrix_file )
df_binary = inf.get_binary( df_continuous )

# FASR requires scaffold data on compound set
df_scaffolds = pd.read_csv('../data/compounds/scaffolds/pkis1_gen_scaffids.csv', index_col='molid')
df_scaffolds.index = df_scaffolds.index.map(str)

s_list = []
#load rankings
for f in rankings_files:
    print( 'evaluating FASR10 on rankings from: {}'.format(f) )
    df_rankings = pd.read_csv( f, index_col='molid' )
    df_rankings.index = df_rankings.index.map(str)
    temp_dict = {}
    for targ in df_rankings.columns:
        try:
            df_rankings[targ].replace( 'informer', -1000.0, inplace=True )
            s_labels = df_binary[targ].rename('labels')
            s_rankings = df_rankings[targ].astype(float).rename('scores')
            df_temp = pd.concat( [s_labels, s_rankings, df_scaffolds], axis=1, sort=False )

            # get scaffids among active informers
            informer_act_scaffids = list( df_temp.sort_values('scores').head(16)['scaffid'][ df_temp['labels'] == True ].unique() )

            # for FASR10 evaluation count only the active informers and non-informers in top 10%
            # (do not consider inactive informers as false positives)
            df_temp = df_temp[ ~((df_temp['scores'] == -1000.0) & (df_temp['labels'] == False)) ]

            # with inactive informers removed, how many cpds in 10% of dataset?
            n_10perc = int(round(len(df_temp)/10.0))
            
            # get scaffid ids for all actives 
            all_act_scaffids = list( df_temp.sort_values('scores')['scaffid'][ df_temp['labels'] == True ].unique() )
            
            # get scaffids for actives in top 10%
            top_act_scaffids = list( df_temp.sort_values('scores').head( n_10perc)['scaffid'][ df_temp['labels'] == True ].unique() )
            '''
            if df_temp['scores'].count() >= 300:
                 top_act_scaffids = list( df_temp.sort_values('scores').head( n_10perc)['scaffid'][ df_temp['labels'] == True ].unique() )
            else:
                 # can only take credit for "active" informers since ranking is random if evaluation failed
                 top_act_scaffids = informer_act_scaffids
            '''

            # fraction active scaffolds recovered in top 10% (includes informers)
            fasr10 = len(top_act_scaffids) * 1.000 / len(all_act_scaffids)
            temp_dict[targ] = fasr10
        
        except:
            temp_dict[targ] = np.nan
            print( 'failed to evaluate: {}'.format(targ) )
            #temp_dict[targ] = 0.000
            continue

    s = pd.Series( temp_dict ).rename(f)
    s_list.append(s)

df = pd.concat( s_list, axis=1 )
df.columns = [ 'BC_l', 'BC_s', 'BC_w', 'BF_l', 'BF_s', 'BF_w', 'RS', 'AS', 'CS' ]
new_order = [ 'BC_s', 'BC_l', 'BC_w', 'BF_s', 'BF_l', 'BF_w', 'RS', 'CS', 'AS' ]
df = df[ new_order ]
df.to_csv('../output_pkis1loto/metrics/pkis1loto_eval_FASR10_ecfp6.csv', index_label='target')
