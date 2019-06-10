
import pandas as pd

for metric in ['ROCAUC', 'NEF10', 'FASR10']:
    df_4 = pd.read_csv('pkis1loto_eval_'+metric+'_v1.2.csv', index_col='target')
    df_6 = pd.read_csv('ecfp6/pkis1loto_eval_'+metric+'_ecfp6_v1.2.csv', index_col='target')
    df = pd.concat( [ df_4.median(), df_6.median() ], axis=1 )
    df.columns = ['ecfp4', 'ecfp6']
    print metric
    print df
    print('')
    print('')


ROCAUC
         ecfp4     ecfp6
BC_s  0.598705  0.586244
BC_l  0.601407  0.584969
BC_w  0.669021  0.554183
BF_s  0.825847  0.807985
BF_l  0.816942  0.807374
BF_w  0.809964  0.792199
RS    0.925259  0.925259
CS    0.832218  0.832218
AS    0.877621  0.877621


NEF10
         ecfp4     ecfp6
BC_s  0.590643  0.555556
BC_l  0.585317  0.555556
BC_w  0.603175  0.525279
BF_s  0.736146  0.736842
BF_l  0.742063  0.742063
BF_w  0.722222  0.733980
RS    0.814815  0.814815
CS    0.795322  0.795322
AS    0.847082  0.847082


FASR10
         ecfp4     ecfp6
BC_s  0.263158  0.176471
BC_l  0.264912  0.176471
BC_w  0.294118  0.125000
BF_s  0.500000  0.500000
BF_l  0.527864  0.500000
BF_w  0.500000  0.500000
RS    0.720779  0.720779
CS    0.641711  0.641711
AS    0.750000  0.750000



for metric in ['F1', 'MCC']:
    df_4 = pd.read_csv('pkis1loto_eval_'+metric+'_estimated_thresh_v1.2.csv', index_col='target')
    df_6 = pd.read_csv('ecfp6/pkis1loto_eval_'+metric+'_estimated_thresh_ecfp6_v1.2.csv', index_col='target')
    df = pd.concat( [ df_4.median(), df_6.median() ], axis=1 )
    df.columns = ['ecfp4', 'ecfp6']
    print metric
    print df
    print('')
    print('')

     
F1
         ecfp4     ecfp6
BC_s  0.219807  0.127717
BC_l  0.213371  0.127717
BC_w  0.208333  0.082483
BF_s  0.468627  0.466667
BF_l  0.495614  0.465891
BF_w  0.442810  0.443152
RS    0.411011  0.411011
CS    0.551724  0.551724
AS    0.616516  0.616516


MCC
         ecfp4     ecfp6
BC_s  0.180706  0.085192
BC_l  0.166572  0.085192
BC_w  0.162309  0.026204
BF_s  0.443153  0.440661
BF_l  0.462040  0.442119
BF_w  0.415248  0.418025
RS    0.440731  0.440731
CS    0.526667  0.526667
AS    0.594085  0.594085


