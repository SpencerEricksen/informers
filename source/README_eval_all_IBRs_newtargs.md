### process rankings data from raw baseline IRB output on new targets
`./newtargs_build_rankings.py 1`\
`./newtargs_build_rankings.py 2`

### merge the individual model rankings for each matrix/target combination
`./newtargs_merge_data.py 1 bglf4` \
`./newtargs_merge_data.py 2 bglf4` \
`./newtargs_merge_data.py 1 pknb` \
`./newtargs_merge_data.py 2 pknb`

### evaluate metrics on merged IBR data for new targets
`for m in 1 2; do for targ in pknb bglf4; do echo "${m},${targ}"; ./newtargs_eval_nef10_trunc.py ${m} ${targ}; echo  ""; done; done > results_newtargs_nef10.csv` \
\
`for m in 1 2; do for targ in pknb bglf4; do echo "${m},${targ}"; ./newtargs_eval_fasr10.py ${m} ${targ}; echo  ""; done; done > results_newtargs_fasr10.csv` \
\
`for m in 1 2; do for targ in pknb bglf4; do echo "${m},${targ}"; ./newtargs_eval_rocauc_trunc.py ${m} ${targ}; echo  ""; done; done > results_newtargs_rocauc.csv`

