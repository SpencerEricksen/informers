# baseline methods

run baseline methods on new targets
```
./newtargs_baseline_inf_sel_rank_v1.2.py 16 1
./newtargs_baseline_inf_sel_rank_v1.2.py 16 2
```

build individual method/matrix data rankings
```
./newtargs_build_rankings_v1.2.py 1
./newtargs_build_rankings_v1.2.py 2

```
merge all methods into data ranking for each target/matrix
```
./newtargs_merge_data_v1.2.py 1 bglf4
./newtargs_merge_data_v1.2.py 2 bglf4
./newtargs_merge_data_v1.2.py 1 pknb
./newtargs_merge_data_v1.2.py 2 pknb
```
Already extracted rop18 from complete ranking of all methods. Used sed to replace "Informers" to "informers" for some of the non-baseline methods for rop18 rankings.
```
./newtarg_rop18_merge_data_v1.2.py
```

run evaluation metrics for each metric-target combination
```
for m in 1 2; do for targ in pknb bglf4 rop18;  do echo "${m},${targ}"; ./newtargs_eval_F1_MCC_pkismedian_thresh_v1.2.py ${m} ${targ}; echo  ""; done; done > ../output_newtargs/results_newtargs_F1_MCC_pkismedian.csv 
for m in 1 2; do for targ in pknb bglf4 rop18;  do echo "${m},${targ}"; ./newtargs_eval_F1_MCC_10perc_thresh_v1.2.py ${m} ${targ}; echo  ""; done; done > ../output_newtargs/results_newtargs_F1_MCC_10perc.csv 
for m in 1 2; do for targ in pknb bglf4 rop18;  do echo "${m},${targ}"; ./newtargs_eval_rocauc_v1.2.py ${m} ${targ}; echo  ""; done; done > ../output_newtargs/results_newtargs_rocauc.csv
for m in 1 2; do for targ in pknb bglf4 rop18;  do echo "${m},${targ}"; ./newtargs_eval_nef10_v1.2.py ${m} ${targ}; echo  ""; done; done > ../output_newtargs/results_newtargs_nef10.csv
for m in 1 2; do for targ in pknb bglf4 rop18;  do echo "${m},${targ}"; ./newtargs_eval_fasr10_v1.2.py ${m} ${targ}; echo  ""; done; done > ../output_newtargs/results_newtargs_fasr10.csv
```

