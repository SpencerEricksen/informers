# baseline methods

run baseline methods on new targets (BGLF4, PknB, and ROP18) using 16 informers on both pkis1 and pkis2 cpd sets
```
python newtargs_baseline_inf_sel_rank.py 16 1
python newtargs_baseline_inf_sel_rank.py 16 2
```

build individual method/matrix data rankings for the baseline methods
```
python newtargs_build_rankings.py 1
python newtargs_build_rankings.py 2

```
merge cpd rankings from all methods (baselines and non-baselines) into complete ranking data set for each target/matrix
```
python newtargs_merge_data.py 1 bglf4
python newtargs_merge_data.py 2 bglf4
python newtargs_merge_data.py 1 pknb
python newtargs_merge_data.py 2 pknb
python newtargs_merge_data.py 1 rop18
python newtargs_merge_data.py 2 rop18
```

run evaluation metrics for each metric-target combination
```
for m in 1 2; do for targ in pknb bglf4 rop18;  do echo "${m},${targ}"; python newtargs_eval_F1_MCC_pkismedian_thresh.py ${m} ${targ}; echo  ""; done; done > ../output_newtargs/results_newtargs_F1_MCC_pkismedian.csv 
for m in 1 2; do for targ in pknb bglf4 rop18;  do echo "${m},${targ}"; python newtargs_eval_rocauc.py ${m} ${targ}; echo  ""; done; done > ../output_newtargs/results_newtargs_rocauc.csv
for m in 1 2; do for targ in pknb bglf4 rop18;  do echo "${m},${targ}"; python newtargs_eval_nef10.py ${m} ${targ}; echo  ""; done; done > ../output_newtargs/results_newtargs_nef10.csv
for m in 1 2; do for targ in pknb bglf4 rop18;  do echo "${m},${targ}"; python newtargs_eval_fasr10.py ${m} ${targ}; echo  ""; done; done > ../output_newtargs/results_newtargs_fasr10.csv
```

