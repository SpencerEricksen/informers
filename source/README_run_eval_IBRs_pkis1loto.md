# Evaluate VS metrics using rankings for baselines and non-baselines.

Run metrics calculations over all targets for Huikun's (CS, AS), Ching-pei's (RS), and my 6 baselines (Bx) using ECFP4.
```
./pkis1loto_eval_rocaucs.py
./pkis1loto_eval_nef10.py
./pkis1loto_eval_fasr10.py
./pkis1loto_eval_F1_MCC_pkis1median_thresh.py
```

> Note: for metrics evaluations, if no informers are returned as "active" (based on inferred active threshold), baseline ranking methods 's' and 'l' will use single most active compound for nearest-neighbors rankings.

Re-evaluate models on PKIS1 LOTO for baseline models using ECFP6 (not significantly different).
```
./pkis1loto_eval_rocaucs_ecfp6.py
./pkis1loto_eval_nef10_ecfp6.py
./pkis1loto_eval_fasr10_ecfp6.py
./pkis1loto_eval_F1_MCC_pkis1median_thresh_ecfp6.py
```

Run non-parametric paired t-tests (Wilcoxon signed-rank test)
```
./pkis1loto_nonparam_ttest_pvals_methods.py ../output_pkis1loto/metrics/pkis1loto_eval_ROCAUC.csv
./pkis1loto_nonparam_ttest_pvals_methods.py ../output_pkis1loto/metrics/pkis1loto_eval_NEF10.csv
./pkis1loto_nonparam_ttest_pvals_methods.py ../output_pkis1loto/metrics/pkis1loto_eval_FASR10.csv
./pkis1loto_nonparam_ttest_pvals_methods.py ../output_pkis1loto/metrics/pkis1loto_eval_F1_pkis1median_thresh.csv
./pkis1loto_nonparam_ttest_pvals_methods.py ../output_pkis1loto/metrics/pkis1loto_eval_MCC_pkis1median_thresh.csv
```

