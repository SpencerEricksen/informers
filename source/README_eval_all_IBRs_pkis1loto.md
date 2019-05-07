# evaluate metrics using rankings for baselines and non-baselines

run metrics calculations over all targets (Huikun's, Ching-pei's, and my baselines using ECFP4
```
./pkis1loto_eval_rocaucs.py
./pkis1loto_eval_nef10.py
./pkis1loto_eval_fasr10.py
```

> Note: for metrics evaluations, if no informers are returned as "active" (based on inferred active threshold), baseline ranking methods 's' and 'l' will fail to rank compounds. In such cases, the target is given a minimal metric score (ROCAUC=0.5, NEF10=0.5, or FASR10=0.0)


re-evaluate models on PKIS1 LOTO, this time baseline models use ECFP6
```
./pkis1loto_eval_rocaucs_ecfp6.py
./pkis1loto_eval_nef10_ecfp6.py
./pkis1loto_eval_fasr10_ecfp6.py
```



run paired Student's t-test
```
./pkis1loto_ttest_pvals_methods.py ../output_pkis1loto/metrics/pkis1loto_eval_ROCAUC.csv
./pkis1loto_ttest_pvals_methods.py ../output_pkis1loto/metrics/pkis1loto_eval_NEF10.csv
./pkis1loto_ttest_pvals_methods.py ../output_pkis1loto/metrics/pkis1loto_eval_FASR10.csv
```

run non-parametric paired t-tests (Wilcoxon signed-rank test)
```
./pkis1loto_nonparam_ttest_pvals_methods.py ../output_pkis1loto/metrics/pkis1loto_eval_ROCAUC.csv
./pkis1loto_nonparam_ttest_pvals_methods.py ../output_pkis1loto/metrics/pkis1loto_eval_NEF10.csv
./pkis1loto_nonparam_ttest_pvals_methods.py ../output_pkis1loto/metrics/pkis1loto_eval_FASR10.csv 
```
