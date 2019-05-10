# Evaluate VS metrics using rankings for baselines and non-baselines.

Run metrics calculations over all targets for Huikun's (CS, AS), Ching-pei's (RS), and my 6 baselines (Bx) using ECFP4.
```
./pkis1loto_eval_rocaucs.py
./pkis1loto_eval_nef10.py
./pkis1loto_eval_fasr10.py
```

> Note: for metrics evaluations, if no informers are returned as "active" (based on inferred active threshold), baseline ranking methods 's' and 'l' will fail to rank compounds. In such cases, the target is given a minimal metric score (ROCAUC=0.5, NEF10=0.5, or FASR10=0.0)


Re-evaluate models on PKIS1 LOTO for baseline models using ECFP6 (not significantly different).
```
./pkis1loto_eval_rocaucs_ecfp6.py
./pkis1loto_eval_nef10_ecfp6.py
./pkis1loto_eval_fasr10_ecfp6.py
```

Run F1-score and MCC metrics

> Note: For these metrics we will need to convert cpd ranking by models into to binary classifications. To do this, we estimated the number of actives (N) for each target based on outcomes from informers and then assigned 'active' labels to top scoring N molecules.
```
./pkis1loto_eval_F1_MCC_estimate_thresh.py
```
> Note: Also looked at using the "ground truth" thresholds (zscore >= +2) for each target for active compound predictions. This ensures that the predictions and real labels have same ratio of positives/negatives.
```
./pkis1loto_eval_F1_MCC_real_thresh.py
```

Run paired 2-sided Student's t-test.
> Note: These are probably not valid based on assumptions.
```
./pkis1loto_ttest_pvals_methods.py ../output_pkis1loto/metrics/pkis1loto_eval_ROCAUC.csv
./pkis1loto_ttest_pvals_methods.py ../output_pkis1loto/metrics/pkis1loto_eval_NEF10.csv
./pkis1loto_ttest_pvals_methods.py ../output_pkis1loto/metrics/pkis1loto_eval_FASR10.csv
```

Run non-parametric paired t-tests (Wilcoxon signed-rank test)
> Note: This is a valid test--though not much different than outcomes from parametric test above.
```
./pkis1loto_nonparam_ttest_pvals_methods.py ../output_pkis1loto/metrics/pkis1loto_eval_ROCAUC.csv
./pkis1loto_nonparam_ttest_pvals_methods.py ../output_pkis1loto/metrics/pkis1loto_eval_NEF10.csv
./pkis1loto_nonparam_ttest_pvals_methods.py ../output_pkis1loto/metrics/pkis1loto_eval_FASR10.csv 
```
