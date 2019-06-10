# PKIS1 LOTO (RDKit Morgan fingerprints with radius=2, 1024 bits)

run chemometric baselines:
```
python pkis1loto_baseline_inf_sel_rank.py 16 clst se
python pkis1loto_baseline_inf_sel_rank.py 16 clst le
python pkis1loto_baseline_inf_sel_rank.py 16 clst we
```

run frequent-hitters baselines:
```
python pkis1loto_baseline_inf_sel_rank.py 16 prom se
python pkis1loto_baseline_inf_sel_rank.py 16 prom le
python pkis1loto_baseline_inf_sel_rank.py 16 prom we
```

# re-run PKIS1 LOTO using Morgan fingerprints with radius=3, 1024 bits

re-run chemometric baselines:
```
python pkis1loto_baseline_inf_sel_rank_ecfp6.py 16 clst se
python pkis1loto_baseline_inf_sel_rank_ecfp6.py 16 clst le
python pkis1loto_baseline_inf_sel_rank_ecfp6.py 16 clst we
```

re-run frequent-hitters baselines:
```
python pkis1loto_baseline_inf_sel_rank_ecfp6.py 16 prom se
python pkis1loto_baseline_inf_sel_rank_ecfp6.py 16 prom le
python pkis1loto_baseline_inf_sel_rank_ecfp6.py 16 prom we
```

# Evaluate VS metrics using rankings for baselines and non-baselines.

Run metrics calculations over all targets for all IBRs--Huikun's (CS, AS), Ching-Pei's (RS), and my 6 baselines (Bx) using ECFP4.
```
python pkis1loto_eval_rocaucs.py
python pkis1loto_eval_nef10.py
python pkis1loto_eval_fasr10.py
python pkis1loto_eval_F1_MCC_pkis1median_thresh.py
```

> Note: for metrics evaluations, if no informers are returned as "active" (based on inferred active threshold), baseline ranking methods 's' and 'l' will use single most active compound for nearest-neighbors rankings (expansion).

Re-evaluate models on PKIS1 LOTO for baseline models using ECFP6 (not significantly different).
```
python pkis1loto_eval_rocaucs_ecfp6.py
python pkis1loto_eval_nef10_ecfp6.py
python pkis1loto_eval_fasr10_ecfp6.py
python pkis1loto_eval_F1_MCC_pkis1median_thresh_ecfp6.py
```

Run non-parametric paired t-tests (Wilcoxon signed-rank test)
```
python pkis1loto_nonparam_ttest_pvals_methods.py ../output_pkis1loto/metrics/pkis1loto_eval_ROCAUC.csv
python pkis1loto_nonparam_ttest_pvals_methods.py ../output_pkis1loto/metrics/pkis1loto_eval_NEF10.csv
python pkis1loto_nonparam_ttest_pvals_methods.py ../output_pkis1loto/metrics/pkis1loto_eval_FASR10.csv
python pkis1loto_nonparam_ttest_pvals_methods.py ../output_pkis1loto/metrics/pkis1loto_eval_F1_pkis1median_thresh.csv
python pkis1loto_nonparam_ttest_pvals_methods.py ../output_pkis1loto/metrics/pkis1loto_eval_MCC_pkis1median_thresh.csv
```

