### run metrics using rankings for baselines and non-baselines \
`./pkis1loto_eval_rocaucs.py` \
`./pkis1loto_eval_nef10.py` \
`./pkis1loto_eval_fasr10.py`

> Note: for metrics evaluations, if no informers are returned as "active" (based on inferred active threshold), baseline ranking methods 's' and 'l' will fail to rank compounds. In such cases, the target is given a minimal metric score (ROCAUC=0.5, NEF10=0.5, or FASR10=0.0)

### run t-tests \
`./pkis1loto_ttest_pvals_methods.py pkis1loto_eval_ROCAUC.csv` \
`./pkis1loto_ttest_pvals_methods.py pkis1loto_eval_NEF10.csv` \
`./pkis1loto_ttest_pvals_methods.py pkis1loto_eval_FASR10.csv`




