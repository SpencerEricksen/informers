### PKIS1 LOTO

run chemometric baselines:\
`./pkis1loto_baseline_inf_sel_rank.py 16 clst se`\
`./pkis1loto_baseline_inf_sel_rank.py 16 clst le`\
`./pkis1loto_baseline_inf_sel_rank.py 16 clst we`

run frequent-hitters baselines\
`./pkis1loto_baseline_inf_sel_rank.py 16 prom se`\
`./pkis1loto_baseline_inf_sel_rank.py 16 prom le`\
`./pkis1loto_baseline_inf_sel_rank.py 16 prom we`


### NEW TARGETS: PknB and BGLF4
run all baselines on new targets using either PKIS1 or PKIS2 matrix:\
`./newtargs_baseline_inf_sel_rank.py 16 1`\
`./newtargs_baseline_inf_sel_rank.py 16 2`

