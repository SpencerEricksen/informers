# Predicting kinase inhibitors using bioactivity matrix derived informer sets
This repository contains baseline Informer-Based Ranking (IBR) methods and procedures for evaluating metrics for IBR performance (including non-baseline IBRs).
It provides the kinase screening data used to evaluate the IBR methods.
See also the repositories for:
- [Regression Selection (RS)](https://github.com/leepei/informer)
- [Coding Selection (CS) and Adaptive Selection (AS)](https://github.com/wiscstatman/esdd/tree/master/informRset)

## Kinase screening data
There are two main data directories
- `inf_newtargs/data/` - new screening data for PknB and BGLF4
- `inf_pkis1loto/data/` - pre-processed PKIS1 and PKIS2 screening data

## Python environment
This code was run in the following conda environment:

Python version: 2.7.15

Packages in environment:

| Name         | Version     |
| ------------ | ----------- |
| matplotlib   | 2.2.2       |
| numpy        | 1.14.5      |
| pandas       | 0.23.3      |
| rdkit        | 2018.03.2   |
| scikit-learn | 0.19.2      |
| scipy        | 1.1.0       |
| seaborn      | 0.9.0       |
