# Predicting kinase inhibitors using bioactivity matrix derived informer sets
This repository contains baseline Informer-Based Ranking (IBR) methods and evaluation metrics for all IBRs, including non-baselines.
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
```
python 2.7.15

# packages in environment at ~/anaconda2:
# Name                    Version                   Build    Channel

matplotlib                2.2.2                    py27_1    conda-forge
numpy                     1.14.5          py27_blas_openblashd3ea46f_201  [blas_openblas]  conda-forge
pandas                    0.23.3                   py27_0    conda-forge
pubchempy                 1.0.4                      py_2    mcs07
rdkit                     2018.03.2        py27h4ce040e_0    conda-forge
scikit-learn              0.19.2          py27_blas_openblas_200  [blas_openblas]  conda-forge
scipy                     1.1.0           py27_blas_openblashd3ea46f_201  [blas_openblas]  conda-forge
seaborn                   0.9.0                      py_0    conda-forge
```
