Zenodo DOI Badge: https://doi.org/10.5281/zenodo.3354432


# Predicting kinase inhibitors using bioactivity matrix derived informer sets
This repository contains baseline Informer-Based Ranking (IBR) methods and procedures for evaluating metrics for IBR performance (including non-baseline IBRs).
It provides the kinase screening data used to evaluate the IBR methods.
See also the repositories for:
- [Regression Selection (RS)](https://github.com/leepei/informer)
- [Coding Selection (CS) and Adaptive Selection (AS)](https://github.com/wiscstatman/esdd/tree/master/informRset)

## Citation

If you use this software or the new chemical screening data, please cite:

Huikun Zhang<sup>+</sup>, Spencer S Ericksen<sup>+</sup>, Ching-pei Lee<sup>+</sup>, Gene E Ananiev, Nathan Wlodarchak, Julie C Mitchell, Anthony Gitter, Stephen J Wright, F Michael Hoffmann, Scott A Wildman, Michael A Newton.
["Predicting kinase inhibitors using bioactivity matrix derived informer sets."](https://doi.org/10.1101/532762)
*bioRxiv* 2019. doi:10.1101/532762.

<sup>+</sup> equal contributions

## File Structure

The `informers` repository comprises 5 main folders:

- `source` - code for running baseline methods and running evaluation metrics for all IBRs (both validation and new targets)

- `data` - new screening data for PknB and BGLF4, pre-processed PKIS1 and PKIS2 screening data
  - `data/compounds` - compound SMILES, Morgan fingerprints, and Morgan Jaccard distance matrices
  - `data/thresholds_2sigma` - inferred target activity thresholds for assigning compound binary activity labels
  - `data/original_data` - original PKIS1 and PKIS2 data sets with descriptions of pre-processing
  - `data/rop18` - PKIS1 activity data from assays on Toxoplasma gondii Rhoptry Kinase ROP18, [Simpson et al. 2016](https://doi.org/10.1021/acsinfecdis.5b00102)
  
- `output_newtargs` - output from all IBR methods on prospective microbial kinase targets (PknB, BGLF, ROP18) and metrics evaluations

- `output_pkis1loto` - output from all IBR methods for 224 PKIS1 targets and metrics evaluations

- `figures` - codes for plotting figures



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


