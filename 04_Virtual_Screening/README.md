# Summary

This directory contains:
- Python [modules](modules) with methods used in this directory
- The [notebook](modelTrain-adora1_nr3c2.ipynb) for training the $\mathrm{A}_{1}\mathrm{AR}$ and $\mathrm{MR}$ models
- The hypermarameter optimization studies, stored in the [optuna_studies_dbs](optuna_studies_dbs) directory
- The optimized $\mathrm{MR}$ and $\mathrm{A}_{1}\mathrm{AR}$ serialized models, found in the [qsprModels](qsprModels) directory
- The [notebook](virtualScreening.ipynb) for performing the virtual screening
- The [notebook](ADPKD-ExplorationScreening-analysis.ipynb) for analyzing the spheroid swelling results on the ADPKD phenotypic assay


# Virtual Screening Results

These are the compounds selected for phenotypic screening based on the virtual screening results. For $\mathrm{A}_{1}\mathrm{AR}$, task variable consisted of solely pKi activity values, while for $\mathrm{MR}$, the training data consisted of both pIC50 and pKi values. These predicted values are represented as pChEMBL, as a placeholder for the -log10 transformed predicted activity values.

| Target   | Compound ID   | Structure                                               | Predicted pChEMBL value |
|:---------|:--------------|:--------------------------------------------------------|------------------------:|
| ADORA1   | 824745        | ![824745](../figures/mol_structures/824745.svg)         |                    8.45 |
| ADORA1   | 1237561       | ![1237561](../figures/mol_structures/1237561.svg)       |                    7.1  |
| ADORA1   | 1249141       | ![1249141](../figures/mol_structures/1249141.svg)       |                    7.27 |
| ADORA1   | 1823372       | ![1823372](../figures/mol_structures/1823372.svg)       |                    7.53 |
| ADORA1   | 22755240      | ![22755240](../figures/mol_structures/22755240.svg)     |                    7.13 |
| ADORA1   | 27070328      | ![27070328](../figures/mol_structures/27070328.svg)     |                    7.75 |
| NR3C2    | Z95680027     | ![Z95680027](../figures/mol_structures/Z95680027.svg)   |                    7.32 |
| NR3C2    | Z318400112    | ![Z318400112](../figures/mol_structures/Z318400112.svg) |                    7.33 |
| NR3C2    | Z90308949     | ![Z90308949](../figures/mol_structures/Z90308949.svg)   |                    7.09 |

# Statistical Test Results

- Results for all the statistical comparisons represented in the plots available in the notebook [ADPKD-ExplorationScreening-analysis](ADPKD-ExplorationScreening-analysis.ipynb).

<div align="center">
  <img src="../figures/adpkd-exploration-screen/boxplot-MannWhitneyU-plate2-A1R_DPCPX_FSK2-5.svg" alt="" width=480>
</div>

## DPCPX - Simulant (FSK) dose: 2.5µM
| group1                     | group2                   | pvalue                | symbol   | test_description                     | target                      |
|:---------------------------|:-------------------------|:----------------------|:---------|:-------------------------------------|:----------------------------|
| FSK 2.5µM, N=12            | DPCPX 0.001µM, N=4       | $2.12 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| DPCPX 0.001µM CPA 1µM, N=4 | DPCPX 0.1µM CPA 1µM, N=4 | $1.14 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| DPCPX 0.001µM, N=4         | DPCPX 0.1µM, N=4         | $6.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| DPCPX 0.1µM, N=4           | DPCPX 1µM, N=4           | $4.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, N=12            | DPCPX 0.1µM, N=4         | $2.97 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| DPCPX 0.001µM CPA 1µM, N=4 | DPCPX 1µM CPA 1µM, N=4   | $2.86 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| DPCPX 0.001µM, N=4         | DPCPX 1µM, N=4           | $1.14 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, N=12            | DPCPX 1µM, N=4           | $1.32 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |

<div align="center">
  <img src="../figures/adpkd-exploration-screen/boxplot-MannWhitneyU-plate2-A1R_Capadenoson_FSK2-5.svg" alt="" width=480>
</div>

## Capadenoson - Simulant (FSK) dose: 2.5µM
| group1                           | group2                         | pvalue                | symbol   | test_description                     | target                      |
|:---------------------------------|:-------------------------------|:----------------------|:---------|:-------------------------------------|:----------------------------|
| FSK 2.5µM, N=12                  | Capadenoson 0.001µM, N=4       | $1.10 \times 10^{-3}$ | **       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| Capadenoson 0.001µM CPA 1µM, N=4 | Capadenoson 0.1µM CPA 1µM, N=4 | $2.00 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| Capadenoson 0.001µM, N=4         | Capadenoson 0.1µM, N=4         | $2.86 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| Capadenoson 0.1µM, N=4           | Capadenoson 1µM, N=4           | $4.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, N=12                  | Capadenoson 0.1µM, N=4         | $1.10 \times 10^{-3}$ | **       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| Capadenoson 0.001µM CPA 1µM, N=4 | Capadenoson 1µM CPA 1µM, N=4   | $2.00 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| Capadenoson 0.001µM, N=4         | Capadenoson 1µM, N=4           | $2.86 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, N=12                  | Capadenoson 1µM, N=4           | $1.10 \times 10^{-3}$ | **       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |

<div align="center">
  <img src="../figures/adpkd-exploration-screen/boxplot-MannWhitneyU-plate2-A1R_MIPS521_FSK2-5.svg" alt="" width=480>
</div>

## MIPS521 - Simulant (FSK) dose: 2.5µM
| group1                       | group2                     | pvalue                | symbol   | test_description                     | target                      |
|:-----------------------------|:---------------------------|:----------------------|:---------|:-------------------------------------|:----------------------------|
| FSK 2.5µM, N=12              | MIPS521 0.001µM, N=4       | $4.46 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| MIPS521 0.001µM CPA 1µM, N=4 | MIPS521 0.1µM CPA 1µM, N=4 | $2.86 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| MIPS521 0.001µM, N=4         | MIPS521 0.1µM, N=4         | $2.86 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| MIPS521 0.1µM, N=4           | MIPS521 1µM, N=4           | $2.86 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, N=12              | MIPS521 0.1µM, N=4         | $1.10 \times 10^{-3}$ | **       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| MIPS521 0.001µM CPA 1µM, N=4 | MIPS521 1µM CPA 1µM, N=4   | $2.86 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| MIPS521 0.001µM, N=4         | MIPS521 1µM, N=4           | $2.86 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, N=12              | MIPS521 1µM, N=4           | $1.10 \times 10^{-3}$ | **       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |

<div align="center">
  <img src="../figures/adpkd-exploration-screen/boxplot-MannWhitneyU-plate3-A1R_824745_FSK2-5.svg" alt="" width=480>
</div>

## 824745 - Simulant (FSK) dose: 2.5µM
| group1                      | group2                    | pvalue                | symbol   | test_description                     | target                      |
|:----------------------------|:--------------------------|:----------------------|:---------|:-------------------------------------|:----------------------------|
| FSK 2.5µM, N=16             | 824745 0.001µM, N=4       | $6.17 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 824745 0.001µM CPA 1µM, N=4 | 824745 0.1µM CPA 1µM, N=4 | $2.00 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 824745 0.001µM, N=4         | 824745 0.1µM, N=4         | $8.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 824745 0.1µM, N=4           | 824745 1µM, N=4           | $1.14 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, N=16             | 824745 0.1µM, N=4         | $4.94 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 824745 0.001µM CPA 1µM, N=4 | 824745 1µM CPA 1µM, N=4   | $4.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 824745 0.001µM, N=4         | 824745 1µM, N=4           | $1.14 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, N=16             | 824745 1µM, N=4           | $8.01 \times 10^{-2}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |

<div align="center">
  <img src="../figures/adpkd-exploration-screen/boxplot-MannWhitneyU-plate3-A1R_1237561_FSK2-5.svg" alt="" width=480>
</div>

## 1237561 - Simulant (FSK) dose: 2.5µM
| group1                       | group2                     | pvalue                | symbol   | test_description                     | target                      |
|:-----------------------------|:---------------------------|:----------------------|:---------|:-------------------------------------|:----------------------------|
| FSK 2.5µM, N=16              | 1237561 0.001µM, N=4       | $2.49 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 1237561 0.001µM CPA 1µM, N=4 | 1237561 0.1µM CPA 1µM, N=4 | $3.43 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 1237561 0.001µM, N=4         | 1237561 0.1µM, N=4         | $4.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 1237561 0.1µM, N=4           | 1237561 1µM, N=4           | $3.43 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, N=16              | 1237561 0.1µM, N=4         | $9.63 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 1237561 0.001µM CPA 1µM, N=4 | 1237561 1µM CPA 1µM, N=4   | $2.86 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 1237561 0.001µM, N=4         | 1237561 1µM, N=4           | $8.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, N=16              | 1237561 1µM, N=4           | $2.90 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |

<div align="center">
  <img src="../figures/adpkd-exploration-screen/boxplot-MannWhitneyU-plate3-A1R_1249141_FSK2-5.svg" alt="" width=480>
</div>

## 1249141 - Simulant (FSK) dose: 2.5µM
| group1                       | group2                     | pvalue                | symbol   | test_description                     | target                      |
|:-----------------------------|:---------------------------|:----------------------|:---------|:-------------------------------------|:----------------------------|
| FSK 2.5µM, N=16              | 1249141 0.001µM, N=4       | $4.37 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 1249141 0.001µM CPA 1µM, N=4 | 1249141 0.1µM CPA 1µM, N=4 | $8.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 1249141 0.001µM, N=4         | 1249141 0.1µM, N=4         | $3.43 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 1249141 0.1µM, N=4           | 1249141 1µM, N=4           | $2.86 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, N=16              | 1249141 0.1µM, N=4         | $9.63 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 1249141 0.001µM CPA 1µM, N=4 | 1249141 1µM CPA 1µM, N=4   | $2.00 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 1249141 0.001µM, N=4         | 1249141 1µM, N=4           | $1.14 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, N=16              | 1249141 1µM, N=4           | $2.89 \times 10^{-3}$ | **       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |

<div align="center">
  <img src="../figures/adpkd-exploration-screen/boxplot-MannWhitneyU-plate3-A1R_1823372_FSK2-5.svg" alt="" width=480>
</div>

## 1823372 - Simulant (FSK) dose: 2.5µM
| group1                       | group2                     | pvalue                | symbol   | test_description                     | target                      |
|:-----------------------------|:---------------------------|:----------------------|:---------|:-------------------------------------|:----------------------------|
| FSK 2.5µM, N=16              | 1823372 0.001µM, N=4       | $2.19 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 1823372 0.001µM CPA 1µM, N=4 | 1823372 0.1µM CPA 1µM, N=4 | $6.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 1823372 0.001µM, N=4         | 1823372 0.1µM, N=4         | $1.14 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 1823372 0.1µM, N=4           | 1823372 1µM, N=4           | $6.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, N=16              | 1823372 0.1µM, N=4         | $6.82 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 1823372 0.001µM CPA 1µM, N=4 | 1823372 1µM CPA 1µM, N=4   | $1.14 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 1823372 0.001µM, N=4         | 1823372 1µM, N=4           | $1.14 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, N=16              | 1823372 1µM, N=4           | $2.90 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |

<div align="center">
  <img src="../figures/adpkd-exploration-screen/boxplot-MannWhitneyU-plate3-A1R_22755240_FSK2-5.svg" alt="" width=480>
</div>

## 22755240 - Simulant (FSK) dose: 2.5µM
| group1                        | group2                      | pvalue                | symbol   | test_description                     | target                      |
|:------------------------------|:----------------------------|:----------------------|:---------|:-------------------------------------|:----------------------------|
| FSK 2.5µM, N=16               | 22755240 0.001µM, N=4       | $1.78 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 22755240 0.001µM CPA 1µM, N=4 | 22755240 0.1µM CPA 1µM, N=4 | $2.00 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 22755240 0.001µM, N=4         | 22755240 0.1µM, N=4         | $8.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 22755240 0.1µM, N=4           | 22755240 1µM, N=4           | $3.43 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, N=16               | 22755240 0.1µM, N=4         | $4.94 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 22755240 0.001µM CPA 1µM, N=4 | 22755240 1µM CPA 1µM, N=4   | $8.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 22755240 0.001µM, N=4         | 22755240 1µM, N=4           | $6.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, N=16               | 22755240 1µM, N=4           | $9.95 \times 10^{-2}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |

<div align="center">
  <img src="../figures/adpkd-exploration-screen/boxplot-MannWhitneyU-plate3-A1R_27070328_FSK2-5.svg" alt="" width=480>
</div>

## 27070328 - Simulant (FSK) dose: 2.5µM
| group1                        | group2                      | pvalue                | symbol   | test_description                     | target                      |
|:------------------------------|:----------------------------|:----------------------|:---------|:-------------------------------------|:----------------------------|
| FSK 2.5µM, N=16               | 27070328 0.001µM, N=4       | $1.78 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 27070328 0.001µM CPA 1µM, N=4 | 27070328 0.1µM CPA 1µM, N=4 | $5.71 \times 10^{-2}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 27070328 0.001µM, N=4         | 27070328 0.1µM, N=4         | $2.86 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 27070328 0.1µM, N=4           | 27070328 1µM, N=4           | $2.86 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, N=16               | 27070328 0.1µM, N=4         | $1.22 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 27070328 0.001µM CPA 1µM, N=4 | 27070328 1µM CPA 1µM, N=4   | $2.86 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| 27070328 0.001µM, N=4         | 27070328 1µM, N=4           | $2.86 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, N=16               | 27070328 1µM, N=4           | $1.11 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | $\mathrm{A}_{1}\mathrm{AR}$ |

<div align="center">
  <img src="../figures/adpkd-exploration-screen/boxplot-MannWhitneyU-plate1-MR_Finerenone_FSK2-5.svg" alt="" width=480>
</div>

## Finerenone - Simulant (FSK) dose: 2.5µM
| group1                                  | group2                                | pvalue                | symbol   | test_description                     | target   |
|:----------------------------------------|:--------------------------------------|:----------------------|:---------|:-------------------------------------|:---------|
| FSK 2.5µM, N=16                         | Finerenone 0.001µM, N=4               | $1.22 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Finerenone 0.001µM Aldosterone 1µM, N=4 | Finerenone 0.1µM Aldosterone 1µM, N=4 | $6.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Finerenone 0.001µM, N=4                 | Finerenone 0.1µM, N=4                 | $8.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Finerenone 0.1µM, N=4                   | Finerenone 1µM, N=4                   | $4.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| FSK 2.5µM, N=16                         | Finerenone 0.1µM, N=4                 | $6.40 \times 10^{-2}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Finerenone 0.001µM Aldosterone 1µM, N=4 | Finerenone 1µM Aldosterone 1µM, N=4   | $2.86 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Finerenone 0.001µM, N=4                 | Finerenone 1µM, N=4                   | $3.43 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| FSK 2.5µM, N=16                         | Finerenone 1µM, N=4                   | $1.11 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | MR       |

<div align="center">
  <img src="../figures/adpkd-exploration-screen/boxplot-MannWhitneyU-plate1-MR_Esaxerenone_FSK2-5.svg" alt="" width=480>
</div>

## Esaxerenone - Simulant (FSK) dose: 2.5µM
| group1                                   | group2                                 | pvalue                | symbol   | test_description                     | target   |
|:-----------------------------------------|:---------------------------------------|:----------------------|:---------|:-------------------------------------|:---------|
| FSK 2.5µM, N=16                          | Esaxerenone 0.001µM, N=4               | $2.90 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Esaxerenone 0.001µM Aldosterone 1µM, N=4 | Esaxerenone 0.1µM Aldosterone 1µM, N=4 | $1.14 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Esaxerenone 0.001µM, N=4                 | Esaxerenone 0.1µM, N=4                 | $6.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Esaxerenone 0.1µM, N=4                   | Esaxerenone 1µM, N=4                   | $2.00 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| FSK 2.5µM, N=16                          | Esaxerenone 0.1µM, N=4                 | $8.92 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Esaxerenone 0.001µM Aldosterone 1µM, N=4 | Esaxerenone 1µM Aldosterone 1µM, N=4   | $5.71 \times 10^{-2}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Esaxerenone 0.001µM, N=4                 | Esaxerenone 1µM, N=4                   | $2.86 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | MR       |
| FSK 2.5µM, N=16                          | Esaxerenone 1µM, N=4                   | $2.93 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | MR       |

<div align="center">
  <img src="../figures/adpkd-exploration-screen/boxplot-MannWhitneyU-plate1-MR_Apararenone_FSK2-5.svg" alt="" width=480>
</div>

## Apararenone - Simulant (FSK) dose: 2.5µM
| group1                                   | group2                                 | pvalue                | symbol   | test_description                     | target   |
|:-----------------------------------------|:---------------------------------------|:----------------------|:---------|:-------------------------------------|:---------|
| FSK 2.5µM, N=16                          | Apararenone 0.001µM, N=4               | $4.94 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Apararenone 0.001µM Aldosterone 1µM, N=4 | Apararenone 0.1µM Aldosterone 1µM, N=4 | $4.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Apararenone 0.001µM, N=4                 | Apararenone 0.1µM, N=4                 | $8.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Apararenone 0.1µM, N=4                   | Apararenone 1µM, N=4                   | $8.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| FSK 2.5µM, N=16                          | Apararenone 0.1µM, N=4                 | $9.63 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Apararenone 0.001µM Aldosterone 1µM, N=4 | Apararenone 1µM Aldosterone 1µM, N=4   | $2.00 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Apararenone 0.001µM, N=4                 | Apararenone 1µM, N=4                   | $1.00$                | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| FSK 2.5µM, N=16                          | Apararenone 1µM, N=4                   | $5.54 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |

<div align="center">
  <img src="../figures/adpkd-exploration-screen/boxplot-MannWhitneyU-plate1-MR_Benidipine_FSK2-5.svg" alt="" width=480>
</div>

## Benidipine - Simulant (FSK) dose: 2.5µM
| group1                                  | group2                                | pvalue                | symbol   | test_description                     | target   |
|:----------------------------------------|:--------------------------------------|:----------------------|:---------|:-------------------------------------|:---------|
| FSK 2.5µM, N=16                         | Benidipine 0.001µM, N=4               | $3.85 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Benidipine 0.001µM Aldosterone 1µM, N=4 | Benidipine 0.1µM Aldosterone 1µM, N=4 | $8.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Benidipine 0.001µM, N=4                 | Benidipine 0.1µM, N=4                 | $8.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Benidipine 0.1µM, N=4                   | Benidipine 1µM, N=4                   | $1.14 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| FSK 2.5µM, N=16                         | Benidipine 0.1µM, N=4                 | $4.37 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Benidipine 0.001µM Aldosterone 1µM, N=4 | Benidipine 1µM Aldosterone 1µM, N=4   | $3.43 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Benidipine 0.001µM, N=4                 | Benidipine 1µM, N=4                   | $2.00 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| FSK 2.5µM, N=16                         | Benidipine 1µM, N=4                   | $2.90 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |

<div align="center">
  <img src="../figures/adpkd-exploration-screen/boxplot-MannWhitneyU-plate1-MR_Z318400112_FSK2-5.svg" alt="" width=480>
</div>

## Z318400112 - Simulant (FSK) dose: 2.5µM
| group1                                  | group2                                | pvalue                | symbol   | test_description                     | target   |
|:----------------------------------------|:--------------------------------------|:----------------------|:---------|:-------------------------------------|:---------|
| FSK 2.5µM, N=16                         | Z318400112 0.001µM, N=4               | $4.94 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Z318400112 0.001µM Aldosterone 1µM, N=4 | Z318400112 0.1µM Aldosterone 1µM, N=4 | $2.86 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Z318400112 0.001µM, N=4                 | Z318400112 0.1µM, N=4                 | $6.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Z318400112 0.1µM, N=4                   | Z318400112 1µM, N=4                   | $1.00$                | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| FSK 2.5µM, N=16                         | Z318400112 0.1µM, N=4                 | $1.48 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Z318400112 0.001µM Aldosterone 1µM, N=4 | Z318400112 1µM Aldosterone 1µM, N=4   | $4.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Z318400112 0.001µM, N=4                 | Z318400112 1µM, N=4                   | $4.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| FSK 2.5µM, N=16                         | Z318400112 1µM, N=4                   | $1.48 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |

<div align="center">
  <img src="../figures/adpkd-exploration-screen/boxplot-MannWhitneyU-plate1-MR_Z90308949_FSK2-5.svg" alt="" width=480>
</div>

## Z90308949 - Simulant (FSK) dose: 2.5µM
| group1                                 | group2                               | pvalue                | symbol   | test_description                     | target   |
|:---------------------------------------|:-------------------------------------|:----------------------|:---------|:-------------------------------------|:---------|
| FSK 2.5µM, N=16                        | Z90308949 0.001µM, N=4               | $8.01 \times 10^{-2}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Z90308949 0.001µM Aldosterone 1µM, N=4 | Z90308949 0.1µM Aldosterone 1µM, N=4 | $8.86 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Z90308949 0.001µM, N=4                 | Z90308949 0.1µM, N=4                 | $2.00 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Z90308949 0.1µM, N=4                   | Z90308949 1µM, N=4                   | $2.00 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| FSK 2.5µM, N=16                        | Z90308949 0.1µM, N=4                 | $9.63 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Z90308949 0.001µM Aldosterone 1µM, N=4 | Z90308949 1µM Aldosterone 1µM, N=4   | $2.00 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Z90308949 0.001µM, N=4                 | Z90308949 1µM, N=4                   | $5.71 \times 10^{-2}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| FSK 2.5µM, N=16                        | Z90308949 1µM, N=4                   | $1.48 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |

<div align="center">
  <img src="../figures/adpkd-exploration-screen/boxplot-MannWhitneyU-plate1-MR_Z95680027_FSK2-5.svg" alt="" width=480>
</div>

## Z95680027 - Simulant (FSK) dose: 2.5µM
| group1                                 | group2                               | pvalue                | symbol   | test_description                     | target   |
|:---------------------------------------|:-------------------------------------|:----------------------|:---------|:-------------------------------------|:---------|
| FSK 2.5µM, N=16                        | Z95680027 0.001µM, N=4               | $5.54 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Z95680027 0.001µM Aldosterone 1µM, N=4 | Z95680027 0.1µM Aldosterone 1µM, N=4 | $2.86 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Z95680027 0.001µM, N=4                 | Z95680027 0.1µM, N=4                 | $1.00$                | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Z95680027 0.1µM, N=4                   | Z95680027 1µM, N=4                   | $2.86 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | MR       |
| FSK 2.5µM, N=16                        | Z95680027 0.1µM, N=4                 | $6.17 \times 10^{-1}$ | ns       | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Z95680027 0.001µM Aldosterone 1µM, N=4 | Z95680027 1µM Aldosterone 1µM, N=4   | $2.86 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | MR       |
| Z95680027 0.001µM, N=4                 | Z95680027 1µM, N=4                   | $2.86 \times 10^{-2}$ | *        | Mann-Whitney-Wilcoxon test two-sided | MR       |
| FSK 2.5µM, N=16                        | Z95680027 1µM, N=4                   | $7.43 \times 10^{-3}$ | **       | Mann-Whitney-Wilcoxon test two-sided | MR       |
