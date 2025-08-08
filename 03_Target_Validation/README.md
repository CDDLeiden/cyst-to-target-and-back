# Summary

This directory contains the results of the target validation screening experiement, where we screened compounds with high affinities to our targets of interest in the ADPKD phenotypic assay. The analysis of the obtained data is performed in the notebook [ADPKD-TargetValidationScreening-analysis](ADPKD-TargetValidationScreening-analysis.ipynb).

Compounds listed herein were prioritized as described in the manuscript and on the notebooks within [02_TargetID_and_Prioritization](../02_TargetID_and_Prioritization/). For further transparency, we also report the affinity values of these screened compounds on our target of interest. The data used here is from ChEMBL.

For the full data, refer to [chembl_data_validation_compounds.csv](../figures/mol_structures/chembl_data_validation_compounds.csv).

# Known affinities of the screened compounds to targets of interest

| Target   | Compound Name   | img                                                         | assay_type   | standard_type   |   pchembl_value_mean |   pchembl_value_std |
|:---------|:----------------|:------------------------------------------------------------|:-------------|:----------------|---------------------:|--------------------:|
| ADORA1   | CPA             | ![CPA](../figures/mol_structures/CPA.svg)                   | B            | Ki              |                 8.41 |                0.46 |
| ADORA1   | CPA             | ![CPA](../figures/mol_structures/CPA.svg)                   | F            | IC50            |                 8.57 |                0    |
| ADORA1   | Capadenoson     | ![Capadenoson](../figures/mol_structures/Capadenoson.svg)   | B            | Ki              |                 8.85 |                0    |
| ADORA1   | DPCPX           | ![DPCPX](../figures/mol_structures/DPCPX.svg)               | B            | Ki              |                 8.55 |                0.5  |
| ADORA1   | DPCPX           | ![DPCPX](../figures/mol_structures/DPCPX.svg)               | B            | IC50            |                 8.46 |                0.66 |
| ADORA1   | DPCPX           | ![DPCPX](../figures/mol_structures/DPCPX.svg)               | F            | Ki              |                 9.19 |                0.14 |
| NR3C2    | Aldosterone     | ![Aldosterone](../figures/mol_structures/Aldosterone.svg)   | B            | IC50            |                 9.52 |                0    |
| NR3C2    | Aldosterone     | ![Aldosterone](../figures/mol_structures/Aldosterone.svg)   | F            | IC50            |                 8.03 |                0    |
| NR3C2    | Esaxerenone     | ![Esaxerenone](../figures/mol_structures/Esaxerenone.svg)   | B            | IC50            |                 8.03 |                0    |
| NR3C2    | Esaxerenone     | ![Esaxerenone](../figures/mol_structures/Esaxerenone.svg)   | F            | IC50            |                 8.62 |                0    |
| NR3C2    | Finerenone      | ![Finerenone](../figures/mol_structures/Finerenone.svg)     | B            | IC50            |                 7.58 |                0.24 |
| NR3C2    | Finerenone      | ![Finerenone](../figures/mol_structures/Finerenone.svg)     | F            | IC50            |                 7.8  |                0    |
| P2X7     | A-804598        | ![A-804598](../figures/mol_structures/A-804598.svg)         | B            | IC50            |                 7.66 |                0.43 |
| P2X7     | A-804598        | ![A-804598](../figures/mol_structures/A-804598.svg)         | B            | Ki              |                 8.05 |                0    |
| P2X7     | JNJ-47965567    | ![JNJ-47965567](../figures/mol_structures/JNJ-47965567.svg) | B            | Ki              |                 7.9  |                0    |
| P2X7     | JNJ-47965567    | ![JNJ-47965567](../figures/mol_structures/JNJ-47965567.svg) | B            | IC50            |                 7.57 |                0.52 |
| SLC2A1   | Z211311146      | ![Z211311146](../figures/mol_structures/Z211311146.svg)     | B            | IC50            |                 7.57 |              nan    |
| SLC2A1   | Z4509024390     | ![Z4509024390](../figures/mol_structures/Z4509024390.svg)   | B            | IC50            |                 8.4  |              nan    |

# Statistical Test Results

- Results for all the statistical comparisons represented in the plots available in the notebook [ADPKD-TargetValidationScreening-analysis](ADPKD-TargetValidationScreening-analysis.ipynb).

<div align="center">
  <img src="../figures/adpkd-validation-screen/boxplot-MannWhitneyU-plate1-GLU_Z4509024390_FSK2-5.svg" alt="" width=480>
</div>

## Z4509024390 - Simulant (FSK) dose: 2.5µM
| group1              | group2                       |   pvalue | symbol   | target   |
|:--------------------|:-----------------------------|---------:|:---------|:---------|
| FSK 2.5µM, ($N=12$) | Z4509024390 0.001µM, ($N=4$) |   0.4462 | ns       | GLUT1    |
| FSK 2.5µM, ($N=12$) | Z4509024390 0.1µM, ($N=4$)   |   0.5209 | ns       | GLUT1    |
| FSK 2.5µM, ($N=12$) | Z4509024390 1.0µM, ($N=4$)   |   0.0132 | *        | GLUT1    |

<div align="center">
  <img src="../figures/adpkd-validation-screen/boxplot-MannWhitneyU-plate1-GLU_Z4509024390_FSK0-79.svg" alt="" width=480>
</div>

## Z4509024390 - Simulant (FSK) dose: 0.79µM
| group1              | group2                       |   pvalue | symbol   | target   |
|:--------------------|:-----------------------------|---------:|:---------|:---------|
| FSK 0.79µM, ($N=8$) | Z4509024390 0.001µM, ($N=4$) |   0.5697 | ns       | GLUT1    |
| FSK 0.79µM, ($N=8$) | Z4509024390 0.1µM, ($N=4$)   |   0.0727 | ns       | GLUT1    |
| FSK 0.79µM, ($N=8$) | Z4509024390 1.0µM, ($N=4$)   |   0.2141 | ns       | GLUT1    |

<div align="center">
  <img src="../figures/adpkd-validation-screen/boxplot-MannWhitneyU-plate1-GLU_Z211311146_FSK2-5.svg" alt="" width=480>
</div>

## Z211311146 - Simulant (FSK) dose: 2.5µM
| group1              | group2                      |   pvalue | symbol   | target   |
|:--------------------|:----------------------------|---------:|:---------|:---------|
| FSK 2.5µM, ($N=12$) | Z211311146 0.001µM, ($N=4$) |   0.1703 | ns       | GLUT1    |
| FSK 2.5µM, ($N=12$) | Z211311146 0.1µM, ($N=4$)   |   0.7703 | ns       | GLUT1    |
| FSK 2.5µM, ($N=12$) | Z211311146 1.0µM, ($N=4$)   |   0.0198 | *        | GLUT1    |

<div align="center">
  <img src="../figures/adpkd-validation-screen/boxplot-MannWhitneyU-plate1-GLU_Z211311146_FSK0-79.svg" alt="" width=480>
</div>

## Z211311146 - Simulant (FSK) dose: 0.79µM
| group1              | group2                      |   pvalue | symbol   | target   |
|:--------------------|:----------------------------|---------:|:---------|:---------|
| FSK 0.79µM, ($N=8$) | Z211311146 0.001µM, ($N=4$) |   0.6828 | ns       | GLUT1    |
| FSK 0.79µM, ($N=8$) | Z211311146 0.1µM, ($N=4$)   |   0.4606 | ns       | GLUT1    |
| FSK 0.79µM, ($N=8$) | Z211311146 1.0µM, ($N=4$)   |   0.6828 | ns       | GLUT1    |

<div align="center">
  <img src="../figures/adpkd-validation-screen/boxplot-MannWhitneyU-plate1-P2RX7_JNJ-47965567_FSK2-5.svg" alt="" width=480>
</div>

## JNJ-47965567 - Simulant (FSK) dose: 2.5µM
| group1              | group2                        |   pvalue | symbol   | target   |
|:--------------------|:------------------------------|---------:|:---------|:---------|
| FSK 2.5µM, ($N=12$) | JNJ-47965567 0.001µM, ($N=4$) |   0.8615 | ns       | P2X7     |
| FSK 2.5µM, ($N=12$) | JNJ-47965567 0.1µM, ($N=4$)   |   0.0132 | *        | P2X7     |
| FSK 2.5µM, ($N=12$) | JNJ-47965567 1.0µM, ($N=4$)   |   0.3165 | ns       | P2X7     |

<div align="center">
  <img src="../figures/adpkd-validation-screen/boxplot-MannWhitneyU-plate1-P2RX7_JNJ-47965567_FSK0-79.svg" alt="" width=480>
</div>

## JNJ-47965567 - Simulant (FSK) dose: 0.79µM
| group1              | group2                        |   pvalue | symbol   | target   |
|:--------------------|:------------------------------|---------:|:---------|:---------|
| FSK 0.79µM, ($N=8$) | JNJ-47965567 0.001µM, ($N=4$) |   0.5697 | ns       | P2X7     |
| FSK 0.79µM, ($N=8$) | JNJ-47965567 0.1µM, ($N=4$)   |   0.2828 | ns       | P2X7     |
| FSK 0.79µM, ($N=8$) | JNJ-47965567 1.0µM, ($N=4$)   |   0.3677 | ns       | P2X7     |

<div align="center">
  <img src="../figures/adpkd-validation-screen/boxplot-MannWhitneyU-plate1-P2RX7_A-804598_FSK2-5.svg" alt="" width=480>
</div>

## A-804598 - Simulant (FSK) dose: 2.5µM
| group1              | group2                    |   pvalue | symbol   | target   |
|:--------------------|:--------------------------|---------:|:---------|:---------|
| FSK 2.5µM, ($N=12$) | A-804598 0.001µM, ($N=4$) |   0.8615 | ns       | P2X7     |
| FSK 2.5µM, ($N=12$) | A-804598 0.1µM, ($N=4$)   |   0.2121 | ns       | P2X7     |
| FSK 2.5µM, ($N=12$) | A-804598 1.0µM, ($N=4$)   |   0.9527 | ns       | P2X7     |

<div align="center">
  <img src="../figures/adpkd-validation-screen/boxplot-MannWhitneyU-plate1-P2RX7_A-804598_FSK0-79.svg" alt="" width=480>
</div>

## A-804598 - Simulant (FSK) dose: 0.79µM
| group1              | group2                    |   pvalue | symbol   | target   |
|:--------------------|:--------------------------|---------:|:---------|:---------|
| FSK 0.79µM, ($N=8$) | A-804598 0.001µM, ($N=4$) |   0.1535 | ns       | P2X7     |
| FSK 0.79µM, ($N=8$) | A-804598 0.1µM, ($N=4$)   |   0.8081 | ns       | P2X7     |
| FSK 0.79µM, ($N=8$) | A-804598 1.0µM, ($N=4$)   |   0.5697 | ns       | P2X7     |

<div align="center">
  <img src="../figures/adpkd-validation-screen/boxplot-MannWhitneyU-plate2-A1R_DPCPX_FSK2-5_CPA.svg" alt="" width=720>
</div>

## DPCPX - Simulant (FSK) dose: 2.5µM
| group1                         | group2                           |   pvalue | symbol   | target                      |
|:-------------------------------|:---------------------------------|---------:|:---------|:----------------------------|
| FSK 2.5µM, ($N=11$)            | CPA 0.001µM, ($N=4$)             |   0.1773 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, ($N=11$)            | CPA 0.1µM, ($N=3$)               |   0.0879 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, ($N=11$)            | CPA 1µM, ($N=4$)                 |   0.0015 | **       | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, ($N=11$)            | DPCPX 0.001µM, ($N=4$)           |   0.2256 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, ($N=11$)            | DPCPX 0.1µM, ($N=4$)             |   0.104  | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, ($N=11$)            | DPCPX 1µM, ($N=4$)               |   0.0015 | **       | $\mathrm{A}_{1}\mathrm{AR}$ |
| DPCPX 0.1µM CPA 0.1µM, ($N=4$) | DPCPX 0.001µM CPA 1µM, ($N=4$)   |   0.6857 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| DPCPX 0.1µM CPA 0.1µM, ($N=4$) | DPCPX 0.1µM CPA 1µM, ($N=4$)     |   0.1143 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| DPCPX 1µM CPA 0.1µM, ($N=4$)   | DPCPX 1µM CPA 1µM, ($N=4$)       |   0.6857 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 0.1µM, ($N=3$)             | DPCPX 0.001µM CPA 0.1µM, ($N=4$) |   0.8571 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 0.1µM, ($N=3$)             | DPCPX 0.1µM CPA 0.1µM, ($N=4$)   |   0.1143 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 0.1µM, ($N=3$)             | DPCPX 1µM CPA 0.1µM, ($N=4$)     |   0.0571 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 1µM, ($N=4$)               | DPCPX 0.001µM CPA 1µM, ($N=4$)   |   0.1143 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 1µM, ($N=4$)               | DPCPX 0.1µM CPA 1µM, ($N=4$)     |   0.8857 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 1µM, ($N=4$)               | DPCPX 1µM CPA 1µM, ($N=4$)       |   0.0286 | *        | $\mathrm{A}_{1}\mathrm{AR}$ |
| DPCPX 0.001µM, ($N=4$)         | DPCPX 0.1µM CPA 0.1µM, ($N=4$)   |   1      | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| DPCPX 0.001µM, ($N=4$)         | DPCPX 0.1µM CPA 1µM, ($N=4$)     |   0.0286 | *        | $\mathrm{A}_{1}\mathrm{AR}$ |
| DPCPX 0.1µM, ($N=4$)           | DPCPX 0.1µM CPA 0.1µM, ($N=4$)   |   0.8857 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| DPCPX 0.1µM, ($N=4$)           | DPCPX 0.1µM CPA 1µM, ($N=4$)     |   0.0286 | *        | $\mathrm{A}_{1}\mathrm{AR}$ |
| DPCPX 1µM, ($N=4$)             | DPCPX 1µM CPA 0.1µM, ($N=4$)     |   0.4857 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| DPCPX 1µM, ($N=4$)             | DPCPX 1µM CPA 1µM, ($N=4$)       |   0.2    | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |

<div align="center">
  <img src="../figures/adpkd-validation-screen/boxplot-MannWhitneyU-plate2-A1R_DPCPX_FSK0-79_CPA.svg" alt="" width=720>
</div>

## DPCPX - Simulant (FSK) dose: 0.79µM
| group1                         | group2                           |   pvalue | symbol   | target                      |
|:-------------------------------|:---------------------------------|---------:|:---------|:----------------------------|
| FSK 0.79µM, ($N=8$)            | CPA 0.001µM, ($N=3$)             |   0.0485 | *        | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 0.79µM, ($N=8$)            | CPA 0.1µM, ($N=4$)               |   0.0283 | *        | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 0.79µM, ($N=8$)            | CPA 1µM, ($N=4$)                 |   0.0283 | *        | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 0.79µM, ($N=8$)            | DPCPX 0.001µM, ($N=4$)           |   0.4606 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 0.79µM, ($N=8$)            | DPCPX 0.1µM, ($N=4$)             |   0.1091 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 0.79µM, ($N=8$)            | DPCPX 1µM, ($N=4$)               |   0.0081 | **       | $\mathrm{A}_{1}\mathrm{AR}$ |
| DPCPX 0.1µM CPA 0.1µM, ($N=4$) | DPCPX 0.001µM CPA 1µM, ($N=4$)   |   0.0286 | *        | $\mathrm{A}_{1}\mathrm{AR}$ |
| DPCPX 0.1µM CPA 0.1µM, ($N=4$) | DPCPX 0.1µM CPA 1µM, ($N=4$)     |   0.1143 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| DPCPX 1µM CPA 0.1µM, ($N=4$)   | DPCPX 1µM CPA 1µM, ($N=4$)       |   0.3429 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 0.1µM, ($N=4$)             | DPCPX 0.001µM CPA 0.1µM, ($N=4$) |   0.0286 | *        | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 0.1µM, ($N=4$)             | DPCPX 0.1µM CPA 0.1µM, ($N=4$)   |   0.1143 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 0.1µM, ($N=4$)             | DPCPX 1µM CPA 0.1µM, ($N=4$)     |   0.0286 | *        | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 1µM, ($N=4$)               | DPCPX 0.001µM CPA 1µM, ($N=4$)   |   0.2    | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 1µM, ($N=4$)               | DPCPX 0.1µM CPA 1µM, ($N=4$)     |   0.6857 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 1µM, ($N=4$)               | DPCPX 1µM CPA 1µM, ($N=4$)       |   0.0286 | *        | $\mathrm{A}_{1}\mathrm{AR}$ |
| DPCPX 0.001µM, ($N=4$)         | DPCPX 0.1µM CPA 0.1µM, ($N=4$)   |   1      | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| DPCPX 0.001µM, ($N=4$)         | DPCPX 0.1µM CPA 1µM, ($N=4$)     |   0.1143 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| DPCPX 0.1µM, ($N=4$)           | DPCPX 0.1µM CPA 0.1µM, ($N=4$)   |   0.3429 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| DPCPX 0.1µM, ($N=4$)           | DPCPX 0.1µM CPA 1µM, ($N=4$)     |   0.0286 | *        | $\mathrm{A}_{1}\mathrm{AR}$ |
| DPCPX 1µM, ($N=4$)             | DPCPX 1µM CPA 0.1µM, ($N=4$)     |   1      | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| DPCPX 1µM, ($N=4$)             | DPCPX 1µM CPA 1µM, ($N=4$)       |   0.1143 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |

<div align="center">
  <img src="../figures/adpkd-validation-screen/boxplot-MannWhitneyU-plate2-A1R_Capadenoson_FSK2-5_CPA.svg" alt="" width=720>
</div>

## Capadenoson - Simulant (FSK) dose: 2.5µM
| group1                               | group2                                 |   pvalue | symbol   | target                      |
|:-------------------------------------|:---------------------------------------|---------:|:---------|:----------------------------|
| FSK 2.5µM, ($N=11$)                  | CPA 0.001µM, ($N=4$)                   |   0.1773 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, ($N=11$)                  | CPA 0.1µM, ($N=3$)                     |   0.0879 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, ($N=11$)                  | CPA 1µM, ($N=4$)                       |   0.0015 | **       | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, ($N=11$)                  | Capadenoson 0.001µM, ($N=4$)           |   0.4894 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, ($N=11$)                  | Capadenoson 0.1µM, ($N=4$)             |   0.0015 | **       | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 2.5µM, ($N=11$)                  | Capadenoson 1µM, ($N=3$)               |   0.0055 | **       | $\mathrm{A}_{1}\mathrm{AR}$ |
| Capadenoson 0.1µM CPA 0.1µM, ($N=4$) | Capadenoson 0.001µM CPA 1µM, ($N=4$)   |   0.0571 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| Capadenoson 0.1µM CPA 0.1µM, ($N=4$) | Capadenoson 0.1µM CPA 1µM, ($N=4$)     |   0.0571 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| Capadenoson 1µM CPA 0.1µM, ($N=4$)   | Capadenoson 1µM CPA 1µM, ($N=4$)       |   0.8857 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 0.1µM, ($N=3$)                   | Capadenoson 0.001µM CPA 0.1µM, ($N=4$) |   0.6286 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 0.1µM, ($N=3$)                   | Capadenoson 0.1µM CPA 0.1µM, ($N=4$)   |   0.6286 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 0.1µM, ($N=3$)                   | Capadenoson 1µM CPA 0.1µM, ($N=4$)     |   0.0571 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 1µM, ($N=4$)                     | Capadenoson 0.001µM CPA 1µM, ($N=4$)   |   0.3429 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 1µM, ($N=4$)                     | Capadenoson 0.1µM CPA 1µM, ($N=4$)     |   0.2    | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 1µM, ($N=4$)                     | Capadenoson 1µM CPA 1µM, ($N=4$)       |   0.0286 | *        | $\mathrm{A}_{1}\mathrm{AR}$ |
| Capadenoson 0.001µM, ($N=4$)         | Capadenoson 0.1µM CPA 0.1µM, ($N=4$)   |   0.0286 | *        | $\mathrm{A}_{1}\mathrm{AR}$ |
| Capadenoson 0.001µM, ($N=4$)         | Capadenoson 0.1µM CPA 1µM, ($N=4$)     |   0.0286 | *        | $\mathrm{A}_{1}\mathrm{AR}$ |
| Capadenoson 0.1µM, ($N=4$)           | Capadenoson 0.1µM CPA 0.1µM, ($N=4$)   |   0.2    | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| Capadenoson 0.1µM, ($N=4$)           | Capadenoson 0.1µM CPA 1µM, ($N=4$)     |   0.0286 | *        | $\mathrm{A}_{1}\mathrm{AR}$ |
| Capadenoson 1µM, ($N=3$)             | Capadenoson 1µM CPA 0.1µM, ($N=4$)     |   0.8571 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| Capadenoson 1µM, ($N=3$)             | Capadenoson 1µM CPA 1µM, ($N=4$)       |   0.8571 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |

<div align="center">
  <img src="../figures/adpkd-validation-screen/boxplot-MannWhitneyU-plate2-A1R_Capadenoson_FSK0-79_CPA.svg" alt="" width=720>
</div>

## Capadenoson - Simulant (FSK) dose: 0.79µM
| group1                               | group2                                 |   pvalue | symbol   | target                      |
|:-------------------------------------|:---------------------------------------|---------:|:---------|:----------------------------|
| FSK 0.79µM, ($N=8$)                  | CPA 0.001µM, ($N=3$)                   |   0.0485 | *        | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 0.79µM, ($N=8$)                  | CPA 0.1µM, ($N=4$)                     |   0.0283 | *        | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 0.79µM, ($N=8$)                  | CPA 1µM, ($N=4$)                       |   0.0283 | *        | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 0.79µM, ($N=8$)                  | Capadenoson 0.001µM, ($N=4$)           |   0.5697 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 0.79µM, ($N=8$)                  | Capadenoson 0.1µM, ($N=4$)             |   0.004  | **       | $\mathrm{A}_{1}\mathrm{AR}$ |
| FSK 0.79µM, ($N=8$)                  | Capadenoson 1µM, ($N=4$)               |   0.004  | **       | $\mathrm{A}_{1}\mathrm{AR}$ |
| Capadenoson 0.1µM CPA 0.1µM, ($N=4$) | Capadenoson 0.001µM CPA 1µM, ($N=4$)   |   0.0286 | *        | $\mathrm{A}_{1}\mathrm{AR}$ |
| Capadenoson 0.1µM CPA 0.1µM, ($N=4$) | Capadenoson 0.1µM CPA 1µM, ($N=4$)     |   0.1143 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| Capadenoson 1µM CPA 0.1µM, ($N=4$)   | Capadenoson 1µM CPA 1µM, ($N=4$)       |   0.1143 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 0.1µM, ($N=4$)                   | Capadenoson 0.001µM CPA 0.1µM, ($N=4$) |   0.2    | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 0.1µM, ($N=4$)                   | Capadenoson 0.1µM CPA 0.1µM, ($N=4$)   |   0.0286 | *        | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 0.1µM, ($N=4$)                   | Capadenoson 1µM CPA 0.1µM, ($N=4$)     |   0.0286 | *        | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 1µM, ($N=4$)                     | Capadenoson 0.001µM CPA 1µM, ($N=4$)   |   0.3429 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 1µM, ($N=4$)                     | Capadenoson 0.1µM CPA 1µM, ($N=4$)     |   0.2    | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| CPA 1µM, ($N=4$)                     | Capadenoson 1µM CPA 1µM, ($N=4$)       |   0.0571 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| Capadenoson 0.001µM, ($N=4$)         | Capadenoson 0.1µM CPA 0.1µM, ($N=4$)   |   0.0286 | *        | $\mathrm{A}_{1}\mathrm{AR}$ |
| Capadenoson 0.001µM, ($N=4$)         | Capadenoson 0.1µM CPA 1µM, ($N=4$)     |   0.1143 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| Capadenoson 0.1µM, ($N=4$)           | Capadenoson 0.1µM CPA 0.1µM, ($N=4$)   |   0.8857 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| Capadenoson 0.1µM, ($N=4$)           | Capadenoson 0.1µM CPA 1µM, ($N=4$)     |   0.1143 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| Capadenoson 1µM, ($N=4$)             | Capadenoson 1µM CPA 0.1µM, ($N=4$)     |   0.4857 | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |
| Capadenoson 1µM, ($N=4$)             | Capadenoson 1µM CPA 1µM, ($N=4$)       |   0.2    | ns       | $\mathrm{A}_{1}\mathrm{AR}$ |

<div align="center">
  <img src="../figures/adpkd-validation-screen/boxplot-MannWhitneyU-plate3-MR_Finerenone_FSK2-5_Aldosterone.svg" alt="" width=720>
</div>

## Finerenone - Simulant (FSK) dose: 2.5µM
| group1                                      | group2                                        |   pvalue | symbol   | target   |
|:--------------------------------------------|:----------------------------------------------|---------:|:---------|:---------|
| FSK 2.5µM, ($N=11$)                         | Aldosterone 0.001µM, ($N=4$)                  |   0.0029 | **       | MR       |
| FSK 2.5µM, ($N=11$)                         | Aldosterone 0.1µM, ($N=4$)                    |   0.0015 | **       | MR       |
| FSK 2.5µM, ($N=11$)                         | Aldosterone 1µM, ($N=4$)                      |   0.0015 | **       | MR       |
| FSK 2.5µM, ($N=11$)                         | Finerenone 0.001µM, ($N=4$)                   |   0.0176 | *        | MR       |
| FSK 2.5µM, ($N=11$)                         | Finerenone 0.1µM, ($N=3$)                     |   0.022  | *        | MR       |
| FSK 2.5µM, ($N=11$)                         | Finerenone 1µM, ($N=4$)                       |   0.1773 | ns       | MR       |
| Finerenone 0.1µM Aldosterone 0.1µM, ($N=4$) | Finerenone 0.001µM Aldosterone 1µM, ($N=4$)   |   0.2    | ns       | MR       |
| Finerenone 0.1µM Aldosterone 0.1µM, ($N=4$) | Finerenone 0.1µM Aldosterone 1µM, ($N=4$)     |   0.1143 | ns       | MR       |
| Finerenone 1µM Aldosterone 0.1µM, ($N=4$)   | Finerenone 1µM Aldosterone 1µM, ($N=4$)       |   0.0286 | *        | MR       |
| Aldosterone 0.1µM, ($N=4$)                  | Finerenone 0.001µM Aldosterone 0.1µM, ($N=4$) |   0.3429 | ns       | MR       |
| Aldosterone 0.1µM, ($N=4$)                  | Finerenone 0.1µM Aldosterone 0.1µM, ($N=4$)   |   0.4857 | ns       | MR       |
| Aldosterone 0.1µM, ($N=4$)                  | Finerenone 1µM Aldosterone 0.1µM, ($N=4$)     |   0.2    | ns       | MR       |
| Aldosterone 1µM, ($N=4$)                    | Finerenone 0.001µM Aldosterone 1µM, ($N=4$)   |   0.4857 | ns       | MR       |
| Aldosterone 1µM, ($N=4$)                    | Finerenone 0.1µM Aldosterone 1µM, ($N=4$)     |   0.6857 | ns       | MR       |
| Aldosterone 1µM, ($N=4$)                    | Finerenone 1µM Aldosterone 1µM, ($N=4$)       |   0.3429 | ns       | MR       |
| Finerenone 0.001µM, ($N=4$)                 | Finerenone 0.1µM Aldosterone 0.1µM, ($N=4$)   |   0.0286 | *        | MR       |
| Finerenone 0.001µM, ($N=4$)                 | Finerenone 0.1µM Aldosterone 1µM, ($N=4$)     |   0.0286 | *        | MR       |
| Finerenone 0.1µM, ($N=3$)                   | Finerenone 0.1µM Aldosterone 0.1µM, ($N=4$)   |   0.0571 | ns       | MR       |
| Finerenone 0.1µM, ($N=3$)                   | Finerenone 0.1µM Aldosterone 1µM, ($N=4$)     |   0.0571 | ns       | MR       |
| Finerenone 1µM, ($N=4$)                     | Finerenone 1µM Aldosterone 0.1µM, ($N=4$)     |   0.0571 | ns       | MR       |
| Finerenone 1µM, ($N=4$)                     | Finerenone 1µM Aldosterone 1µM, ($N=4$)       |   0.0286 | *        | MR       |

<div align="center">
  <img src="../figures/adpkd-validation-screen/boxplot-MannWhitneyU-plate3-MR_Finerenone_FSK0-79_Aldosterone.svg" alt="" width=720>
</div>

## Finerenone - Simulant (FSK) dose: 0.79µM
| group1                                      | group2                                        |   pvalue | symbol   | target   |
|:--------------------------------------------|:----------------------------------------------|---------:|:---------|:---------|
| FSK 0.79µM, ($N=8$)                         | Aldosterone 0.001µM, ($N=4$)                  |   0.6828 | ns       | MR       |
| FSK 0.79µM, ($N=8$)                         | Aldosterone 0.1µM, ($N=4$)                    |   0.004  | **       | MR       |
| FSK 0.79µM, ($N=8$)                         | Aldosterone 1µM, ($N=4$)                      |   0.004  | **       | MR       |
| FSK 0.79µM, ($N=8$)                         | Finerenone 0.001µM, ($N=4$)                   |   0.3677 | ns       | MR       |
| FSK 0.79µM, ($N=8$)                         | Finerenone 0.1µM, ($N=4$)                     |   0.1535 | ns       | MR       |
| FSK 0.79µM, ($N=8$)                         | Finerenone 1µM, ($N=4$)                       |   0.5697 | ns       | MR       |
| Finerenone 0.1µM Aldosterone 0.1µM, ($N=3$) | Finerenone 0.001µM Aldosterone 1µM, ($N=4$)   |   0.0571 | ns       | MR       |
| Finerenone 0.1µM Aldosterone 0.1µM, ($N=3$) | Finerenone 0.1µM Aldosterone 1µM, ($N=4$)     |   0.0571 | ns       | MR       |
| Finerenone 1µM Aldosterone 0.1µM, ($N=4$)   | Finerenone 1µM Aldosterone 1µM, ($N=4$)       |   0.0286 | *        | MR       |
| Aldosterone 0.1µM, ($N=4$)                  | Finerenone 0.001µM Aldosterone 0.1µM, ($N=4$) |   0.2    | ns       | MR       |
| Aldosterone 0.1µM, ($N=4$)                  | Finerenone 0.1µM Aldosterone 0.1µM, ($N=3$)   |   1      | ns       | MR       |
| Aldosterone 0.1µM, ($N=4$)                  | Finerenone 1µM Aldosterone 0.1µM, ($N=4$)     |   0.0571 | ns       | MR       |
| Aldosterone 1µM, ($N=4$)                    | Finerenone 0.001µM Aldosterone 1µM, ($N=4$)   |   0.2    | ns       | MR       |
| Aldosterone 1µM, ($N=4$)                    | Finerenone 0.1µM Aldosterone 1µM, ($N=4$)     |   0.2    | ns       | MR       |
| Aldosterone 1µM, ($N=4$)                    | Finerenone 1µM Aldosterone 1µM, ($N=4$)       |   0.1143 | ns       | MR       |
| Finerenone 0.001µM, ($N=4$)                 | Finerenone 0.1µM Aldosterone 0.1µM, ($N=3$)   |   0.0571 | ns       | MR       |
| Finerenone 0.001µM, ($N=4$)                 | Finerenone 0.1µM Aldosterone 1µM, ($N=4$)     |   0.0286 | *        | MR       |
| Finerenone 0.1µM, ($N=4$)                   | Finerenone 0.1µM Aldosterone 0.1µM, ($N=3$)   |   0.0571 | ns       | MR       |
| Finerenone 0.1µM, ($N=4$)                   | Finerenone 0.1µM Aldosterone 1µM, ($N=4$)     |   0.0286 | *        | MR       |
| Finerenone 1µM, ($N=4$)                     | Finerenone 1µM Aldosterone 0.1µM, ($N=4$)     |   0.0286 | *        | MR       |
| Finerenone 1µM, ($N=4$)                     | Finerenone 1µM Aldosterone 1µM, ($N=4$)       |   0.0286 | *        | MR       |

<div align="center">
  <img src="../figures/adpkd-validation-screen/boxplot-MannWhitneyU-plate3-MR_Esaxerenone_FSK2-5_Aldosterone.svg" alt="" width=720>
</div>

## Esaxerenone - Simulant (FSK) dose: 2.5µM
| group1                                       | group2                                         |   pvalue | symbol   | target   |
|:---------------------------------------------|:-----------------------------------------------|---------:|:---------|:---------|
| FSK 2.5µM, ($N=11$)                          | Aldosterone 0.001µM, ($N=4$)                   |   0.0029 | **       | MR       |
| FSK 2.5µM, ($N=11$)                          | Aldosterone 0.1µM, ($N=4$)                     |   0.0015 | **       | MR       |
| FSK 2.5µM, ($N=11$)                          | Aldosterone 1µM, ($N=4$)                       |   0.0015 | **       | MR       |
| FSK 2.5µM, ($N=11$)                          | Esaxerenone 0.001µM, ($N=4$)                   |   0.7531 | ns       | MR       |
| FSK 2.5µM, ($N=11$)                          | Esaxerenone 0.1µM, ($N=4$)                     |   0.0777 | ns       | MR       |
| FSK 2.5µM, ($N=11$)                          | Esaxerenone 1µM, ($N=4$)                       |   0.0015 | **       | MR       |
| Esaxerenone 0.1µM Aldosterone 0.1µM, ($N=4$) | Esaxerenone 0.001µM Aldosterone 1µM, ($N=4$)   |   0.1143 | ns       | MR       |
| Esaxerenone 0.1µM Aldosterone 0.1µM, ($N=4$) | Esaxerenone 0.1µM Aldosterone 1µM, ($N=4$)     |   0.0286 | *        | MR       |
| Esaxerenone 1µM Aldosterone 0.1µM, ($N=4$)   | Esaxerenone 1µM Aldosterone 1µM, ($N=4$)       |   0.1143 | ns       | MR       |
| Aldosterone 0.1µM, ($N=4$)                   | Esaxerenone 0.001µM Aldosterone 0.1µM, ($N=4$) |   0.1143 | ns       | MR       |
| Aldosterone 0.1µM, ($N=4$)                   | Esaxerenone 0.1µM Aldosterone 0.1µM, ($N=4$)   |   0.4857 | ns       | MR       |
| Aldosterone 0.1µM, ($N=4$)                   | Esaxerenone 1µM Aldosterone 0.1µM, ($N=4$)     |   0.0286 | *        | MR       |
| Aldosterone 1µM, ($N=4$)                     | Esaxerenone 0.001µM Aldosterone 1µM, ($N=4$)   |   0.6857 | ns       | MR       |
| Aldosterone 1µM, ($N=4$)                     | Esaxerenone 0.1µM Aldosterone 1µM, ($N=4$)     |   0.2    | ns       | MR       |
| Aldosterone 1µM, ($N=4$)                     | Esaxerenone 1µM Aldosterone 1µM, ($N=4$)       |   0.0286 | *        | MR       |
| Esaxerenone 0.001µM, ($N=4$)                 | Esaxerenone 0.1µM Aldosterone 0.1µM, ($N=4$)   |   0.0286 | *        | MR       |
| Esaxerenone 0.001µM, ($N=4$)                 | Esaxerenone 0.1µM Aldosterone 1µM, ($N=4$)     |   0.0286 | *        | MR       |
| Esaxerenone 0.1µM, ($N=4$)                   | Esaxerenone 0.1µM Aldosterone 0.1µM, ($N=4$)   |   0.0571 | ns       | MR       |
| Esaxerenone 0.1µM, ($N=4$)                   | Esaxerenone 0.1µM Aldosterone 1µM, ($N=4$)     |   0.0286 | *        | MR       |
| Esaxerenone 1µM, ($N=4$)                     | Esaxerenone 1µM Aldosterone 0.1µM, ($N=4$)     |   0.0286 | *        | MR       |
| Esaxerenone 1µM, ($N=4$)                     | Esaxerenone 1µM Aldosterone 1µM, ($N=4$)       |   0.0286 | *        | MR       |

<div align="center">
  <img src="../figures/adpkd-validation-screen/boxplot-MannWhitneyU-plate3-MR_Esaxerenone_FSK0-79_Aldosterone.svg" alt="" width=720>
</div>

## Esaxerenone - Simulant (FSK) dose: 0.79µM
| group1                                       | group2                                         |   pvalue | symbol   | target   |
|:---------------------------------------------|:-----------------------------------------------|---------:|:---------|:---------|
| FSK 0.79µM, ($N=8$)                          | Aldosterone 0.001µM, ($N=4$)                   |   0.6828 | ns       | MR       |
| FSK 0.79µM, ($N=8$)                          | Aldosterone 0.1µM, ($N=4$)                     |   0.004  | **       | MR       |
| FSK 0.79µM, ($N=8$)                          | Aldosterone 1µM, ($N=4$)                       |   0.004  | **       | MR       |
| FSK 0.79µM, ($N=8$)                          | Esaxerenone 0.001µM, ($N=4$)                   |   0.8081 | ns       | MR       |
| FSK 0.79µM, ($N=8$)                          | Esaxerenone 0.1µM, ($N=4$)                     |   0.0485 | *        | MR       |
| FSK 0.79µM, ($N=8$)                          | Esaxerenone 1µM, ($N=4$)                       |   0.004  | **       | MR       |
| Esaxerenone 0.1µM Aldosterone 0.1µM, ($N=4$) | Esaxerenone 0.001µM Aldosterone 1µM, ($N=4$)   |   0.1143 | ns       | MR       |
| Esaxerenone 0.1µM Aldosterone 0.1µM, ($N=4$) | Esaxerenone 0.1µM Aldosterone 1µM, ($N=4$)     |   0.0286 | *        | MR       |
| Esaxerenone 1µM Aldosterone 0.1µM, ($N=4$)   | Esaxerenone 1µM Aldosterone 1µM, ($N=4$)       |   0.0286 | *        | MR       |
| Aldosterone 0.1µM, ($N=4$)                   | Esaxerenone 0.001µM Aldosterone 0.1µM, ($N=4$) |   0.8857 | ns       | MR       |
| Aldosterone 0.1µM, ($N=4$)                   | Esaxerenone 0.1µM Aldosterone 0.1µM, ($N=4$)   |   0.2    | ns       | MR       |
| Aldosterone 0.1µM, ($N=4$)                   | Esaxerenone 1µM Aldosterone 0.1µM, ($N=4$)     |   0.0286 | *        | MR       |
| Aldosterone 1µM, ($N=4$)                     | Esaxerenone 0.001µM Aldosterone 1µM, ($N=4$)   |   0.8857 | ns       | MR       |
| Aldosterone 1µM, ($N=4$)                     | Esaxerenone 0.1µM Aldosterone 1µM, ($N=4$)     |   0.1143 | ns       | MR       |
| Aldosterone 1µM, ($N=4$)                     | Esaxerenone 1µM Aldosterone 1µM, ($N=4$)       |   0.0286 | *        | MR       |
| Esaxerenone 0.001µM, ($N=4$)                 | Esaxerenone 0.1µM Aldosterone 0.1µM, ($N=4$)   |   0.0286 | *        | MR       |
| Esaxerenone 0.001µM, ($N=4$)                 | Esaxerenone 0.1µM Aldosterone 1µM, ($N=4$)     |   0.0286 | *        | MR       |
| Esaxerenone 0.1µM, ($N=4$)                   | Esaxerenone 0.1µM Aldosterone 0.1µM, ($N=4$)   |   0.0286 | *        | MR       |
| Esaxerenone 0.1µM, ($N=4$)                   | Esaxerenone 0.1µM Aldosterone 1µM, ($N=4$)     |   0.0286 | *        | MR       |
| Esaxerenone 1µM, ($N=4$)                     | Esaxerenone 1µM Aldosterone 0.1µM, ($N=4$)     |   0.0286 | *        | MR       |
| Esaxerenone 1µM, ($N=4$)                     | Esaxerenone 1µM Aldosterone 1µM, ($N=4$)       |   0.0286 | *        | MR       |
