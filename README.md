# cyst-to-target-and-back

Repository with the scripts used for the paper XXXX

Install:
```shell
micromamba env create -n cystToTarget python=3.10 -c conda-forge 
micromamba activate cystToTarget
python -m pip install qsprpred==3.0.2 xgboost==2.1.0
python -m pip install statsmodels venn chembl_webresource_client pubchempy statannotations plotly
python -m pip install papyrus-scripts==1.0.2 papyrus-structure-pipeline==0.0.4
python -m pip install git+https://github.com/valence-labs/mood-experiments.git
python -m pip install -e .
```

# Structure

## Data

All data used and generated in our paper are stored under the `data/` directory.

For information on the available files, please refer to the [data/README.md](data/README.md) file.

## 01_ADPKD_ScreeningProcess

This directoriy contains the scripts used to process the ADPKD screening data coming from (Booij et al. 2017; Booij et al. 2020).

The scripts are:

### 1. screening_process.py

Process the raw screening data, combine with the standardized structures, normalize and save the datasets. The generated files are the following:

- `not_normalized_<dataset_name>.csv` - The dataset with the raw data.
- `npi_median_normalized_<dataset_name>.csv` - The dataset with the raw data normalized by NPI (normalized percent inhibition) as in the [KNIME node](https://nodepit.com/node/de.mpicbg.knime.hcs.base.nodes.norm.npi.NpiNormalizerNodeFactory).
- `z-score_median_normalized_<dataset_name>.csv` - The dataset with the raw data normalized by Z-score.
- `z-prime_median_stats_<dataset_name>.csv` - Z-prime statistics for each feature in the dataset.
 
Additionally, the script will perform hit identification analysis, with figures saved under the [identified_hits](figures/hit_analysis/identified_hits/) directory.

### 2. hit_identification.py

Script to identify the hits in the screening data using the NPI-normalized data produced by script `1`. The generated file is used for section [02_TargetID_and_Prioritization](#02_targetid_and_prioritization). The output file is stored in the [identified_hits](figures/hit_analysis/identified_hits/) directory:

- `pkd_HitCompounds_NPI-median-DefaultDistance-hitflag_as-isSMILES_<date>-<hour>.csv`: hits identified per-plate on the NPI-normalized data, according to the `default` distance threshold, with mutual distances of $median(DMSO+FSK) \pm MAD(DMSO+FSK) \times 1.5 $

## 02_TargetID_and_Prioritization

This directory contains the notebook used to identify the targets of the ADPKD-screened compounds, processed according to [01_ADPKD_ScreeningProcess](#01_adpkd_screeningprocess). A prioritization scheme is also performed for target selection.

### 1. papyrus_data_linker.ipynb

This notebooks performs the target identification and prioritization of the identified targets according to the ADPKD bioactivity. Three different plots are generated:

1. [A Venn diagram](figures/ADPKD-Bioactivity_and_TargetSpace.svg) showing the overlap between targets and the four bioactivity classes (cyst-swelling enhancers, reducers, inactives and antineoplastic);
2. [Bar plot with N(compounds) per bioactivity class](figures/ADPKD-Compound_Bioactivity_Counts-Has_Papyrus_Data.svg);
3. [Bar plot with the prioritized targets](figures/ADPKD-Prioritized_Target_Set-prioritized.svg). On the top part, the number of actives is shown ($N_{active}$), and on the bottom part, the Cyst Swelling (CS) bioactivity ratio: $CS_{ratio} = N_{active}/(N_{total})$, where $N_{total}$ is the total number of ADPKD-screened compounds active for that target on the given threshold (default: $pchembl\_value < 6.5$)

<!-- 
<div align="left">
    <img src="figures/ADPKD-Bioactivity_and_TargetSpace.png" alt="ADPKD Bioactivity and Target Space" width="10%">
</div>
<div align="left">
    <img src="figures/ADPKD-Compound_Bioactivity_Counts-Has_Papyrus_Data.png" alt="ADPKD Bioactivity and Target Space" width="10%">
</div>
<div align="left">
    <img src="figures/ADPKD-Prioritized_Target_Set.png" alt="ADPKD Bioactivity and Target Space" width="10%">
</div> -->

To reproduce the results and analysis, follow the steps in the notebook [papyrus_data_linker.ipynb](02_TargetID_and_Prioritization/papyrus_data_linker.ipynb).

## 03_Target_Validation

This directory contains the notebook used to visualize the results from the ADPKD cystic spheroid target validation experiment. The input data for this notebook is the `data/target_validation/ADPKD-TargetValidationScreen_Batch3791_and_Batch3753.csv` file. The generated figures are the following:

1. Boxplots with Mann-Whitney U test results for the different distributions of the treatment and double treatment groups along the three concentration points (0.001, 0.1, and 1 µM).
2. Point plots with the 95% confidence intervals of the treatment groups along the three concentration points (0.001, 0.1, and 1 µM).

<div align="left">
    <img src="https://www.frontiersin.org/files/Articles/1397864/fphar-15-1397864-HTML-r1/image_m/fphar-15-1397864-g001.jpg" alt="ADPKD Target Validation Boxplots" width="50%">
</div>

## 04_Virtual_Screening

TODO

- Add the data and the notebook for the data analysis on the ADPKD screening data (phenotypical)


<!-- Notes on the data analysis:
- Boxplots are preferred (then they don't "clash" too much with the style of the A1R activity plots);
- Adenosine deaminse and MIPS521:
  - When testing positive allosteric modulators, adenonsine deaminase is usually added to the assay to prevent the accumulation of adenosine;
  - From the gene expression data we got, we can assume there's *some* adenosine deaminase in the system, but not enough to break down all the adenosine present in the system. Therefore, we should acknowledge and not rule out the possibility of A1R activity in the ADPKD screening data;
- Activity observed for 1249141 is *not* due to A1R activity (too low);
- 2707038 ADPKD activity can be explained by A1R activity;
  - At the hit concentration of 1µM we have full receptor ocupancy;
  - Wouldn't call the activity allosteric - it's more like an agonist;
    - On this point, we can still try to plot CPA from the target validation experiment and see if it's an agonist or an allosteric modulator;
- Compound 1237561 shows a lower occupancy for A1R at 1µM;
  - Even though it doesn't show activity for A2A at such concentration, the lower occupancy (~50%) might not be enough for provoke activity;
 -->

## 05_A1R_Screening

TODO