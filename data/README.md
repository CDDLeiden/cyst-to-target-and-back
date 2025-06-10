# Summary

## 1. adpkd_screening

`chemical_structures/`

- [SPECTRUM library](adpkd_screening/chemical_structures/combined_SPECTRUM_structures.sdf) containing the chemical structures of the compounds in the SPECTRUM library.
- [SelleckChem library](adpkd_screening/chemical_structures/sel_chem_structures.tsv) containing the chemical structures of the compounds in the SelleckChem library.

`identified_hits/`

- [Hit Compounds](adpkd_screening/identified_hits/pkd_HitCompounds_NPI-median-DefaultDistance-hitflag_as-isSMILES_20221222-150846.csv) - hits identified per-plate on the NPI-normalized data, according to the `default` distance threshold, with mutual distances of $median(DMSO+FSK) \pm MAD(DMSO+FSK) \times 1.5 $.

`screening_data/`

*SelleckChem library*:
- [Re-analyzed data from Booij et al. 2017](Booij-reanalyzed_selleckchem_Batch2650_2022-11-06.csv.gz) - All non-normalized features extracted from the image analysis workflow. For the manuscript, refer to the original publication: Booij et al. 2017, [DOI: 10.1177/2472555217716056](https://doi.org/10.1177/2472555217716056).

*SPECTRUM library*
- [Re-analyzed data from Booij et al. 2020](Booij-reanalyzed_spectrum_Batch2686_2022-11-07.csv.gz) - All non-normalized features extracted from the image analysis workflow. For the manuscript, refer to the original publication: Booij et al. 2020, [DOI: 10.1093/jmcb/mjz029](https://doi.org/10.1093/jmcb/mjz029).
- Validation Screening [Re-analyzed data from Booij et al. 2020](Booij-reanalyzed_spectrum-validation_Batch2880_2022-11-17.csv.gz) - All non-normalized features extracted from the image analysis workflow. For the manuscript, refer to the original publication: Booij et al. 2020, [DOI: 10.1093/jmcb/mjz029](https://doi.org/10.1093/jmcb/mjz029). This validation screening was performed to validate some of the hits identified in the original screening.

**Note**: We only used readouts from 1µM testing concentration on downstream analysis. The reason for this is that it is was the only concentration used across all the different experiments.

## target_validation

- [Scored compounds](target_validation/Scored_th65_Papyrus_targets_of_ADPKD-screened_compounds.csv) - Contains all the targets identified through the [Papyrus Data Linker notebook](../02_TargetID_and_Prioritization/papyrus_data_linker.ipynb) as well as their $CS_{ratio} = N_{active}/(N_{total})$ reported to all bioactivity types: CS reducer, enhancer, inactive and antineoplastic. The targets are scored according to the following criteria:
- [Validation experiment results](target_validation/ADPKD-TargetValidationScreen_Batch3791_and_Batch3753.csv) - This data is used by the notebook [ADPKD-TargetValidationScreening-analysis](../03_Target_Validation/ADPKD-TargetValidationScreening-analysis.ipynb) to show the results of our validation screening experiment. The data contains the following readouts:

    - obj.Sum(area).meas - the sum of areas occupied by the cystic spheroids in the respective well.
    - obj.Count.meas - the number of cystic spheroids in the respective well.
    - obj.Mean(area).um2.meas - the mean area occupied by the cystic spheroids in the respective well. **Used for hit profiling**.
    - Fraction_dead_cells - a proxy to detect cytotoxicity calculated based on nuclei signal intensity in the absence of a respective spheroid segmentation, indicating apoptosis. Staurosporine (STS) is used as positive control for cytotoxicity.


## compound_exploration

- [Compound exploration experiment results](compound_exploration/-ADPKD-CpdExplorationScreen_Batch4042_and_Batch4064.csv) - Contains the results for compounds selected as part of the `04_Virtual_Screening` section of the main [README.md](../README.md) file. Different from the data reported in `screening_data/`, this file only contains a few readouts, which are:

    - obj.Sum(area).meas - the sum of areas occupied by the cystic spheroids in the respective well.
    - obj.Count.meas - the number of cystic spheroids in the respective well.
    - obj.Mean(area).um2.meas - the mean area occupied by the cystic spheroids in the respective well. **Used for hit profiling**.
    - Fraction_dead_cells - a proxy to detect cytotoxicity calculated based on nuclei signal intensity in the absence of a respective spheroid segmentation, indicating apoptosis. Staurosporine (STS) is used as positive control for cytotoxicity.
