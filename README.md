# cyst-to-target-and-back
Repository with the scripts used for the paper XXXX

Install:
```shell
micromamba env create -n cystToTarget python=3.10 -c conda-forge 
micromamba activate cystToTarget
python -m pip install qsprpred==3.0.2
python -m pip install statsmodels venn chembl_webresource_client pubchempy

python -m pip install papyrus-scripts==1.0.2 papyrus-structure-pipeline==0.0.4
```

# Structure

## Data

```text
data/
├── adpkd_screening/
│   ├── chem_structurs/
│   ├── screening_data/
│   └── standardized_structures/
│
├── anotherpart/
│   ├── dummy_file1.txt
│   └── dummy_file2.txt
```

## 01_ADPKD_ScreeningProcess

This directoriy contains the scripts used to process the ADPKD screening data coming from (Booij et al. 2017; Booij et al. 2020).

The scripts are:

1. `screening_process.py` - Script to process the "raw" screening data, combine with the standardized structures, normalize and save the datasets.
2. `hit_identification.py` - Script to identify the hits in the screening data using the NPI-normalized data produced by script `1`.

