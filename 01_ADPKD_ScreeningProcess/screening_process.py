"""Module to process the data from the ADPKD screening results from the SelleckChem and SPECTRUM libraries.

Input files patterns:
    - SelleckChem: data/adpkd_screening/screening_data/*selleckchem_Batch*.csv
    - SPECTRUM: data/adpkd_screening/screening_data/*SPECTRUM_Batch*.csv
    - SPECTRUM validation: data/adpkd_screening/screening_data/*spectrum-validation_Batch*.csv

Output files:
    # TODO...
"""

import json
import re
import time
from functools import partial
from itertools import combinations
from multiprocessing import Pool
from pathlib import Path
from typing import Union

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import pubchempy as pcp
import seaborn as sns
import statsmodels.api as sm
from feature_manipulation import discrete_transform, get_spearman_corrs, scale_features
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from rdkit import Chem
from rdkit.Chem import DataStructs, PandasTools
from scipy.stats import mannwhitneyu, median_abs_deviation, norm
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

try:
    from chembl_webresource_client.new_client import new_client
except:  # noqa: E722
    pass
from data_prepare import (
    chembl_mol_standardizer,
    chembl_smi_standardizer,
    rd_shut_the_hell_up,
    sanitize_smiles,
    smi_to_fp,
)


def get_mad(x):
    try:
        return median_abs_deviation(x, nan_policy="omit")
    except RuntimeWarning:
        return np.nan


def plate_to_heatmap(
    df: pd.DataFrame,
    plate_n: int,
    measurement: str,
    categorical: bool = False,
    cmap: str = None,
):
    """
    Function to pivot input data into row/col formatting & plotting
    the resulting dataframe as a heatmap.

    Params:
    df -> dataframe with columns `PlateID` identifying plate number,
          `PlateRow` & `PlateColumn` identifying well location.
    plate_n -> `PlateID` to be visualized.
    measurement -> Name of the column containg the desired measurement.
    categorical -> Whether the `measurement` parameter is categorical.
    cmap -> keyword argument with the name of the desired color map.

    For checking available color maps, check:
    https://matplotlib.org/stable/tutorials/colors/colormaps.html
    """
    subset_df = df[df["PlateID"] == plate_n]
    heat_df = subset_df.pivot(index="PlateRow", columns="PlateColumn", values=measurement)

    fig, ax = plt.subplots(figsize=(15, 6))

    if categorical:
        # Creating an interger mapping for the compounds
        value_to_int = {j: i for i, j in enumerate(subset_df[measurement].unique())}
        n = len(value_to_int)

        if cmap is None:
            heatmap = sns.heatmap(
                heat_df.replace(value_to_int),
                annot=True,
                linewidths=0.5,
                ax=ax,
                cmap="tab20b",
            )
        else:
            heatmap = sns.heatmap(
                heat_df.replace(value_to_int),
                annot=True,
                linewidths=0.5,
                ax=ax,
                cmap=cmap,
            )

        # Updating color bar legend
        colorbar = ax.collections[0].colorbar
        r = colorbar.vmax - colorbar.vmin
        colorbar.set_ticks([colorbar.vmin + r / n * (0.5 + i) for i in range(n)])
        colorbar.set_ticklabels(list(value_to_int.keys()))

    else:
        if cmap is None:
            heatmap = sns.heatmap(
                heat_df,
                annot=True,
                linewidths=0.5,
                ax=ax,
                cmap="bwr",
                fmt=".2f",
            )
        else:
            heatmap = sns.heatmap(heat_df, annot=True, linewidths=0.5, ax=ax, cmap=cmap, fmt=".2f")

        ax.set_title(f"Heatmap of plate {1}, measurement = {measurement}")

    fig = heatmap.get_figure()
    return fig, ax


def visualize_distribution(data, xlabel, bins=30):
    """
    Function to visualize the data distribution of a given feature.
    It plots the histogram and the probability density function and
    a qq-plot to check for normality.

    Args:
        data: vector with the data to be plotted.
        xlabel: label for the x-axis.
        bins: number of bins to the histogram. Defaults to 30.

    Returns:
        fig, ax: figure and axis objects.
    """
    fig, axs = plt.subplots(ncols=2, figsize=(12, 4))

    # Fitting a normal distribution to the data:
    (mu, sigma) = norm.fit(data)

    axs[0].hist(data, bins=bins, density=True, alpha=0.6, color="b")
    xmin, xmax = axs[0].get_xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, sigma)

    axs[0].plot(x, p, "k", linewidth=2)
    axs[0].set_title(f"Fit Values: $\mu = {mu:.2f}$ and $\sigma = {sigma:.2f}$")
    axs[0].set_xlabel(xlabel)
    axs[0].set_ylabel("Frequency")

    sm.qqplot(data, line="s", ax=axs[1])
    axs[1].set_title("Q-Q plot of the data")
    return fig, axs


def molname_clean(molname: str) -> str:
    """
    Takes a string and removes anything inside parenthesis. This function is used
    as there are duplicates with the same molecule name but with redundant
    information inside parenthesis. Some names are also customly fixed to allow
    later structure retrieval from pubchempy.

    Args:
        molname: name of the molecule.

    Returns:
        cleaned name of the molecule (without parenthesis).
    """
    pattern = re.compile(r" \(([^\)]+)\)| \[([^\)]+)\]")
    # Some molecule names that I need to fix for pubchempy on get_smiles_from_name()
    mol_dictionary = {
        "DMSO+FSK": "DMSO+FSK",
        "Rapamycin-FSK": "Rapamycin",
        "Roscovitine-FSK": "Roscovitine",
        "NVP-BEZ-235": "Bez-235",
        "Staurosporine": "Staurosporin",
        # "Tofacitinib citrate": "Tofacitinib",
        # "Dovitinib Dilactic acid": "Dovitinib",
        # "Osi-420":
    }
    new_string = pattern.sub("", molname)
    if molname in mol_dictionary.keys():
        new_string = mol_dictionary[molname]

    # if all caps, then lowercase, but not for controls
    if all([new_string.isupper(), new_string not in ["DMSO", "DMSO+FSK"]]):
        new_string = new_string.capitalize()
    return new_string


def get_smiles_from_name(
    comp_name: str, smi_type: str, sleep_time: int = 0, verbose: bool = False
) -> str | bool:
    """
    Gets the list of compounds fetched with PubChemPy for the compound
    name provided. Will return false if no compounds are retrieved. Smiles
    from the resulting list then go through sanitization, depending on the
    `smi_type` argument.

    Args:
        comp_name: name of the molecule to be searched.
        smi_type: type of smiles to return. Either `non-isomeric` or `parent`.
        verbose: turn on/off print statements. Defaults to False.

    Raises:
        AttributeError: if `smi_type` is not `non-isomeric` or `parent`.

    Returns:
        `False` if no compounds are found, otherwise the obtained smiles.
    """

    if smi_type not in ["non-isomeric", "parent"]:
        raise AttributeError("Invalid smi_type attribute. Check documentation " "for the available types.")

    comp_smiles = []
    try:
        if comp_name == "DMSO+FSK":
            comp_list = pcp.get_compounds("Forskolin", "name")
        else:
            comp_list = pcp.get_compounds(comp_name, "name")
            comp_list[0]
    except:
        if verbose:
            print("No compounds in PubChem with this name", comp_name)
        return False

    for c in comp_list:
        if smi_type == "non-isomeric":
            comp_smiles.append(c.canonical_smiles)
        elif smi_type == "parent":
            comp_smiles.append(c.isomeric_smiles)

    # Making a set to have only unique smiles
    candidates = set()
    if smi_type == "non-isomeric":
        for smi in comp_smiles:
            sanit_smi = sanitize_smiles(smi)
            candidates.add(sanit_smi)

    if smi_type == "parent":
        for smi in comp_smiles:
            parent_smi, error = chembl_smi_standardizer(smi)
            if error:  # Don't append the ones that couldn't be standardized
                continue
            candidates.add(parent_smi)
        # chembl smi_standardizer might fail for all molecules...
        if len(candidates) == 0:
            if verbose:
                print(f"chembl_smi_standardizer() failed for all smiles: {comp_name}")
            for smi in comp_smiles:
                parent_smi, error = chembl_smi_standardizer(smi)
                if error:
                    parent_smi = sanitize_smiles(parent_smi)
                candidates.add(parent_smi)

    # From the sanitized molecules, pick the smallest
    final_candidate = min(candidates, key=len)
    time.sleep(sleep_time)
    return final_candidate


def identifier_func():
    return lambda x: x["Compound"] + "~" + x["Concentration"].astype(str) + "~" + x["PlateID"].astype(str)


class selleck_chem:
    """
    Class for processing and visualizing the data
    from the SelleckChem screening.

    Standard identification columns are:
    `Compound`, `Concentration`, `PlateID`, `PlateColumn`, `PlateRow`
    """

    def __init__(self, name: str = "selleckchem") -> None:
        self.name = name
        self.root_dir = Path(__file__).absolute().parents[1]
        self.file_root = self.root_dir / "data/adpkd_screening/screening_data"
        self.file_path = list(self.file_root.glob(f"*{self.name}_Batch*.csv.gz"))[0]
        print(f"    Loading {self.name} screening under: {self.file_path}")
        self.chemstructs = None
        self.chemstructs_path = (
            self.root_dir / "data/adpkd_screening/chemical_structures/sel_chem_structures.tsv"
        )
        self.dropped_comps = dict()
        self.comp_mapping = dict()
        self.df = pd.read_csv(self.file_path)
        self.norm_df = None  # Holds normalized dataframe
        self.normalization_type = None  # Normalization method -> NPI / Z'
        self.normalization_method = None  # Whether used mean or median
        self.stats_df = None  # Holds dataframe after statistics
        self.stats_type = None  # type of statistical testing
        self.stats_method = None  # Whether used mean or median
        self.stats_test_type = None  # Type of statistical test for p-values
        self.multitest_results = None  # Corrected p-values
        self.control_treatments = [
            "bez-235",
            "rapamycin",
            "roscovitine",
            "sorafenib",
            "metformin",
        ]
        self.toxic_treatments = ["staurosporin"]
        self.sd_val_cols = None
        self.desired_cols = None
        self.id_cols = [
            "TreatmentType",
            "Concentration",
            "PlateID",
            "PlateColumn",
            "PlateRow",
            "QC",
        ]
        self.comp_info_cols = [
            "SMILES",
            "Antineoplastic",
        ]
        self.outlier_wells = [
            "A1_7",
            "A2_8",
            "K2_2",
            "E14_7",
            "H15_7",
            "H15_8",
            "H16_7",
            "H16_8",
        ]

    def drop_cols(self) -> None:
        """
        Drops undersired columns from self.df
        """
        remove_pattern = re.compile(r"path|name|\.tif|\.csv|code|folder|^row$|root")
        to_drop_cols = [c for c in self.df.columns if remove_pattern.findall(c)]
        self.df.drop(columns=to_drop_cols, inplace=True)

        rename_pattern = re.compile(r"plate\.layout\.info\.")
        to_rename = {  # Removing this prefix as we don't need it
            key: value
            for key, value in zip(
                self.df.columns,
                [rename_pattern.sub("", c) for c in self.df.columns],
            )
        }
        self.df.rename(columns=to_rename, inplace=True)
        return

    def update_stimulation(self) -> None:
        """
        Funtion to update compound name based on stimulation values.
        Negative control -> DMSO; positive control -> DMSO+FSK
        """
        comp_stim = self.df["Compound"] + "+" + self.df["Stimulation"]

        def keep_only(string):
            if string == "DMSO+FSK":
                return string
            elif string == "DMSO+DMSO":
                return "DMSO"
            else:
                return string.split("+")[0]

        updated_comp = comp_stim.apply(keep_only)
        self.df = self.df.assign(Compound=updated_comp)
        return

    def update_id_attributes(self) -> None:
        """
        1) Updates the indentification column names from self.df for
        further normalization steps.
        2) Sorts the dataframe by the identification columns
        PlateID, PlateColumn, PlateRow
        """

        def plateID_to_number(plate_id_str):
            splitted = plate_id_str.split("_")
            return int(splitted[1])

        # Standardizing compound names with molname_clean
        standard_names = {name: molname_clean(name) for name in self.df["Compound"].unique()}
        self.df.replace({"Compound": standard_names}, inplace=True)

        # mapping the types of treatment within the screening
        numerics = ["int16", "int32", "int64", "float16", "float32", "float64"]
        conditions = [
            (self.df["Compound"] == "DMSO"),
            (self.df["Compound"] == "DMSO+FSK"),
            (self.df["Compound"].str.lower().isin(self.control_treatments)),
            (
                (self.df["Compound"].str.lower().isin(self.toxic_treatments))
                # & (self.df["Concentration_micromolar"] >= 1)
            ),
        ]
        choices = ["neg_control", "pos_control", "treatment_control", "treatment_toxic"]
        self.df["TreatmentType"] = np.select(conditions, choices, default="treatment")

        newdf = (
            self.df.select_dtypes(include=numerics)
            .assign(
                Compound=self.df["Compound"],
                TreatmentType=self.df["TreatmentType"],
                PlateID=self.df["Plate_layout"],
                PlateRow=self.df["row.char"],
                PlateColumn=self.df["column"],
                QC=self.df["QC"],
            )
            .rename(
                columns={
                    "Concentration_micromolar": "Concentration",
                }
            )
            .drop(
                columns=["column"],
                # `column` dropped as it's numerical and remains after select_dtypes
            )
            .replace("NO CELLS", "Empty")
        )

        # Fixing column order, changing PlateID to numerical and resetting index
        self.df = (
            newdf[["Compound"] + self.id_cols + newdf.columns.tolist()[1:-6]]
            .assign(PlateID=newdf["PlateID"].apply(plateID_to_number))
            .sort_values(by=["PlateID", "PlateRow", "PlateColumn"])
            .reset_index(drop=True)
        )
        return

    def update_features(self) -> None:
        """
        1) Drops features that are not used for our analysis & renames
        them according to regular expressions.
        2) Checks for wells with nan values in the desired features and
        flags them on the QC column if `QC == 'OK'`.
        """
        # Part 1)
        columns = self.df.columns
        # Dropping features in micrometer as I won't use those:
        um_pattern = re.compile("\.um")
        rm_features_pattern = re.compile("^Fraction")  # |num\.obj")
        # Fraction_dead_cells removed as it it doesn't translate well in our images
        # num.obj removed as it may tend to infinity after normalization
        todrop_cols = [c for c in columns if any([um_pattern.findall(c), rm_features_pattern.findall(c)])]
        self.df.drop(columns=todrop_cols, inplace=True)

        remove_meas = re.compile("\.meas|\*")  # There's a feature with Count*(area)
        rename_sd = re.compile("Standard\sdeviation")  # Rename for abbreviation
        zorder = re.compile("zernike_order")  # Rename for abbreviation

        update_features = dict()
        for c in self.df.columns:
            oldstring = c
            if remove_meas.findall(c):
                c = re.sub(pattern=remove_meas, repl="", string=c)
            if rename_sd.findall(c):
                c = re.sub(pattern=rename_sd, repl="SD", string=c)
            if zorder.findall(c):
                c = re.sub(pattern=zorder, repl="ze_order", string=c)
            update_features[oldstring] = c

        # Renaming columns according to the updated feature names
        self.df.rename(columns=update_features, inplace=True)

        # Saving value columns as object attributes for later use
        dev_pattern = re.compile("\.SD")
        val_cols = [c for c in self.df.columns if c not in ["Compound"] + self.id_cols]

        self.sd_val_cols = [
            c for c in val_cols if dev_pattern.findall(c)
        ]  # standard deviations, not used in further analysis..
        self.desired_cols = [c for c in val_cols if not dev_pattern.findall(c)]  # mean values

        # Printing information about the wells and doing QC flagging
        empty_flagged = [
            "empty" if comp == "Empty" else qc for comp, qc in zip(self.df["Compound"], self.df["QC"])
        ]
        self.df = self.df.assign(QC=empty_flagged)

        nan_idxs = self.df[self.desired_cols][(self.df.isna().any(axis=1))].index
        not_ok_idxs = self.df[self.desired_cols][self.df["QC"] != "OK"].index
        updated_QC = [
            "has_nan" if all([idx in nan_idxs, idx not in not_ok_idxs]) else qc
            for idx, qc in enumerate(self.df["QC"])
        ]
        print(f"Detected {len(nan_idxs)} wells with nan values. Flagged as 'has_nan'.")
        print(f"Detected {len(not_ok_idxs)} wells with QC flagging as not 'OK'.")
        self.df = self.df.assign(QC=updated_QC)
        return todrop_cols

    def flag_outliers(
        self,
    ):
        """
        Function to flag outliers on the QC. This was based on the analysis
        I've done on 20221212 (see notes on vljournal). Wrote this function
        to make it complatible with the hit_identification.py module.
        """
        if len(self.outlier_wells) == 0:
            return
        wellid_lambda_func = (
            lambda x: x["PlateRow"] + x["PlateColumn"].astype(str) + "_" + x["PlateID"].astype(str)
        )

        df = self.df.copy()
        df = df.assign(PlateWellID=wellid_lambda_func)
        df_idxs = df.index

        condition = np.array([df["PlateWellID"].isin(self.outlier_wells)]).ravel()
        outlier_idxs = np.compress(condition, df_idxs)

        qc_flag = ["outlier" if idx in outlier_idxs else qc for qc, idx in zip(df["QC"], df_idxs)]
        self.df = self.df.assign(QC=qc_flag)

    def npi_normalize(
        self,
        method: str = "median",
        stack_norm: bool = False,
        dropna: bool = False,
    ) -> pd.DataFrame:
        """Function inspired on the KNIME node Normalize Plates (NPI).
        https://nodepit.com/node/de.mpicbg.knime.hcs.base.nodes.norm.npi.NpiNormalizerNodeFactory

        Original formula is set for antagonist assays, following:
        >>> x = (mean(x[sub{pos}]) - x) / (mean(x[sub{pos}]) - mean(x[sub{neg}])) * 100

        As our essay is not antagonistic, we implemented instead:
        >>> x = (mean(x[sub{neg}]) - x) / (mean(x[sub{neg}]) - mean(x[sub{pos}])) * 100

        Args:
            method: whether to perform the normalization using the mean
                or the median. Defaults to "median".
            stack_norm: if desired to stack the normaliztion on top of `self.norm_df`.
                Defaults to False.
            dropna: whether to drop nan values. Defaults to False.

        Raises:
            AttributeError: if method is not "median" or "mean".

        Returns:
            norm_df(pd.DataFrame) -> normalized dataframe.
        """
        if method not in ["median", "mean"]:
            raise AttributeError("method unavailable. Should be either `mean` or `median`")
        start_len = len(self.df)
        if stack_norm:  # applied to self.norm_df
            norm_df = pd.DataFrame(columns=self.norm_df.columns)
            if dropna:
                grouped_plates = self.norm_df.dropna().groupby("PlateID")
            else:
                grouped_plates = self.norm_df.groupby("PlateID")
        else:
            norm_df = pd.DataFrame(columns=self.df.columns)
            if dropna:
                grouped_plates = self.df.dropna().groupby("PlateID")
            else:
                grouped_plates = self.df.groupby("PlateID")

        # Method performed per plate:
        for pnumber, plate_df in grouped_plates:
            # Saving values that will be dropped
            t_type = plate_df["TreatmentType"].values
            conc = plate_df["Concentration"].values
            plate = plate_df["PlateID"].values
            col = plate_df["PlateColumn"].values
            row = plate_df["PlateRow"].values
            qc = plate_df["QC"].values

            # Dropping to have only numerical values to normalize
            plate_df.drop(columns=self.id_cols, inplace=True)

            if method == "mean":
                neg_control_m = plate_df.groupby("Compound").mean().loc["DMSO"]
                pos_control_m = plate_df.groupby("Compound").mean().loc["DMSO+FSK"]

            elif method == "median":
                neg_control_m = plate_df.groupby("Compound").median().loc["DMSO"]
                pos_control_m = plate_df.groupby("Compound").median().loc["DMSO+FSK"]

            # Function formula:
            # x = (mean(x[sub{pos}]) - x) / (mean(x[sub{pos}]) - mean(x[sub{neg}])) * 100
            denominator = neg_control_m - pos_control_m
            # Switched around so that we have percentage of activation instead of inhibition
            plate_df = (
                plate_df.set_index("Compound")
                .apply(lambda x: x * -1)
                .add(neg_control_m, axis="columns")
                .div(denominator, axis="columns")
                .apply(lambda x: x * 100)
                .reset_index()
                .assign(
                    TreatmentType=t_type,
                    Concentration=conc,
                    PlateID=plate,
                    PlateColumn=col,
                    PlateRow=row,
                    QC=qc,
                )
            )

            # Appending normalized data to new pd.DataFrame
            norm_df = pd.concat([norm_df, plate_df])

        # Fixing column order & resetting index
        norm_df = norm_df[["Compound"] + self.id_cols + plate_df.columns.tolist()[1:-6]].reset_index(
            drop=True
        )
        if dropna:
            final_len = len(norm_df)
            print(f"Dropped {start_len - final_len} rows with nan values.")

        self.norm_df = norm_df
        self.normalization_type = "npi"
        self.normalization_method = method
        return norm_df

    def minmax_normalize(
        self,
        stack_norm: bool = False,
        dropna: bool = False,
    ) -> pd.DataFrame:
        """
        Function to perform the minmax normalization on the dataframe.

        Args:
            stack_norm: if desired to stack the normaliztion on top of `self.norm_df`.
                Defaults to False.
            dropna: whether to drop nan values. Defaults to False.

        Returns:
            norm_df(pd.DataFrame) -> normalized dataframe.
        """
        start_len = len(self.df)
        if stack_norm:  # applied to self.norm_df
            norm_df = pd.DataFrame(columns=self.norm_df.columns)
            if dropna:
                grouped_plates = self.norm_df.dropna().groupby("PlateID")
            else:
                grouped_plates = self.norm_df.groupby("PlateID")
        else:
            norm_df = pd.DataFrame(columns=self.df.columns)
            if dropna:
                grouped_plates = self.df.dropna().groupby("PlateID")
            else:
                grouped_plates = self.df.groupby("PlateID")

        for pnumber, plate_df in grouped_plates:
            # Saving values that will be dropped
            t_type = plate_df["TreatmentType"].values
            conc = plate_df["Concentration"].values
            plate = plate_df["PlateID"].values
            col = plate_df["PlateColumn"].values
            row = plate_df["PlateRow"].values
            qc = plate_df["QC"].values

            # Dropping to have only numerical values to normalize
            plate_df.drop(columns=self.id_cols, inplace=True)
            plate_df = (
                plate_df.set_index("Compound")
                .apply(scale_features, axis=1)
                .reset_index()
                .assign(
                    TreatmentType=t_type,
                    Concentration=conc,
                    PlateID=plate,
                    PlateColumn=col,
                    PlateRow=row,
                    QC=qc,
                )
            )
            # Appending normalized data to new pd.DataFrame
            norm_df = pd.concat([norm_df, plate_df])

        # Fixing column order & resetting index
        norm_df = norm_df[["Compound"] + self.id_cols + plate_df.columns.tolist()[1:-6]].reset_index(
            drop=True
        )
        if dropna:
            final_len = len(norm_df.dropna())
            print(f"Dropped {start_len - final_len} rows with nan values.")

        self.norm_df = norm_df
        self.normalization_type = "minmax"
        return norm_df

    def z_score_normalize(
        self,
        method: str = "median",
        dropna: bool = True,
        as_booij: bool = True,
    ) -> pd.DataFrame:
        """
        Function to perform the Z-score normalization of self.df dataframe. The default
        implementation of the function follows the description given by T.
        Booij et al.
        in their paper "High-Throughput Phenotypic Screening of Kinase Inhibitors to
        Identify Drug Targets for Polycystic Kidney Disease".

        However, the standard implementation of the Z-score normalization is as follows:
        Z = (x - mean(x)) / std(x). This normalization can be performed by setting the
        `as_booij` parameter to False.

        Args:
            method: whether to perform the normalization using
                  the mean or the median (for robust z-score).
                  Defaults to "median".
            dropna: whether to drop nan values. Defaults to False.
            as_booij: implementation of z-score as described in the paper by Booij et al.
            If false, will use default implementation. Defaults to True.

        Raises:
            AttributeError: if method is not "median" or "mean".

        Returns:
            norm_df(pd.DataFrame) -> normalized dataframe.
        """

        numerics = ["int16", "int32", "int64", "float16", "float32", "float64"]

        if method not in ["median", "mean"]:
            raise AttributeError("method unavailable. Should be either `mean` or `median`")

        norm_df = pd.DataFrame(columns=self.df.columns)

        # z-score normalization should be done per plate:
        if dropna:
            start_len = len(self.df)
            grouped_plates = self.df.dropna().groupby("PlateID")
            final_len = len(self.df.dropna())
            print(f"Dropped {start_len - final_len} rows with nan values.")
        else:
            grouped_plates = self.df.groupby("PlateID")
        # Performing the z-score normalizaton per plate.
        for pnumber, plate_df in grouped_plates:
            # Saving values that will be dropped
            t_type = plate_df["TreatmentType"].values
            conc = plate_df["Concentration"].values
            plate = plate_df["PlateID"].values
            col = plate_df["PlateColumn"].values
            row = plate_df["PlateRow"].values
            qc = plate_df["QC"].values

            # Dropping to have only numerical values to normalize
            plate_df.drop(columns=self.id_cols, inplace=True)

            if as_booij:
                if method == "mean":
                    neg_control_m = plate_df.groupby("Compound").mean().loc["DMSO"]
                    neg_control_dev = plate_df.groupby("Compound").std().loc["DMSO"]

                elif method == "median":
                    neg_control_m = plate_df.groupby("Compound").median().loc["DMSO"]
                    neg_control_dev = (
                        plate_df.query("Compound == 'DMSO'")
                        .select_dtypes(numerics)
                        .apply(median_abs_deviation, nan_policy="omit")
                    )

                assert (neg_control_m.index == neg_control_dev.index).all()

                plate_df = (
                    plate_df.set_index("Compound")
                    .sub(neg_control_m, axis="columns")
                    .div(neg_control_dev, axis="columns")
                    .reset_index()
                    .assign(
                        TreatmentType=t_type,
                        Concentration=conc,
                        PlateID=plate,
                        PlateColumn=col,
                        PlateRow=row,
                        QC=qc,
                    )
                )
            else:
                if method == "mean":
                    control_m = plate_df.mean()
                    control_dev = plate_df.std()

                elif method == "median":
                    control_m = plate_df.groupby("Compound").median().loc["DMSO"]
                    control_dev = plate_df.select_dtypes(numerics).apply(
                        median_abs_deviation, nan_policy="omit"
                    )

                assert (control_m.index == control_dev.index).all()

                plate_df = (
                    plate_df.set_index("Compound")
                    .sub(control_m, axis="columns")
                    .div(control_dev, axis="columns")
                    .reset_index()
                    .assign(
                        TreatmentType=t_type,
                        Concentration=conc,
                        PlateID=plate,
                        PlateColumn=col,
                        PlateRow=row,
                        QC=qc,
                    )
                )

            # Appending normalized data to new pd.DataFrame
            norm_df = pd.concat([norm_df, plate_df])

        # Fixing column order & resetting index
        norm_df = norm_df[["Compound"] + self.id_cols + plate_df.columns.tolist()[1:-6]].reset_index(
            drop=True
        )
        self.norm_df = norm_df
        if as_booij:
            self.normalization_type = "z-score_booij"
        else:
            self.normalization_type = "z-score"
        self.normalization_method = method
        return norm_df

    def b_score_normalize(
        self, max_iter: int, tol: float, which: str = "df", hide_progress: bool = True
    ) -> pd.DataFrame:
        """
        Function to perform the B-score normalization of the feature
        values. This normalization is performed using the iterative
        algorithm described in the paper "Improved Statistical Methods for
        Hit Selection in High-Throughput Screening" by Brideau et al. (2003).

        For the calculation of the score, a two-way fitted median polish
        is applied to the values. The obtained residuals are then divided
        by the median absolute deviation of the plate.

        The following implementation was inspired on the code from:
        https://github.com/borisvish/Median-Polish/blob/master/AdditiveModelFitByMedianPolish.py
        https://sparkrma.readthedocs.io/en/latest/_modules/spark_rma/median_polish.html


        Args:
            max_iter: maximun number of median polish iterations.
            tol: tolerance value for the median polish algorithm.
            which: whether to apply it to `self.df` or `self.norm_df`.
            Defaults to "df".
            hide_progress: hides the progress bar. Defaults to True.

        Raises:
            AttributeError: value `which` is not "df" or "norm_df".

        Returns:
            plate_dict(dict), new_df(pd.Dataframe), containing the b-scored
            feature values.
        """
        if which not in ["norm_df", "df"]:
            raise AttributeError("Attribute unavailable. Should be either `norm_df` or `df`")

        if which == "norm_df":
            df = self.norm_df.copy()
        elif which == "df":
            df = self.df.copy()

        # Need to think about what to do with the nans
        df = df.fillna(0).sort_values(by=["PlateID", "PlateRow", "PlateColumn"]).reset_index(drop=True)
        vals = self.desired_cols

        def update_row(series, grouped_median):
            subset = grouped_median.loc[series["PlateRow"], :]
            return series[vals] - subset[vals]

        def update_col(series, grouped_median):
            subset = grouped_median.loc[series["PlateColumn"], :]
            return series[vals] - subset[vals]

        plate_dict = {p: {} for p in df["PlateID"].unique()}

        # Labels for initializing the median polish
        row_labels = sorted(df["PlateRow"].unique())
        col_labels = sorted(df["PlateColumn"].unique())

        for pnumber in df["PlateID"].unique():
            plate_df = df.query(f"PlateID == {pnumber}").copy()

            grand_effect = pd.Series(0, index=vals, dtype=np.float64)
            # Initializing median row effect -> one value per feature!
            median_row_effects = pd.Series(0, index=vals, dtype=np.float64)
            # Initializing median column effect -> one value per feature!
            median_col_effects = pd.Series(0, index=vals, dtype=np.float64)

            # Initializing row effect -> one value per row
            row_effects = pd.DataFrame(0, columns=vals, index=row_labels, dtype=np.float64)
            # Initializing column effect -> one value per column
            col_effects = pd.DataFrame(0, columns=vals, index=col_labels, dtype=np.float64)

            for i in tqdm(range(0, max_iter), disable=hide_progress):
                row_medians = plate_df.loc[:, vals + ["PlateRow"]].groupby("PlateRow").median()[vals]
                # return row_effects, row_medians
                row_effects += row_medians
                median_row_effects = row_effects.median()

                grand_effect += median_row_effects
                row_effects -= median_row_effects

                # Perform the median scrape on each row and update the residuals
                plate_df.loc[:, vals] = plate_df.apply(
                    partial(update_row, grouped_median=row_medians), axis=1
                )

                col_medians = plate_df.loc[:, vals + ["PlateColumn"]].groupby("PlateColumn").median()[vals]
                col_effects += col_medians
                median_col_effects = col_effects.median()

                # Perform the median scrape on each column and update the residuals
                plate_df.loc[:, vals] = plate_df.apply(
                    partial(update_col, grouped_median=col_medians), axis=1
                )
                grand_effect += median_col_effects

                # if any of the values in the row/col effects are less than the tolerance, stop
                conditions = [
                    row_effects.abs().max().any() < tol,
                    col_effects.abs().max().any() < tol,
                ]
                if any(conditions):
                    print(f"Plate {pnumber} converged after {i} iterations")
                    break

            # Calculating the median absolute deviation of each feature within the plate
            feature_mads = plate_df[vals].apply(median_abs_deviation, axis=0)
            b_scores = plate_df[vals].div(feature_mads)
            plate_dict[pnumber].update(
                {
                    "grand_effect": grand_effect,
                    "row_effects": row_effects,
                    "col_effects": col_effects,
                    "residuals": plate_df[vals],
                    "original": df.query(f"PlateID == {pnumber}"),
                    "feature_mads": feature_mads,
                    "b_score": b_scores,
                }
            )

        new_df = pd.DataFrame(columns=["Compound"] + self.id_cols + self.desired_cols)
        for n in range(1, len(plate_dict) + 1):
            new_plate_df = pd.concat(
                [
                    plate_dict[n]["original"].loc[:, ["Compound"] + self.id_cols],
                    plate_dict[n]["b_score"],
                ],
                axis=1,
            )
            new_df = pd.concat([new_df, new_plate_df], axis=0)
        new_df.reset_index(drop=True, inplace=True)
        self.norm_df = new_df
        self.normalization_type = "b-score"
        self.normalization_method = "median-polish"
        return plate_dict, new_df

    def save_data(self, which: str) -> None:
        """
        Function to save the desired dataframe

        Params:
        which -> which dataframe to save. Either `df` for the
                 not normalized, `norm_df` for the normalized data,
                 or `stats_df` for the statistics data.
        """

        if which not in ["df", "norm_df", "stats_df"]:
            raise AttributeError("Which parameter should be either `df`, `norm_df`, `stats_df`.")
        if which == "df":
            pd.DataFrame.to_csv(
                self.df,
                self.file_root / f"not_normalized_{self.name}.csv",
                index=False,
            )

        elif which == "norm_df":
            id_string = f"{self.normalization_type}_{self.normalization_method}"
            pd.DataFrame.to_csv(
                self.norm_df,
                self.file_root / f"{id_string}_normalized_{self.name}.csv",
                index=False,
            )

        elif which == "stats_df":
            id_string = f"{self.stats_type}_{self.stats_method}"
            pd.DataFrame.to_csv(
                self.stats_df,
                self.file_root / f"{id_string}_stats_{self.name}.csv",
                index=False,
            )

    def plot_data_distribution(
        self,
        value_var: str,
        plot_type: str,
        per_plate: bool,
        which: str = "norm_df",
    ):
        """
        Method to plot the self.dataframe `value_var` sample distributions.

        Params:
        value_var -> column name of the feature to be plotted.
        plot_type -> Available: `violinplot`, `boxplot`.
        per_plate -> whether to plot the distribution per plate or not.
        which -> Either `df` for the not normalized, `norm_df` for the normalized data.
        """

        if which not in ["df", "norm_df"]:
            raise AttributeError("Which parameter should be either `df`, or `norm_df`.")

        if plot_type not in ["violinplot", "boxplot"]:
            raise AttributeError("plot_type not available.")

        # Subsetting desired bioactivity
        if which == "df":
            df = self.df
        elif which == "norm_df":
            df = self.norm_df

        controls_only = df[df["TreatmentType"].isin(["neg_control", "pos_control"])]
        id_vars = self.id_cols + ["Compound"]
        plot_df = controls_only.melt(id_vars=id_vars, value_vars=value_var)

        if per_plate:
            fig, ax = plt.subplots(figsize=(5, 3))
            ax.set_title(f"{value_var} Distribution (per plate)")
        else:
            fig, ax = plt.subplots(figsize=(3, 3))
            ax.set_title(f"{value_var} Distribution (aggregated)")

        if plot_type == "violinplot":
            figure = sns.violinplot(data=plot_df, x="PlateID", y="value", hue="Compound", ax=ax)
        elif plot_type == "boxplot":
            figure = sns.boxplot(data=plot_df, x="PlateID", y="value", hue="Compound", ax=ax)

        ax.legend(bbox_to_anchor=(1.04, 0), loc="lower left", borderaxespad=0)
        fig = figure.get_figure()
        return fig, ax

    def compare_normalizations(
        self,
        value_var: str,
        plot_type: str,
        query_col: str = "Compound",
        query_items: list = ["DMSO", "DMSO+FSK"],
        per_plate: bool = True,
        use_saved: bool = True,
        figsize=None,
    ):
        """
        Method to plot sample distributions after applying the different
        normalization methods.

        Params:
        value_var -> column name of the feature to be plotted.
        plot_type -> Available: `violinplot`, `boxplot`.
        query_col -> column name to query the dataframe.
        query_items -> list of items to query the dataframe (from query_col).
        per_pate -> whether to plot the whole dataframe or not.
        use_saved -> whether the dataframes are already saved.
        figsize -> tuple of the figure size if desired.

        For average size of the cystic objects, use `value_var=obj.Mean(area)`.
        """

        if plot_type not in ["violinplot", "boxplot"]:
            raise AttributeError("plot_type not available.")

        if use_saved:
            npi_df = pd.read_csv(list(self.file_root.glob(f"npi*{self.name}.csv"))[0])
            z_df = pd.read_csv(list(self.file_root.glob(f"z-score*{self.name}.csv"))[0])

        else:
            # Calling and initializing the class to allow inheritance from `spectrum` class
            screening = self.__class__()
            screening.__init__()
            screening.drop_cols()
            screening.update_stimulation()
            screening.update_id_attributes()
            screening.npi_normalize()
            npi_df = screening.norm_df
            screening.z_score_normalize()
            z_df = screening.norm_df

        if per_plate:
            if figsize is None:
                fig, axs = plt.subplots(ncols=2, figsize=(10, 3))
            else:
                fig, axs = plt.subplots(ncols=2, figsize=(figsize[0], figsize[1]))
        else:
            if figsize is None:
                fig, axs = plt.subplots(ncols=2, figsize=(6, 3))
            else:
                fig, axs = plt.subplots(ncols=2, figsize=(figsize[0], figsize[1]))

        for ax, df in zip(axs, [z_df, npi_df]):
            # Subsetting desired bioactivity
            controls_only = df[df[query_col].isin(query_items)]
            id_vars = self.id_cols + ["Compound"]
            plot_df = controls_only.melt(id_vars=id_vars, value_vars=value_var)

            if per_plate:
                if plot_type == "violinplot":
                    sns.violinplot(
                        data=plot_df,
                        x="PlateID",
                        y="value",
                        hue=query_col,
                        hue_order=query_items,
                        ax=ax,
                    )
                elif plot_type == "boxplot":
                    sns.boxplot(
                        data=plot_df,
                        x="PlateID",
                        y="value",
                        hue=query_col,
                        hue_order=query_items,
                        ax=ax,
                    )
                ax.legend(bbox_to_anchor=(1.04, 0), loc="lower left", borderaxespad=0)
                # Setting plot title
                if df is z_df:
                    ax.set_title("Z-Score normalization (per plate)")
                else:
                    ax.set_title("NPI normalization (per plate)")

            else:  # Aggregated values
                if plot_type == "violinplot":
                    sns.violinplot(data=plot_df, x=query_col, y="value", ax=ax)
                elif plot_type == "boxplot":
                    sns.boxplot(data=plot_df, x=query_col, y="value", ax=ax)
                # Setting plot title
                if df is z_df:
                    ax.set_title("Z-Score normalization (aggregated)")
                else:
                    ax.set_title("NPI normalization (aggregated)")

        fig.tight_layout()
        return fig, axs

    def z_prime_statistics(self, value_var: str, method: str = "median") -> pd.DataFrame:
        """Function to calculate the z'-factor for each plate.

        Args:
            value_var: column name for calculating the z-prime statistics.
            Use `all` for calculating z-scores of all features.
            method: whether to use mean (z' calculation) or median (robust z').
            Defaults to "median".

        Raises:
            AttributeError: when method is not "median" or "mean".

        Returns:
            stats_df -> dataframe with the z-prime statistics.
        """
        numerics = ["int16", "int32", "int64", "float16", "float32", "float64"]

        if method not in ["median", "mean"]:
            raise AttributeError("method unavailable. Should be either `mean` or `median`")

        norm_df = self.norm_df.drop(columns=["Concentration"])
        stats_df = pd.DataFrame(columns=["Measurement", "Z-prime"])

        grouped_plates = norm_df.groupby("PlateID")
        for pnumber, plate_df in grouped_plates:
            if method == "median":
                neg_control_m = (
                    plate_df.set_index("Compound")
                    .select_dtypes(numerics)
                    .reset_index(names=["Compound"])
                    .groupby("Compound")
                    .median()
                    .loc["DMSO"]
                )
                neg_control_dev = plate_df.query("Compound == 'DMSO'").select_dtypes(numerics).apply(get_mad)
                assert (neg_control_m.index == neg_control_dev.index).all()

                pos_control_m = (
                    plate_df.set_index("Compound")
                    .select_dtypes(numerics)
                    .reset_index(names=["Compound"])
                    .groupby("Compound")
                    .median()
                    .loc["DMSO+FSK"]
                )
                pos_control_dev = (
                    plate_df.query("Compound == 'DMSO+FSK'").select_dtypes(numerics).apply(get_mad)
                )
                assert (pos_control_m.index == pos_control_dev.index).all()

            if method == "mean":
                neg_control_m = plate_df.groupby("Compound").mean().loc["DMSO"]
                neg_control_dev = plate_df.groupby("Compound").std().loc["DMSO"]
                pos_control_m = plate_df.groupby("Compound").mean().loc["DMSO+FSK"]
                pos_control_dev = plate_df.groupby("Compound").std().loc["DMSO+FSK"]

            denominator = pos_control_m.sub(neg_control_m).abs()
            z_prime = pos_control_dev.add(neg_control_dev).apply(lambda x: x * 3).div(denominator)
            z_prime = 1 - z_prime

            # If to get all z-scores
            if value_var == "all":
                measurements = [name for name in z_prime.index]
                values = [z_prime[name] for name in z_prime.index]
                data = pd.DataFrame(
                    {"Measurement": measurements, "Z-prime": values},
                    index=[pnumber] * len(z_prime),
                )
                stats_df = pd.concat([stats_df, data], axis=0)

            else:
                data = pd.DataFrame(
                    {"Measurement": value_var, "Z-prime": z_prime[value_var]},
                    index=[pnumber],
                )
                stats_df = pd.concat([stats_df, data], axis=0)

        # Updating index and ordering the dataframe
        stats_df = stats_df.reset_index().rename(columns={"index": "PlateID"})
        order = ["PlateID", "Measurement", "Z-prime"]
        stats_df = stats_df[order]

        if value_var == "all":
            stats_df = stats_df.sort_values(by=["Measurement", "PlateID"], ascending=True).reset_index(
                drop=True
            )

        self.stats_df = stats_df
        self.stats_type = "z-prime"
        self.stats_method = method
        return stats_df

    def get_chem_structs_df(self, verbose: bool = True) -> pd.DataFrame:
        """
        Reads the chemical structure data from `self.chem_structs_path`
        and outputs a dataframe with the parent smiles (both with and
        without isomeric information).

        Args:
            verbose: Show error messages. Defaults to True.

        Returns:
            dataframe containing the columns: [Compound, Antineoplastic,
            iso_ParentSmiles, non_iso_ParentSmiles, Screening]

        """
        selchem_df = pd.read_csv(self.chemstructs_path, sep="\t")

        parent_smiles = list()
        idx_problem_mols = list()

        # Get parent smiles of all molecules, if error, save index
        smi_arr = selchem_df["Smiles"].values
        for idx, smi in enumerate(smi_arr):
            parent_smi, error = chembl_smi_standardizer(smi)
            if error:
                idx_problem_mols.append(idx)
            parent_smiles.append(parent_smi)
        if verbose:
            print(
                f"Parent smiles of {len(idx_problem_mols)} out of {len(smi_arr)}"
                " molecules could not be obtained"
            )
        # Parent smiles without isomeric information
        mol_arr = [Chem.MolFromSmiles(smi) for smi in parent_smiles]
        noniso_parents = [Chem.MolToSmiles(mol, isomericSmiles=False) for mol in mol_arr]

        # Adding antineoplastic activity flag
        key_words = re.compile(r"(antimitotic|antineoplastic)")
        antineo_idx = list()
        for idx, c in enumerate(selchem_df["Brief Description"]):
            if key_words.findall(c):
                antineo_idx.append(idx)
        flags = [True if idx in antineo_idx else False for idx in selchem_df.index]

        # Renaming and dropping some of the columns
        selchem_df = selchem_df.rename(
            columns={"Item Name": "Compound"},
        ).drop(
            columns=[
                "Catalog Number",
                "Concentration",
                "Plate Location",
                "Target",
                "Brief Description",
                "CAS Number",
            ],
        )

        # Updating the Compound names
        selchem_df["Compound"] = selchem_df["Compound"].apply(molname_clean)

        selchem_df = selchem_df.assign(
            Antineoplastic=flags,
            iso_ParentSmiles=parent_smiles,
            non_iso_ParentSmiles=noniso_parents,
            Screening=[self.name.capitalize()] * len(selchem_df),
        ).rename(columns={"Smiles": "SMILES"})
        return selchem_df

    def structure_mapping(
        self,
        structs_df: pd.DataFrame,
        smiles_type: str = "isomeric",
        use_pcp: bool = True,
        map_inplace: bool = True,
        verbose: bool = True,
        n_jobs: int = 5,
    ) -> dict:
        """
        Function to write a Compound name -> [SMILES, Antineo_flag] mapping
        as a dictionary and save that as a json file. Most of the SMILES are
        given by `structs_df`, but others are retrieved by querying PubChem.

        Note: This function should be called after `update_id_attributes()`

        Args:
            structs_df: output from `get_chem_structs_df`
            smiles_type: either `isomeric`, `non-isomeric` or `as-is`.
                Beware the server accepts max 5 requests/min. Defaults to 4.
            use_pcp: if True, also applies chemical structures retrieved from
                PubChem. Defaults to True.
            map_inplace: if True, maps to `self.df` and `self.norm_df`.
                Defaults to True.
            n_jobs: Retrieve smiles from pubchempy in parallel.

        Raises:
            ValueError: if n_jobs is greater than 5.
            ValueError: if smiles_type is not one of the three options.

        Returns:
            dict: with keys ['Structure_map', 'Antineo_map'], both holding
            names as keys and the respective information [SMILES, bool_flag].

        Note: The `Structure_map` is written using `structs_df`, which contains
        the columns `SMILES`, queried by the `smiles_type==as-is`, `iso_ParentSmiles`
        queried by `isomeric` and `non_iso_ParentSmiles`, queried by `non-isomeric`.
        """

        def map_identifiers(df, smiles_type: str, struct_dict: dict, use_pcp: bool = True):
            """
            Function to map the SMILES and the antineoplastic flagging
            based on compound names. Antineoplastic flagging will be np.Nan
            if the compound was retrieved from PubChem.

            Args:
                df: e.g. self.norm_df or self.df
                smiles_type: type of smiles structures to map to the dataframe.
                struct_dict: dictionary from `structure_mapping()`
                use_pcp: if true, will use pubchempy to retrieve smiles that did
                not have structures initially. Defaults to True.

            Returns:
                updated dataframe with mapped identifiers
            """
            allcomp_dict = {}
            assert smiles_type in ["isomeric", "non-isomeric", "as-is"], (
                "smiles_type must be one of the following: " "'isomeric', 'non-isomeric', 'as-is'"
            )
            sdf_struct_key = smiles_type + "_map"
            if smiles_type in ["isomeric", "as-is"]:
                pcp_smiles_type = "isomeric"
            else:
                pcp_smiles_type = "non-isomeric"

            allcomp_dict.update(struct_dict[sdf_struct_key])
            if use_pcp:
                allcomp_dict.update(struct_dict["pubchem_map"][pcp_smiles_type])

            df = df.assign(
                SMILES=df["Compound"].apply(lambda x: allcomp_dict[x] if x in allcomp_dict else np.NaN),
                Antineoplastic=df["Compound"].apply(
                    lambda x: allcomp_dict[x] if x in allcomp_dict else False
                ),
            )
            return df

        def rm_iso_info(smi: str):
            if isinstance(smi, str):
                try:
                    mol = Chem.MolFromSmiles(smi)
                    new_smi = Chem.MolToSmiles(mol, isomericSmiles=False)
                    return new_smi
                except:
                    return np.NaN
            else:
                return np.NaN

        # Script won't retrive SMILES again if .json is already saved
        json_path = self.file_root / f"{self.name}_structs.json"
        if json_path.exists():
            with json_path.open("r") as json_file:
                comp_info_dict = json.load(json_file)
            if verbose:
                print(f"Loaded structures saved under {self.file_root}")
            # Writing message with loaded structures
            pcp_structures = comp_info_dict["pubchem_map"]["isomeric"]
            sdf_structures = comp_info_dict["isomeric_map"]
            n_nan_pcp_structs = sum([1 for smi in pcp_structures.values() if isinstance(smi, float)])
            if verbose:
                print(
                    f"No structures found on the .SDF for {len(pcp_structures)} "
                    f"out of {len(pcp_structures) + len(sdf_structures)} compounds"
                )
                print(
                    f"From the {len(pcp_structures)} retrieved from pubchem, "
                    f"{n_nan_pcp_structs} were not found."
                )
            if map_inplace:
                if self.norm_df is not None:
                    self.norm_df = map_identifiers(self.norm_df, smiles_type, comp_info_dict)
                self.df = map_identifiers(self.df, smiles_type, comp_info_dict)
            self.comp_mapping = comp_info_dict
            return comp_info_dict

        if n_jobs > 5:
            raise ValueError("n_jobs must be less than 5 to avoid blocking by the server")

        all_names = structs_df["Compound"].values
        comp_with_struct = self.df[self.df["Compound"].isin(all_names)]["Compound"].unique()
        comp_without_struct = self.df[~self.df["Compound"].isin(all_names)]["Compound"].unique()

        subset_structs = structs_df[structs_df["Compound"].isin(comp_with_struct)].copy()

        print(subset_structs.columns)

        comp_info_dict = {
            "non-isomeric_map": {},
            "isomeric_map": {},
            "as-is_map": {},
            "Antineo_map": {},
        }
        comp_info_dict["non-isomeric_map"] = dict(
            zip(structs_df["Compound"], structs_df["non_iso_ParentSmiles"])
        )
        comp_info_dict["isomeric_map"] = dict(zip(structs_df["Compound"], structs_df["iso_ParentSmiles"]))
        comp_info_dict["as-is_map"] = dict(zip(structs_df["Compound"], structs_df["SMILES"]))
        comp_info_dict["Antineo_map"] = dict(zip(structs_df["Compound"], structs_df["Antineoplastic"]))
        # Getting structures of the compounds not in the dataframe in parallel
        if verbose:
            print("Retrieving structures not given by the .SDF from PubChem")
        with Pool(n_jobs) as pool:
            pcp_smiles = list(
                tqdm(
                    pool.imap(
                        # Sleep time = 1 to avoid blocking by the server
                        partial(
                            get_smiles_from_name,
                            smi_type="parent",
                            sleep_time=1,
                        ),
                        comp_without_struct,
                    ),
                    total=len(comp_without_struct),
                )
            )
        # Dictionary with mapping for smiles from pubchempy
        pcp_struct_dict = {
            name: (smi if isinstance(smi, str) else np.NaN)
            for name, smi in zip(comp_without_struct, pcp_smiles)
        }
        pcp_noniso_dict = {name: rm_iso_info(smi) for name, smi in zip(comp_without_struct, pcp_smiles)}
        comp_info_dict.update(
            {
                "pubchem_map": {
                    "isomeric": pcp_struct_dict,
                    "non-isomeric": pcp_noniso_dict,
                }
            }
        )
        # Giving the use information about the retrieved structures
        pcp_structures = comp_info_dict["pubchem_map"]["isomeric"]
        sdf_structures = comp_info_dict["isomeric_map"]
        n_nan_pcp_structs = sum([1 for smi in pcp_structures.values() if isinstance(smi, float)])
        if verbose:
            print(
                f"No structures found on the .SDF for {len(pcp_structures)} "
                f"out of {len(pcp_structures) + len(sdf_structures)} compounds"
            )
            print(
                f"From the {len(pcp_structures)} retrieved from pubchem, "
                f"{n_nan_pcp_structs} were not found."
            )
        # Write the mapping as a json file
        with json_path.open("w") as json_file:
            json.dump(comp_info_dict, json_file, indent=2)
        if map_inplace:
            if self.norm_df is not None:
                self.norm_df = map_identifiers(self.norm_df, smiles_type, comp_info_dict, use_pcp=use_pcp)
            self.df = map_identifiers(self.df, smiles_type, comp_info_dict, use_pcp=use_pcp)
        self.comp_mapping = comp_info_dict
        return comp_info_dict

    def npi_hitflagging(self, method="median", distance: float = 1.5):
        """
        Method to perform hitflagging on the npi-normalized dataframe. Hits will be
        flagged as cyst swelling-reducing when their `median + distance * MAD` is below
        the `median - distance * MAD` of the stimulant control (DMSO+FSK). The opposite
        opperation is used to flag cyst swelling-enhancing compounds, but with
        `distance * 2`.

        Args:
            method: Which method to apply the standard deviation. Defaults to "median".
            distance: Distance threshold from  the mean or the median to flag hit
                compounds. Defaults to 1.5.

        Raises:
            AttributeError: If the npi normalization was not calculated before

        Returns:
            pd.Dataframe containing hit_flagging columns, which is also assigned
            to self.norm_df
        """
        assert self.normalization_type == "npi", "Should apply npi normalization first."

        if self.norm_df is None:
            raise AttributeError("self.norm_df is not defined." "Hit flagging requires normalization")
        df = (
            self.norm_df.query('QC == "OK"')
            .reset_index(drop=True)
            .copy()
            .assign(Identifier=identifier_func())
        )
        controls_dict = dict()
        _idcols = ["Compound", "PlateID"]
        controls = (
            df.query('TreatmentType in ["pos_control", "neg_control"]')
            .loc[:, _idcols + ["obj.Mean(area)"]]
            .groupby(_idcols)
        )
        if method == "median":
            m_dict = controls.median()["obj.Mean(area)"].to_dict()
            dev_dict = controls["obj.Mean(area)"].apply(median_abs_deviation).to_dict()
        elif method == "mean":
            m_dict = controls.mean()["obj.Mean(area)"].to_dict()
            dev_dict = controls["obj.Mean(area)"].apply(np.mean).to_dict()
        controls_dict.update({"m": m_dict, "dev": dev_dict})

        # Hitflagging based solely on the cyst size
        median_vals = (
            df.loc[:, ["Identifier", "obj.Mean(area)"]]
            .groupby("Identifier")
            .median()
            .to_dict()["obj.Mean(area)"]
        )
        mad_vals = (
            df.loc[:, ["Identifier", "obj.Mean(area)"]]
            .groupby("Identifier")["obj.Mean(area)"]
            .apply(median_abs_deviation)
            .to_dict()
        )

        df = df.assign(
            Median_cystsize=df["Identifier"].map(median_vals),
            MAD_cystsize=df["Identifier"].map(mad_vals),
        )

        flagged_dfs = []
        for plate, sub_df in df.groupby("PlateID"):
            reducer_thresh = (
                controls_dict["m"][("DMSO+FSK", plate)] - distance * controls_dict["dev"][("DMSO+FSK", plate)]
            )
            enhancer_thresh = (
                controls_dict["m"][("DMSO+FSK", plate)]
                + distance * 2 * controls_dict["dev"][("DMSO+FSK", plate)]
            )
            # Bioactivity of the compounds
            reducer_boolean = (sub_df["Median_cystsize"] + distance * sub_df["MAD_cystsize"]) < reducer_thresh
            enhancer_boolean = (
                sub_df["Median_cystsize"] - distance * 2 * sub_df["MAD_cystsize"]
            ) > enhancer_thresh
            conditions = [np.array(reducer_boolean), np.array(enhancer_boolean)]
            choices = ["reducer", "enhancer"]
            sub_df = sub_df.assign(NPIScore_hitflag=np.select(conditions, choices, default="inactive"))
            flagged_dfs.append(sub_df)
        final_df = pd.concat(flagged_dfs, ignore_index=True).reset_index(drop=True)
        self.norm_df = final_df
        return final_df

    def zscore_hitflagging(self, method="median", distance: float = 1.5):
        """
        Method to perform hitflagging on the npi-normalized dataframe. Hits will be
        flagged as cyst swelling-reducing when their `median + distance * MAD` is below
        the `median - distance * MAD` of the stimulant control (DMSO+FSK). The opposite
        opperation is used to flag cyst swelling-enhancing compounds, but with
        `distance * 2`.

        Args:
            method: Which method to apply the standard deviation. Defaults to "median".
            distance: Distance threshold from  the mean or the median to flag hit
                compounds. Defaults to 1.5.

        Raises:
            AttributeError: If the z-score normalization was not calculated before

        Returns:
            pd.Dataframe containing hit_flagging columns, which is also assigned to self.norm_df
        """
        assert self.normalization_type == "z-score", "Should apply z-score normalization first."

        if self.norm_df is None:
            raise AttributeError("self.norm_df is not defined." "Hit flagging requires normalization")
        df = (
            self.norm_df.query('QC == "OK"')
            .reset_index(drop=True)
            .copy()
            .assign(Identifier=identifier_func())
        )
        controls_dict = dict()
        _idcols = ["Compound", "PlateID"]
        controls = (
            df.query('TreatmentType in ["pos_control", "neg_control"]')
            .loc[:, _idcols + ["obj.Mean(area)"]]
            .groupby(_idcols)
        )
        if method == "median":
            m_dict = controls.median()["obj.Mean(area)"].to_dict()
            dev_dict = controls["obj.Mean(area)"].apply(median_abs_deviation).to_dict()
        elif method == "mean":
            m_dict = controls.mean()["obj.Mean(area)"].to_dict()
            dev_dict = controls["obj.Mean(area)"].apply(np.mean).to_dict()
        controls_dict.update({"m": m_dict, "dev": dev_dict})

        # Hitflagging based solely on the cyst size
        median_vals = (
            df.loc[:, ["Identifier", "obj.Mean(area)"]]
            .groupby("Identifier")
            .median()
            .to_dict()["obj.Mean(area)"]
        )
        mad_vals = (
            df.loc[:, ["Identifier", "obj.Mean(area)"]]
            .groupby("Identifier")["obj.Mean(area)"]
            .apply(median_abs_deviation)
            .to_dict()
        )

        df = df.assign(
            Median_cystsize=df["Identifier"].map(median_vals),
            MAD_cystsize=df["Identifier"].map(mad_vals),
        )

        flagged_dfs = []
        for plate, sub_df in df.groupby("PlateID"):
            reducer_thresh = (
                controls_dict["m"][("DMSO+FSK", plate)] - distance * controls_dict["dev"][("DMSO+FSK", plate)]
            )
            enhancer_thresh = (
                controls_dict["m"][("DMSO+FSK", plate)]
                + distance * 2 * controls_dict["dev"][("DMSO+FSK", plate)]
            )
            # Bioactivity of the compounds
            reducer_boolean = (sub_df["Median_cystsize"] + distance * sub_df["MAD_cystsize"]) < reducer_thresh
            enhancer_boolean = (
                sub_df["Median_cystsize"] - distance * 2 * sub_df["MAD_cystsize"]
            ) > enhancer_thresh
            conditions = [np.array(reducer_boolean), np.array(enhancer_boolean)]
            choices = ["reducer", "enhancer"]
            sub_df = sub_df.assign(ZScore_hitflag=np.select(conditions, choices, default="inactive"))
            flagged_dfs.append(sub_df)
        final_df = pd.concat(flagged_dfs, ignore_index=True).reset_index(drop=True)
        self.norm_df = final_df
        return final_df

    def bscore_hitflagging(self, hit_rate: float = 0.01, perplate=False):
        """
        Function to perform the hit flagging on the results from `bscore_normalize`.
        Note that differently from NPI and Z-Score hitpicking, this function will take
        the x% top or bottom scoring compounds based on the `hit_rate` parameter. Since
        we observe lower amount of cyst-enhancing compounds in this assay, the hit_rate
        for the identification of enhancers is set to 10% of the actual hit_rate.

        Args:
            method: Which method to apply the standard deviation. Defaults to "median".
            hit_rate: Top or bottom % of compounds to be flagged as hits.
                Defaults to 0.01.

        Raises:
            AttributeError: If the b-score normalization was not calculated before

        Returns:
            pd.Dataframe containing hit_flagging columns, which is also assigned to self.norm_df
        """

        def b_hitpicking(df, hit_type: str, n_compounds: int, hit_rate: float):
            """
            Function to organize the dataframe and extract the top or the
            lowest % of compounds based on the hit rate.

            Returns the `Identifier` of the compounds (list)
            """
            # Dataframe will have top scoring compounds as first indexes
            if hit_type == "enhancer":
                ascending = False
            elif hit_type == "reducer":
                ascending = True
            identifiers = (
                df.sort_values("Median_cystsize", ascending=ascending)
                .reset_index(drop=True)
                .loc[: round(n_compounds * hit_rate)]["Identifier"]
                .unique()
                .tolist()
            )
            return identifiers

        assert self.normalization_type == "b-score", "Should apply b-score normalization first."

        if self.norm_df is None:
            raise AttributeError("self.norm_df is not defined." "Hit flagging requires normalization")
        df = (
            self.norm_df.query('QC == "OK"')
            .reset_index(drop=True)
            .copy()
            .assign(Identifier=identifier_func())
        )
        # Aggregate values observed for replicates using the median
        median_vals = (
            df.loc[:, ["Identifier", "obj.Mean(area)"]]
            .groupby("Identifier")
            .median()
            .to_dict()["obj.Mean(area)"]
        )
        df = df.assign(Median_cystsize=df["Identifier"].map(median_vals))

        comps_to_flag = {
            "reducer": [],
            "enhancer": [],
        }
        if perplate:
            for plate, sub_df in df.groupby("PlateID"):
                n_compounds = sub_df.query('TreatmentType == "treatment"')["Compound"].nunique()
                treats = sub_df.query('TreatmentType == "treatment"').drop_duplicates("Identifier")
                reducing_treats = b_hitpicking(treats, "reducer", n_compounds, hit_rate)
                enhancing_treats = b_hitpicking(treats, "enhancer", n_compounds, hit_rate / 10)
                comps_to_flag["reducer"].extend(reducing_treats)
                comps_to_flag["enhancer"].extend(enhancing_treats)
        else:
            n_compounds = df.query('TreatmentType == "treatment"')["Compound"].nunique()
            treats = df.query('TreatmentType == "treatment"').drop_duplicates("Identifier")
            reducing_treats = b_hitpicking(treats, "reducer", n_compounds, hit_rate)
            enhancing_treats = b_hitpicking(treats, "enhancer", n_compounds, hit_rate)
            comps_to_flag["reducer"].extend(reducing_treats)
            comps_to_flag["enhancer"].extend(enhancing_treats)

        conditions = [
            df["Identifier"].isin(comps_to_flag["reducer"]),
            df["Identifier"].isin(comps_to_flag["enhancer"]),
        ]
        choices = ["reducer", "enhancer"]
        df = df.assign(BScore_hitflag=np.select(conditions, choices, default="inactive"))
        self.norm_df = df
        return df

    def top_percent_hitpicking(self, hit_rate: float = 0.01, perplate=False):
        """
        Function to perform the same hit flagging as the `bscore_hitflagging` function,
        but implemented to work with the top % of compounds based on the `obj.Mean(area)`
        regardless of the normalization method.

        This function will take the x% top or bottom scoring compounds based on the
        `hit_rate` parameter. Since we observe lower amount of cyst-enhancing compounds
        in this assay, the hit_rate for the identification of enhancers is set to 10% of
        the actual hit_rate.

        Args:
            method: Which method to apply the standard deviation. Defaults to "median".
            hit_rate: Top or bottom % of compounds to be flagged as hits.
                Defaults to 0.01.

        Raises:
            AttributeError: If the b-score normalization was not calculated before

        Returns:
            pd.Dataframe containing hit_flagging columns, which is also assigned to self.norm_df
        """

        def top_percent_picking(df, hit_type: str, n_compounds: int, hit_rate: float):
            """
            Function to organize the dataframe and extract the top or the
            lowest % of compounds based on the hit rate.

            Returns the `Identifier` of the compounds (list)
            """
            # Dataframe will have top scoring compounds as first indexes
            if hit_type == "enhancer":
                ascending = False
            elif hit_type == "reducer":
                ascending = True
            identifiers = (
                df.sort_values("Median_cystsize", ascending=ascending)
                .reset_index(drop=True)
                .loc[: round(n_compounds * hit_rate)]["Identifier"]
                .unique()
                .tolist()
            )
            return identifiers

        assert self.normalization_type in [
            "z-score",
            "z-score_booij",
            "npi",
        ], "Make sure that the dataset is already normalized."

        if self.norm_df is None:
            raise AttributeError("self.norm_df is not defined." "Hit flagging requires normalization")
        df = (
            self.norm_df.query('QC == "OK"')
            .reset_index(drop=True)
            .copy()
            .assign(Identifier=identifier_func())
        )
        # Aggregate values observed for replicates using the median
        median_vals = (
            df.loc[:, ["Identifier", "obj.Mean(area)"]]
            .groupby("Identifier")
            .median()
            .to_dict()["obj.Mean(area)"]
        )
        df = df.assign(Median_cystsize=df["Identifier"].map(median_vals))

        comps_to_flag = {
            "reducer": [],
            "enhancer": [],
        }
        if perplate:
            for plate, sub_df in df.groupby("PlateID"):
                n_compounds = sub_df.query('TreatmentType == "treatment"')["Compound"].nunique()
                treats = sub_df.query('TreatmentType == "treatment"').drop_duplicates("Identifier")
                reducing_treats = top_percent_picking(treats, "reducer", n_compounds, hit_rate)
                enhancing_treats = top_percent_picking(treats, "enhancer", n_compounds, hit_rate / 10)
                comps_to_flag["reducer"].extend(reducing_treats)
                comps_to_flag["enhancer"].extend(enhancing_treats)
        else:
            n_compounds = df.query('TreatmentType == "treatment"')["Compound"].nunique()
            treats = df.query('TreatmentType == "treatment"').drop_duplicates("Identifier")
            reducing_treats = top_percent_picking(treats, "reducer", n_compounds, hit_rate)
            enhancing_treats = top_percent_picking(treats, "enhancer", n_compounds, hit_rate)
            comps_to_flag["reducer"].extend(reducing_treats)
            comps_to_flag["enhancer"].extend(enhancing_treats)

        conditions = [
            df["Identifier"].isin(comps_to_flag["reducer"]),
            df["Identifier"].isin(comps_to_flag["enhancer"]),
        ]
        choices = ["reducer", "enhancer"]
        norm_name = self.normalization_type.replace("_", "")
        df = df.assign(Top_percent=np.select(conditions, choices, default="inactive"))
        self.norm_df = df
        return df

    def get_antineoplastic_from_chembl(self):
        """
        Returns a dictionary with the names of antineoplastic compounds with
        the keys [iso_smiles] and [noniso_smiles].
        """

        def get_noniso_smiles(smi):
            mol = Chem.MolFromSmiles(smi)
            noniso_smi = Chem.MolToSmiles(mol, kekuleSmiles=False, canonical=True, isomericSmiles=False)
            return noniso_smi

        json_path = self.file_root / "chembl_antineostructs.json"
        if json_path.exists():
            with json_path.open("r") as json_file:
                antineo_drugs = json.load(json_file)
            print(f"Loaded antineoplastic structures saved under {self.file_root}")

            # Writing message with loaded structures
            print("Total number of antineoplastic molecules: ", len(antineo_drugs))
            return antineo_drugs

        molecule = new_client.molecule
        approved_drugs = molecule.filter(max_phase=4).order_by("molecule_properties__mw_freebase")
        approved_drugs.set_format("json")

        # L01 referes to the atc antineoplastic class, such as described in:
        # https://www.whocc.no/atc_ddd_index/?code=L01&showdescription=yes
        antineoplastic = re.compile("(L01)")
        drug_dict = dict()

        for drug in approved_drugs:
            name = drug["pref_name"]
            try:  # Get only approved drugs that have known SMILES
                drug["molecule_structures"]["canonical_smiles"]
            except TypeError:
                continue
            drug_dict[name] = {
                "chembl_id": drug["molecule_chembl_id"],
                "usan_definition": drug["usan_stem_definition"],
                "smiles": drug["molecule_structures"]["canonical_smiles"],
                "antineoplastic": False,
            }

            try:
                drug_dict[name].update(parent_chembl_id=drug["parent_chembl_id"])
            except KeyError:
                drug_dict[name].update(parent_chembl_id=None)

            try:
                drug_dict[name].update(atc_class=drug["atc_classifications"])
            except KeyError:
                drug_dict[name].update(atc_class=None)

            if drug_dict[name]["atc_class"] != []:
                # Some of the compounds are assigned multiple atc classes
                atc_class = " ".join(drug_dict[name]["atc_class"])
                if antineoplastic.findall(atc_class):
                    drug_dict[name]["antineoplastic"] = True

        antineo_drugs = dict()

        for key in drug_dict:
            if drug_dict[key]["antineoplastic"]:
                antineo_drugs[key] = dict()
                antineo_drugs[key]["chembl_id"] = drug_dict[key]["chembl_id"]
                antineo_drugs[key]["iso_smiles"] = drug_dict[key]["smiles"]
                antineo_drugs[key]["noniso_smiles"] = get_noniso_smiles(drug_dict[key]["smiles"])
        if not json_path.exists():
            with json_path.open("w") as json_file:
                json.dump(antineo_drugs, json_file, indent=1)
            print(f"Saved antineoplastic structures saved under {self.file_root}")
            print(  # Writing message with saved structures
                "Total number of antineoplastic molecules: ",
                len(antineo_drugs),
            )
        return antineo_drugs

    def flag_antineoplastic_compounds(
        self,
        which: str,
        antineo_drugs: dict,
        simi_threshold: float,
        fingerprint="ecfp4",
        n_jobs=1,
    ):
        """Adds a columns "Antineoplastic" to the dataframe with a boolean
        flag indicating if the compound is within the defined similarity threshold
        with approved antineoplastic (L01) drugs.

        Args:
            which: Which dataframe to flag. Options are "df", "norm_df" or "both".
            antineo_drugs: Dictionary with the antineoplastic compounds. Input from
            `get_antineoplastic_from_chembl`.
            simi_threshold: Similarity threshold to flag compounds as antineoplastic
            fingerprint: Type of fingerprint to use. Defaults to "ecfp4".
            n_jobs: Number of jobs for the similarity comparison. Defaults to 1.

        Raises:
            ValueError: If the dataframe is not available.

        Returns:
            None
        """
        if which not in ["df", "norm_df", "both"]:
            raise ValueError("which should be either 'df', 'norm_df' or 'both'")

        df = self.df.copy().dropna(subset=["SMILES"])
        known_flags = {  # Only take into consideration the true flags
            comp: flag
            for comp, flag in self.comp_mapping["Antineo_map"].items()
            if self.comp_mapping["Antineo_map"][comp]
        }
        antineo_smiles = [  # Antineoplastic compounds from ChEMBL
            antineo_drugs[key]["noniso_smiles"] for key in antineo_drugs.keys()
        ]
        data_smiles = df["SMILES"].unique().tolist()

        with Pool(n_jobs) as pool:
            antineo_fp_arr = pool.map(
                partial(smi_to_fp, fp_name=fingerprint), antineo_smiles
            )  # -> Fingerprints from antineoplastic compounds
            data_fp_arr = pool.map(
                partial(smi_to_fp, fp_name=fingerprint), data_smiles
            )  # -> Fingerprints from my own dataset

        idx_to_flag = list()
        print("Flagging antineoplastic compounds...")
        for antineo_fp in tqdm(antineo_fp_arr):
            for idx, fp in enumerate(data_fp_arr):
                # When we find a True bioactivity flag, we don't need to assess again
                # if known_flags[idx]:
                #     idx_to_flag.append(idx)
                #     continue
                if idx in idx_to_flag:
                    continue  # Similar to one compound? Don't assess again
                else:
                    tani = DataStructs.FingerprintSimilarity(antineo_fp, fp)
                    if tani >= simi_threshold:
                        idx_to_flag.append(idx)

        smi_toflag = np.take(data_smiles, idx_to_flag).tolist()
        print("Antineoplastic compounds within dataset: ", len(smi_toflag))
        smi_flag_mapping = {smi: True for smi in smi_toflag}

        def assign_antineo_flag(df: pd.DataFrame, smitoflag_map):
            """
            Args:
                df: self.df
                smitoflag_map: dictionary with smiles: flag
            Returns:
                flagged df
            """
            df = df.assign(Antineoplastic=df["SMILES"].map(smitoflag_map))
            df["Antineoplastic"] = df["Antineoplastic"].fillna(False)
            return df

        if which == "df":
            self.df = assign_antineo_flag(self.df, smi_flag_mapping)
        elif all([which in ["norm_df", "both"], self.norm_df is None]):
            self.df = assign_antineo_flag(self.df, smi_flag_mapping)
            print("norm_df is None. Function won't run on self.norm_df")
        elif which == "norm_df":
            self.norm_df = assign_antineo_flag(self.norm_df, smi_flag_mapping)
        elif which == "both":
            self.df = assign_antineo_flag(self.df, smi_flag_mapping)
            self.norm_df = assign_antineo_flag(self.norm_df, smi_flag_mapping)
        return

    def drop_smi_nan(self, which: str) -> None:
        """
        Drops rows without smiles

        Args:
            which: in which dataframe to drop the rows.
        """
        if which not in ["df", "norm_df", "both"]:
            raise ValueError("which should be either 'df', 'norm_df' or 'both'")

        starting_comps = len(self.df["Compound"].unique())

        if which == "df":
            self.df.drop(subset=["SMILES"], inplace=True).reset_index(drop=True)
        elif all([which in ["norm_df", "both"], self.norm_df is None]):
            self.df.drop(subset=["SMILES"], inplace=True).reset_index(drop=True)
            print("norm_df is None. Function won't run on self.norm_df")
        elif which == "norm_df":
            self.norm_df.drop(subset=["SMILES"], inplace=True).reset_index(drop=True)
        elif which == "both":
            self.df.drop(subset=["SMILES"], inplace=True).reset_index(drop=True)
            self.norm_df.drop(subset=["SMILES"], inplace=True).reset_index(drop=True)

        final_comps = len(self.df["Compound"].unique())
        print(f"Dropped {starting_comps - final_comps} compounds without smiles")
        return

    def drop_bad_qc(self, which: str) -> None:
        """
        Drops rows with QC != 'OK'

        Args:
            which: in which dataframe to drop the rows.
        """
        if which not in ["df", "norm_df", "both"]:
            raise ValueError("which should be either 'df', 'norm_df' or 'both'")

        starting_comps = len(self.df["Compound"].unique())

        if which == "df":
            self.df = self.df.query('QC == "OK"').reset_index(drop=True)
        elif all([which in ["norm_df", "both"], self.norm_df is None]):
            self.df = self.df.query('QC == "OK"').reset_index(drop=True)
            print("norm_df is None. Function won't run on self.norm_df")
        elif which == "norm_df":
            self.norm_df = self.norm_df.query('QC == "OK"').reset_index(drop=True)
        elif which == "both":
            self.df = self.df.query('QC == "OK"').reset_index(drop=True)
            self.norm_df = self.norm_df.query('QC == "OK"').reset_index(drop=True)

        final_comps = len(self.df["Compound"].unique())
        print(f"Dropped {starting_comps - final_comps} compounds with bad QC")
        return

    def barplot_z_prime(self, value_var):
        """
        Function for plotting the z-prime statistics of the selected feature
        (value_var).

        Params:
        value_var -> value_var for which to plot the z-prime statistics.

        For average size of the cystic objects, use `value_var=obj.Mean(area)`.
        """
        if self.z_prime_statistics is None:
            raise AttributeError("No z-prime statistics available")
        fig, ax = plt.subplots(figsize=(6, 3))
        self.stats_df[self.stats_df["Measurement"] == value_var].plot.bar(x="PlateID", ax=ax)

        for tick in ax.get_xticklabels():
            tick.set_rotation(0)
            tick.set_ha("right")
        return fig, ax

    def lineplot_z_prime(self, feature_type: str, ylim: tuple = None):
        """
        Function for plotting the z-prime statistics of (almost) all
        the features available.

        Params:
        feature_type -> either `morphology` or `texture`
        """
        if self.stats_df is None:
            raise AttributeError("No z-prime statistics available")
        elif feature_type not in ["morphology", "texture"]:
            raise AttributeError("feature_type unavailable. Should be either `morphology` or `texture`")

        # Removing information for the standard deviations:
        df = self.stats_df
        df = df[~(df["Measurement"].str.contains("\.SD", case=True, regex=True))]

        # Subset df for feature type and removing pattern matching cols
        avoid_pattern = "Count|Sum|\.corr\.|maximum|minimum|connection"
        morpho_pattern = "_order_"

        if feature_type == "morphology":
            titles = [
                f"Z' values for {feature_type} features per plate\nNuclei mask",
                f"Z' values for {feature_type} features per plate\nOrganoid mask",
            ]
            nc_subset_df = df[
                ~(df["Measurement"].str.contains(avoid_pattern))
                & ~(df["Measurement"].str.contains(morpho_pattern))
                & (df["Measurement"].str.startswith("nc."))
            ]
            obj_subset_df = df[
                ~(df["Measurement"].str.contains(avoid_pattern))
                & ~(df["Measurement"].str.contains(morpho_pattern))
                & (df["Measurement"].str.startswith("obj."))
            ]

        elif feature_type == "texture":
            titles = [
                f"Z' values for {feature_type} features per plate\nNuclei mask",
                f"Z' values for {feature_type} features per plate\nOrganoid mask",
            ]
            nc_subset_df = df[
                ~(df["Measurement"].str.contains(avoid_pattern))
                & (df["Measurement"].str.contains(morpho_pattern))
                & (df["Measurement"].str.startswith("nc."))
            ]
            obj_subset_df = df[
                ~(df["Measurement"].str.contains(avoid_pattern))
                & (df["Measurement"].str.contains(morpho_pattern))
                & (df["Measurement"].str.startswith("obj."))
            ]

        fig, axs = plt.subplots(nrows=2, figsize=(5, 6))
        for ax, title, df in zip(axs, titles, [nc_subset_df, obj_subset_df]):
            lineplot = sns.lineplot(
                data=df,
                y="Z-prime",
                hue="Measurement",
                x="PlateID",
                palette="hsv",
                ax=ax,
            )
            ax.legend(
                bbox_to_anchor=(1.04, 1),
                borderaxespad=0,
                ncol=2,
                prop={"size": 7},
            )
            ax.xaxis.set_tick_params(labelsize=8)
            ax.yaxis.set_tick_params(labelsize=8)
            ax.set_ylabel("Z-prime", fontsize=8)
            ax.set_xlabel("Plate number", fontsize=8)
            ax.set_title(title, size=10)
            if ylim is not None:
                ax.set_ylim(ylim)

        fig.subplots_adjust(hspace=0.5)
        return fig, ax

    def calculate_spearman_corr(
        self,
        comp_subset: list = ["DMSO", "DMSO+FSK"],
        discrete: bool = True,
        njobs=5,
    ) -> Union[pd.DataFrame, pd.DataFrame]:
        """
        Function for calculating the spearman correlation within the available
        features. This function should be called after applying the z-score
        normalization on the data (takes `self.norm_df` as input).

        >>> # Example of usage:
        >>> corr_df, signif_df = selchem.calculate_spearman_corr(comp_subset=None)
        >>> # Results can be visualized with:
        >>> fig, ax = selchem.plot_spearman_corr_heatmap(corr_df, signif_df)

        Args:
            comp_subset: Subset of compounds to calculate feature correlation.
            If `comp_subset is None`, all compounds will be used.
            Defaults to ["DMSO", "DMSO+FSK"].
            discrete: Boolean for whether to return discrete significance
                [yes, no] or the actual p-value. Defaults to True.
            njobs: number of jobs to run spearman correlation in parallel.
                Defaults to 5.

        Raises:
            AttributeError: when self.norm_df is None.
                Should be done on normalized data..

        Returns:
            spearman_corr -> Spearman correlation
            [corr_df | pv_df] -> [1] as discrete signif. or continuous p-values.
        """
        if self.norm_df is None:
            raise AttributeError("No normalized data available. Use z-score normalization")

        # Removing information about standard deviations
        dev_pattern = re.compile("\.SD")
        final_features = [c for c in self.norm_df.columns if not dev_pattern.findall(c)]

        norm_df = self.norm_df[final_features]
        val_cols = [c for c in norm_df.columns if c not in ["Compound"] + self.id_cols]

        # Subset df for comp_subset & drop nan values
        if comp_subset is None:
            subset = norm_df.dropna()
        else:
            subset = norm_df[norm_df["Compound"].isin(comp_subset)].dropna()

        # Calculate spearman correlation for all non-redundant feature combinations
        combi = list(combinations(val_cols, 2))
        with Pool(njobs) as pool:
            results = pool.starmap(partial(get_spearman_corrs, subset_df=subset), combi)

        # Nested dictionaries to hold the values for the dataframe
        pvals_dict = {val: {v: np.NaN for v in val_cols} for val in val_cols}
        corrs_dict = {val: {v: np.NaN for v in val_cols} for val in val_cols}

        pvals_to_correct = list()
        spear_corrs = list()
        key_names = list()

        # unpacking results from multiprocessing
        # values are held by a key in the form of (feature1~feature2)
        for d in results:
            key = list(d.keys())[0]
            col1, col2 = key.split("~")
            corr, pval = [i for i in d[key]]
            pvals_dict[col1][col2] = pval
            corrs_dict[col1][col2] = corr
            pvals_to_correct.append(pval)
            spear_corrs.append(corr)
            key_names.append(key)

        # Correcting pvalues for multiple testing and updating original dictionary
        multitest_result = multipletests(pvals=pvals_to_correct, method="fdr_bh")
        self.multitest_results = multitest_result
        self.stats_test_type = "Spearman Correlation"
        print(
            "Number of null hypothesis that were rejected:",
            np.count_nonzero(multitest_result[0] == True),
        )
        print(
            "Number of null hypothesis that were accepted:",
            np.count_nonzero(multitest_result[0] == False),
        )

        # Updating the dictionary with the corrected pvalues
        corrected_pvals = multitest_result[1]
        for pv, k in zip(corrected_pvals, key_names):
            col1, col2 = k.split("~")
            pvals_dict[col1][col2] = pv

        pv_heatmap_df = pd.DataFrame.from_dict(pvals_dict, orient="index")
        corr_heatmap_df = pd.DataFrame.from_dict(corrs_dict, orient="index")
        signif_heatmap_df = pv_heatmap_df.apply(discrete_transform, axis=0)

        if discrete:
            return corr_heatmap_df, signif_heatmap_df
        else:
            return corr_heatmap_df, pv_heatmap_df

    @staticmethod
    def plot_spearman_corr_heatmap(
        corr_heatmap_df: pd.DataFrame,
        signif_heatmap_df: pd.DataFrame,
        fig_path: str = None,
    ) -> plt.Figure:
        """
        Function to plot the spearman correlation heatmap.
        Takes as input the two outputs from the `calculate_spearman_corr` method.
        Usage of this function is also documented there.

        Note:
        Plotting the continuous p-values was avoided as it makes it harder to visualize.

        Args:
            corr_heatmap_df: df with the spearman correlation.
            signif_heatmap_df: df with the significance (yes, no) of the spearman correlation.
            fig_path: path for saving the figure. Defaults to None. (won't save it)

        Returns:
            fig, ax -> figure and axis of the plot
        """
        toplot_signif = signif_heatmap_df.to_numpy()
        # toplot_pvals = pv_heatmap_df.to_numpy().T # if we want to plot the pvalues
        toplot_corr = corr_heatmap_df.to_numpy().T

        fig, ax = plt.subplots(figsize=(20, 15))

        """
        Adding custom axis for the color bar 
        Credits https://stackoverflow.com/questions/67035996/seaborn-heatmap-colorbar-custom-location
        
        First part of the plot:
        -> binary (significance) heatmap on the top right of the figure
        """
        cax = inset_axes(
            ax,
            width="1%",
            height="40%",
            loc="upper left",
            bbox_to_anchor=(1.01, 0, 1, 1),  # Custom position
            bbox_transform=ax.transAxes,
            borderpad=0,
        )

        heatmap1 = sns.heatmap(
            toplot_signif,
            ax=ax,
            xticklabels=False,
            yticklabels=False,
            vmin=0,
            vmax=1,
            cmap="bwr",
            cbar_ax=cax,
            cbar_kws={
                "label": "Significance\n(p-value < 0.05)",
                "fraction": 0.02,
                "pad": 0.10,
            },
        )

        """
        Customizing the tick locators to suit the binary label: significant
            [True, False] For further usage reference, check:
        https://matplotlib.org/3.4.3/gallery/ticks_and_spines/tick-locators.html
        """

        cax.yaxis.set_major_locator(ticker.FixedLocator([0, 1]))
        cax.set_yticklabels(["No", "Yes"])

        """
        Second part of the plot:
        -> Spearman correlation heatmap on the bottom left of the figure
         """
        cax2 = inset_axes(
            ax,
            width="1%",
            height="40%",
            loc="lower left",
            bbox_to_anchor=(1.01, 0, 1, 1),
            bbox_transform=ax.transAxes,
            borderpad=0,
        )
        heatmap2 = sns.heatmap(
            toplot_corr,
            ax=ax,
            xticklabels=corr_heatmap_df.columns,
            yticklabels=corr_heatmap_df.index,
            vmin=-1,
            vmax=1,
            cmap="bwr",
            cbar_ax=cax2,
            cbar_kws={
                "label": "Spearman Correlation",
                "fraction": 0.02,
                "pad": 0.01,
            },
        )

        for tick in ax.get_yticklabels():
            tick.set_size(6)
        for tick in ax.get_xticklabels():
            tick.set_rotation(45)
            tick.set_horizontalalignment("right")
            tick.set_size(6)

        ax.tick_params("both")
        ax.set_title(
            "Heatmap of the Spearman correlation within the different features (bottom)\n \
            and their respective statistical significance (top)",
            fontsize=14,
        )

        if fig_path is not None:
            heatmap2.get_figure().savefig(
                fig_path / "heatmap_significance.png",
                dpi=600,
                bbox_inches="tight",
            )

        return fig, ax

    def spearman_most_signif_features(self, signif_heatmap_df, feature_type: None):
        """
        Function to get the most significant features from the spearman
        correlation heatmap. Takes as input discrete output from the
        `calculate_spearman_corr` method.

        Params:
        signif_heatmap_df -> df with the significance (yes, no) of the spearman
            correlation.
        feature_type -> either `morphology` or `texture`.
        """
        if feature_type not in ["morphology", "texture"]:
            raise AttributeError("feature_type unavailable. Should be either " "`morphology` or `texture`")
        # Summing over two axis as the dataframe contains only non-redundant comparisons
        compressed_signif = signif_heatmap_df.sum(axis=0) + signif_heatmap_df.sum(axis=1)

        fig, ax = plt.subplots(figsize=(20, 3))

        if feature_type == "morphology":
            signif_df = compressed_signif[~compressed_signif.index.str.contains("_order_")].sort_values(
                ascending=True
            )

        elif feature_type == "texture":
            signif_df = compressed_signif[compressed_signif.index.str.contains("_order_")].sort_values(
                ascending=True
            )

        ax = signif_df.plot.bar()
        ax.set_title(
            f"Number of non-redundant significant spearman's correlations within "
            f"{len(signif_df)} {feature_type} features:\n{self.name} screening"
        )

        ax.set_ylabel("N signif. correlations")
        for tick in ax.get_xticklabels():
            tick.set_rotation(45)
            tick.set_horizontalalignment("right")
            tick.set_fontsize(8)
        ax.set_xlabel("Features")
        return fig, ax

    def plot_corrected_pvalues(self):
        """
        Function that outputs a lineplot with two `axs`, one containing
        the p-values that rejected the null-hypothesis and another that
        contains the p-values that were not rejected.
        """
        if self.multitest_results is None:
            raise AttributeError(
                "Please run a function that performs " "FDR correction first. Check documentation"
            )

        multitests = self.multitest_results

        fig, axs = plt.subplots(ncols=2, figsize=(9, 3))

        accepted = np.compress(multitests[0], multitests[1])
        rejected = np.compress(~multitests[0], multitests[1])

        for ax, idx in zip(axs, range(2)):
            if idx == 0:
                ax.plot(accepted, color="tab:cyan")
                ax.hlines(
                    y=0.05,
                    xmin=0,
                    xmax=len(accepted),
                    colors=["red"],
                    linestyles="dashed",
                )
                ax.set_title(
                    f"{len(accepted)} FDR-corrected p-values with rejected null-hypothesis:"
                    f"\n{self.stats_test_type} - {self.name.capitalize()} screening",
                    fontsize=10,
                )
                ax.set_ylim(0, 0.06)
                ax.set_xlim(0, len(accepted))
                ax.set_ylabel("p-value")
                ax.set_xlabel("Statistical tests (number)")
            else:
                ax.plot(rejected, color="tab:cyan")
                ax.hlines(
                    y=0.05,
                    xmin=0,
                    xmax=len(rejected),
                    colors=["red"],
                    linestyles="dashed",
                    label="p-value = 0.05",
                )
                ax.set_title(
                    f"{len(rejected)} FDR-corrected p-values with accepted null-hypothesis:"
                    f"\n{self.stats_test_type} - {self.name.capitalize()} screening",
                    fontsize=10,
                )
                ax.set_ylim(0, 1)
                ax.set_xlim(0, len(rejected))
                ax.set_ylabel("p-value")
                ax.set_xlabel("Statistical tests (number)")
                ax.legend(
                    bbox_to_anchor=(1, -0.07),
                    loc="lower right",
                    bbox_transform=fig.transFigure,
                    ncol=3,
                )

        # Removing grid from the x axis
        plt.grid(False, axis="x")
        plt.tight_layout()
        fig.subplots_adjust(wspace=0.5)
        return fig, axs

    def calculate_mann_whitney_u_signif(
        self, id_col: str = "Compound", identifiers: list = ["DMSO", "DMSO+FSK"]
    ):
        """
        Calculates the statistical significance for the data distributions
        of the two groups. Statistical method for testing this hypothesis is the
        Mann-Whitney U test.

        Example:
        >>> selchem.calculate_mann_whitney_u_signif(id_col='TreatmentType',
        >>>                 identifiers=['pos_control', 'treatment_control'])

        Args:
            id_col: column name of the identifier column. Defaults to "Compound".
            identifiers: list of identifiers to be tested. Defaults to ["DMSO", "DMSO+FSK"].

        Returns:
            mann_whit_df -> dataframe with the p-values and the statistical significance.
        """
        dev_pattern = re.compile("\.SD")
        final_features = [c for c in self.norm_df.columns if not dev_pattern.findall(c)]

        norm_df = self.norm_df[final_features].copy().dropna()

        # Selecting all the measurements that are not identifiers
        features_only = [c for c in norm_df.columns if c not in self.id_cols + ["Compound"]]

        # Creating a nested dictionary of all features per plate
        stats_results = {p: {c: np.NaN for c in features_only} for p in norm_df["PlateID"].unique()}
        pvals_to_correct = dict()

        perplate_grouped = norm_df.groupby(by="PlateID")
        for plate, df in perplate_grouped:
            solvent = df[df[id_col] == identifiers[0]]
            stimulant = df[df[id_col] == identifiers[1]]
            for featur in features_only:
                result = mannwhitneyu(solvent[featur], stimulant[featur])
                pvals_to_correct[f"{plate}~{featur}"] = result.pvalue

        pval_arr = list(pvals_to_correct.values())
        multitest_result = multipletests(pvals=pval_arr, method="fdr_bh")
        self.multitest_results = multitest_result
        self.stats_test_type = "Mann-Whitney U"
        print(
            "Number of hypothesis that were rejected:",
            np.count_nonzero(multitest_result[0] == True),
        )
        print(
            "Number of hypothesis that were accepted:",
            np.count_nonzero(multitest_result[0] == False),
        )
        corrected_pvals = multitest_result[1]

        for idx, k in enumerate(pvals_to_correct.keys()):
            plate, feature = k.split("~")
            plate = int(plate)
            stats_results[plate][feature] = corrected_pvals[idx]

        mann_whit_df = pd.DataFrame.from_dict(stats_results, orient="index")
        return mann_whit_df

    def plot_mann_whitney_u_heatmap(
        self,
        mann_whit_df,
        figsize=(18, 3),
        title: str = False,
        feature_type: str = None,
        signal_type: str = None,
        fig_path=None,
        discrete=False,
        cmap=None,
    ):
        """
        Function to plot the heatmap of the statistical significance of the difference
        within data distributions tested by the `test_mann_whitney_u_signif` method.

        Example of usage:
        >>> fig, ax = selchem.plot_mann_whitney_u_heatmap(mann_whit_df,
        >>>                            figsize=(16,3),
        >>>                            feature_type='texture',
        >>>                            signal_type='objects',
        >>>                            discrete=True,
        >>>                            cmap='hot'
        >>>                            )

        Args:
            mann_whit_df: output from the method `calculate_mann_whitney_u_signif`.
            figsize: size of the figure. Defaults to (18, 3).
            title: title of the heatmap. Defaults to False.
            feature_type: either `morphology` or `texture`. Defaults to None.
            signal_type: either `nuclei` or `objects`. Defaults to None.
            fig_path:  path to save the figure. Defaults to None. (won't save)
            discrete: whether heatmap should be discrete or not. Defaults to False.
            cmap: colormap to use for the heatmap. Defaults to None. -> 'Reds'

        Raises:
            AttributeError: feature_type not in ["morphology", "texture"]
            AttributeError: signal_type not in ["nuclei", "objects"]

        Returns:
            fig, ax -> matplotlib figure and axis
        """
        if feature_type not in ["morphology", "texture"]:
            raise AttributeError("feature_type unavailable. Should be either `morphology` or `texture`")
        if signal_type not in ["nuclei", "objects"]:
            raise AttributeError("signal_type unavailable. Should be either `nuclei` or `objects`")

        # Subsetting the data according to input parameters:
        if feature_type == "texture":
            featur_signif = mann_whit_df.loc[:, (mann_whit_df.columns.str.contains("_order_"))]
        elif feature_type == "morphology":
            featur_signif = mann_whit_df.loc[:, ~(mann_whit_df.columns.str.contains("_order_"))]

        if signal_type == "nuclei":
            toplot_df = featur_signif.loc[:, featur_signif.columns.str.startswith("nc.")]
        elif signal_type == "objects":
            toplot_df = featur_signif.loc[:, featur_signif.columns.str.startswith("obj.")]

        fig, ax = plt.subplots(figsize=figsize)

        if cmap is None:
            cmap = "Reds"

        if cmap.endswith("_r"):
            center = 0.15
        else:
            center = 0.85

        if discrete:
            toplot_df = toplot_df.apply(discrete_transform)
            heatmap = sns.heatmap(toplot_df, cmap=cmap, center=center, linewidths=0.1)
        else:
            heatmap = sns.heatmap(
                toplot_df,
                annot=True,
                fmt=".2f",
                cmap=cmap,
                center=center,
                annot_kws={"size": 6},
            )

        # heatmap = sns.heatmap(df, annot=True, fmt=".2f", cmap="Blues_r", annot_kws={'size': 8})
        if title:
            ax.set_title(title)

        # Set tick label sizes
        for tick in ax.get_yticklabels():
            tick.set_size(8)
            tick.set_ha("right")
        for tick in ax.get_xticklabels():
            tick.set_size(8)
            tick.set_rotation(45)
            tick.set_ha("right")

        ax.title.set_fontsize(10)
        ax.set_xlabel("Extracted Features", fontsize=8)
        ax.set_ylabel("Plate Number", fontsize=8)
        fig = heatmap.get_figure()

        if fig_path is not None:
            fig.savefig(Path(fig_path), format="png", bbox_inches="tight")

        return fig, ax

    def mann_whitney_u_signif_perplate(
        self,
        mann_whit_df,
        figsize=(12, 3),
        title: str = False,
        feature_type: str = None,
        signal_type: str = None,
        fig_path=None,
    ):
        """
        Function to plot the summary of the heatmap. Each feature signifiance is summed by
        plate numbers and the result is shown as a barplot.

        Params:
        mann_whit_df -> output from the method `calculate_mann_whitney_u_signif`.
        figsize -> size of the figure.
        title -> title of the barplot.
        feature_type -> either `morphology` or `texture`.
        signal_type -> either `nuclei` or `objects`.
        fig_path -> path to save the figure. [default: `None` -> not saved]
        """

        if feature_type not in ["morphology", "texture"]:
            raise AttributeError("feature_type unavailable. Should be either `morphology` or `texture`")
        if signal_type not in ["nuclei", "objects"]:
            raise AttributeError("signal_type unavailable. Should be either `nuclei` or `objects`")

        # Subsetting the data according to input parameters:
        if feature_type == "texture":
            featur_signif = mann_whit_df.loc[:, (mann_whit_df.columns.str.contains("_order_"))]
        elif feature_type == "morphology":
            featur_signif = mann_whit_df.loc[:, ~(mann_whit_df.columns.str.contains("_order_"))]

        if signal_type == "nuclei":
            toplot_df = featur_signif.loc[:, featur_signif.columns.str.startswith("nc.")]
        elif signal_type == "objects":
            toplot_df = featur_signif.loc[:, featur_signif.columns.str.startswith("obj.")]

        fig, ax = plt.subplots(figsize=figsize)
        ax = toplot_df.apply(discrete_transform).sum(axis=0).sort_values(ascending=True).plot.bar()
        # Grabbing distributions that were signficant for 25% of the plates:
        # mask = toplot_df.apply(discrete_transform).sum(axis=0) >= len(toplot_df/4)
        # self.most_significant = toplot_df.apply(discrete_transform).sum(axis=0)[mask].index
        if title:
            ax.set_title(title)
        ax.set_ylabel("N significant tests\np-value<0.05")
        for tick in ax.get_xticklabels():
            tick.set_rotation(45)
            tick.set_horizontalalignment("right")
            tick.set_fontsize(8)

        ax.set_xlabel("Features")
        return fig, ax


class spectrum(selleck_chem):
    """
    Class for processing the spectrum screening data
    """

    def __init__(self, name="spectrum") -> None:
        super().__init__(name=name)  # __init__ content from selleck_chem class
        self.file_path = list(self.file_root.glob(f"*{self.name}_Batch*.csv.gz"))[0]
        self.chemstructs_path = (
            self.root_dir / "data/adpkd_screening/chemical_structures/combined_SPECTRUM_structures.sdf"
        )
        self.df = pd.read_csv(self.file_path)
        # id_cols are not the same; there's only available compound concentration.
        self.control_treatments = [
            "bez-235",
            "rapamycin",
            "roscovitine",
            "sorafenib",
            "metformin",
        ]
        self.toxic_treatments = [
            "daunorubicin",
            "doxorubicin",
            "gambogic acid",
            "epirubicin hydrochloride",
            # Toxic according to the publication by Tijmen Booij et al.
        ]
        self.outlier_wells = [
            "N20_2",
            "A1_15",
            "A1_4",
            "I2_10",
            "L23_3",
            "B23_7",
            "J24_13",
            "B24_6",
        ]

    def drop_cols(self) -> None:
        """
        Drops undersired columns from self.df
        """
        # Need to keep the column plate.folder
        remove_pattern = re.compile(r"path|name|\.tif|\.csv|code|folder|^row$|root")
        to_drop_cols = [c for c in self.df.columns if all([remove_pattern.findall(c), c != "plate.folder"])]
        self.df.drop(columns=to_drop_cols, inplace=True)

        rename_pattern = re.compile(r"plate\.layout\.info\.")
        to_rename = {  # Removing this prefix as we don't need it
            key: value
            for key, value in zip(
                self.df.columns,
                [rename_pattern.sub("", c) for c in self.df.columns],
            )
        }
        self.df.rename(columns=to_rename, inplace=True)
        return

    def update_stimulation(self) -> None:
        """
        Overwriting method as the stimulation already follows the pattern we want
        """
        pass

    def update_id_attributes(self) -> None:
        """
        1) Updates the indentification column names from self.df for
        further normalization steps.
        2) Sorts the dataframe by the identification columns
        PlateID, PlateColumn, PlateRow
        """
        # Changing attributes that are different from selleck_chem
        self.df["plate.folder"] = self.df["plate.folder"].str.split("\\").str[-1]
        self.df.rename(columns={"MoleculeName": "Compound"}, inplace=True)
        # Concentration of control treatments not registered on the plate layouts
        # I assign 1 to all, as it's the concentration for the treatment compounds
        self.df["Concentration"] = [1 for _ in range(len(self.df))]

        # Standardizing compound names with molname_clean
        standard_names = {name: molname_clean(name) for name in self.df["Compound"].unique()}
        self.df.replace({"Compound": standard_names}, inplace=True)

        # mapping the types of treatment within the screening
        numerics = ["int16", "int32", "int64", "float16", "float32", "float64"]
        conditions = [
            (self.df["Compound"] == "DMSO"),
            (self.df["Compound"] == "DMSO+FSK"),
            (self.df["Compound"].str.lower().isin(self.control_treatments)),
            (
                (self.df["Compound"].str.lower().isin(self.toxic_treatments))
                # & (self.df["Concentration"] >= 1)
            ),
        ]
        choices = [
            "neg_control",
            "pos_control",
            "treatment_control",
            "treatment_toxic",
        ]
        self.df["TreatmentType"] = np.select(conditions, choices, default="treatment")

        newdf = (
            self.df.select_dtypes(include=numerics)
            .assign(
                Compound=self.df["Compound"],
                TreatmentType=self.df["TreatmentType"],
                PlateID=self.df["plate.folder"],
                PlateRow=self.df["row.char"],
                PlateColumn=self.df["column"],
                QC=self.df["QC"],
            )
            .drop(
                columns=["column", "Volume"],
                # `column` dropped as it's numerical and remains after select_dtypes
            )
        )

        # Fixing column order, changing PlateID to numerical and resetting index
        self.df = (
            newdf[["Compound"] + self.id_cols + newdf.columns.tolist()[1:-6]]
            .assign(PlateID=newdf["PlateID"].astype(int))
            .sort_values(by=["PlateID", "PlateRow", "PlateColumn"])
            .reset_index(drop=True)
        )
        return

    def get_chem_structs_df(self, verbose: bool = True) -> pd.DataFrame:
        """
        Reads the chemical structure data from `self.chem_structs_path`
        and outputs a dataframe with the parent smiles (both with and
        without isomeric information).

        Args:
            verbose: Show error messages. Defaults to True.

        Returns:
            dataframe containing the columns: [Compound, Antineoplastic,
            iso_ParentSmiles, non_iso_ParentSmiles, Screening]

        """
        # Reading the utf-8 file
        spectrum_df: pd.DataFrame = PandasTools.LoadSDF(str(self.chemstructs_path))

        # Checking for molecules that weren't parsed:
        comp_names = list()
        with open(self.chemstructs_path, "r") as source_file:
            for line in source_file:
                if line.startswith(">  <MOLENAME>"):
                    comp_names.append(next(source_file).rstrip("\n"))

        all_idxs = np.arange(spectrum_df.index[-1] + 1)
        failed_idxs = np.setdiff1d(all_idxs, spectrum_df.index.values)

        for idx in failed_idxs:
            self.dropped_comps[comp_names[idx]] = "Smiles not parsed by RDKit"

        # Adding antineoplastic activity flag
        key_words = re.compile(r"(antimitotic|antineoplastic)")

        spectrum_df = spectrum_df.reset_index(drop=True).fillna("", downcast=None)
        antineo_idxs = [idx for idx, c in enumerate(spectrum_df["therap"]) if key_words.findall(c)]

        flags = [True if idx in antineo_idxs else False for idx in spectrum_df.index]

        # Editing the dataframe:
        spectrum_df.drop(
            columns=["ID", "ref", "cas#", "source", "therap", "tradename"],
            inplace=True,
        )

        # Get parent smiles of all molecules, if error, save index
        mol_arr = np.array(spectrum_df["ROMol"])
        idx_problem_mols = list()
        parent_smiles = list()
        for idx, mol in enumerate(mol_arr):
            parent_smi, error = chembl_mol_standardizer(mol)
            if error:
                idx_problem_mols.append(idx)
            parent_smiles.append(parent_smi)
        if verbose:
            print(
                f"Parent smiles of {len(idx_problem_mols)} out of {len(mol_arr)}"
                " molecules could not be obtained"
            )
        noniso_parents = [Chem.MolToSmiles(mol, isomericSmiles=False) for mol in mol_arr]

        spectrum_df = (
            # updating dataframe and returning result
            spectrum_df.assign(
                MOLENAME=spectrum_df["MOLENAME"].apply(molname_clean),
                Antineoplastic=flags,
                iso_ParentSmiles=parent_smiles,
                non_iso_ParentSmiles=noniso_parents,
                Screening=[self.name.capitalize()] * len(spectrum_df),
            )
            .rename(columns={"MOLENAME": "Compound"})
            .drop(columns=["ROMol", "status"])
        )
        return spectrum_df


class spectrum_validation(spectrum):
    """
    Class for processing the spectrum screening data
    """

    def __init__(self, name="spectrum-validation") -> None:
        super().__init__(name=name)  # __init__ content from selleck_chem class
        self.chemstructs_path = (
            self.root_dir / "data/adpkd_screening/chemical_structures/combined_SPECTRUM_structures.sdf"
        )
        self.file_path = list(self.file_root.glob(f"*{self.name}_Batch*.csv.gz"))[0]
        self.df = pd.read_csv(self.file_path)
        # id_cols are not the same; there's only available compound concentration.
        self.control_treatments = [
            "bez-235",
            "rapamycin",
            "roscovitine",
            "sorafenib",
        ]
        self.toxic_treatments = [
            "staurosporin",
            "daunorubicin",
            "doxorubicin",
            "gambogic acid",
            "epirubicin hydrochloride",
            # Toxic according to the publication by Tijmen Booij et al.
        ]
        self.outlier_wells = ["J4_1"]

    def update_id_attributes(self) -> None:
        """
        1) Updates the indentification column names from self.df for
        further normalization steps.
        2) Sorts the dataframe by the identification columns
        PlateID, PlateColumn, PlateRow
        """
        self.df["Compound"].fillna("DMSO", inplace=True, downcast=None)  # BAD QC wells named as DMSO
        # Changing attributes that are different from selleck_chem
        self.df["plate.folder"] = self.df["plate.folder"].str.split("\\").str[-1]
        # Standardizing compound names with molname_clean
        standard_names = {name: molname_clean(name) for name in self.df["Compound"].unique()}
        self.df.replace({"Compound": standard_names}, inplace=True)
        self.df["Compound"] = self.df.apply(
            lambda x: x["Compound"] if x["Compound"] + x["Exposure"] != "DMSOForskolin" else "DMSO+FSK",
            axis=1,  # Updating positive control name
        )

        # mapping the types of treatment within the screening
        numerics = ["int16", "int32", "int64", "float16", "float32", "float64"]
        conditions = [
            (self.df["Compound"] == "DMSO"),
            (self.df["Compound"] == "DMSO+FSK"),
            (self.df["Compound"].str.lower().isin(self.control_treatments)),
            (
                (self.df["Compound"].str.lower().isin(self.toxic_treatments))
                # & (self.df["Concentration"] >= 1)
            ),
        ]
        choices = [
            "neg_control",
            "pos_control",
            "treatment_control",
            "treatment_toxic",
        ]
        self.df["TreatmentType"] = np.select(conditions, choices, default="treatment")

        newdf = (
            self.df.select_dtypes(include=numerics)
            .assign(
                Compound=self.df["Compound"],
                TreatmentType=self.df["TreatmentType"],
                PlateID=self.df["plate.folder"],
                PlateRow=self.df["row.char"],
                PlateColumn=self.df["column"],
                QC=self.df["QC"],
            )
            .drop(
                columns=["column"],
                # `column` dropped as it's numerical and remains after select_dtypes
            )
        )
        # Fixing column order, changing PlateID to numerical and resetting index
        self.df = (
            newdf[["Compound"] + self.id_cols + newdf.columns.tolist()[1:-6]]
            .assign(PlateID=newdf["PlateID"].astype(int))
            .sort_values(by=["PlateID", "PlateRow", "PlateColumn"])
            .reset_index(drop=True)
        )
        return


if __name__ == "__main__":
    rd_shut_the_hell_up()
    print("Processing selleck_chem data")
    selchem = selleck_chem()
    selchem.drop_cols()
    selchem.update_stimulation()
    selchem.update_id_attributes()
    selchem.update_features()
    selchem.save_data(which="df")
    selchem.npi_normalize(method="median")
    selchem.save_data(which="norm_df")
    selchem.z_score_normalize(method="median")
    selchem.save_data(which="norm_df")
    selleck_df = selchem.get_chem_structs_df()
    struct_dict = selchem.structure_mapping(selleck_df, n_jobs=5)
    antineo_drugs = selchem.get_antineoplastic_from_chembl()
    selchem.flag_antineoplastic_compounds("both", antineo_drugs, 1, n_jobs=5)
    selchem.z_prime_statistics(value_var="all", method="median")
    selchem.save_data("stats_df")

    print("Processing spectrum data")
    spectr = spectrum()
    spectr.drop_cols()
    spectr.update_id_attributes()
    spectr.update_features()
    spectr.save_data(which="df")
    spectr.npi_normalize(method="median")
    spectr.save_data(which="norm_df")
    spectr.z_score_normalize(method="median")
    spectr.save_data(which="norm_df")
    spectrum_df = spectr.get_chem_structs_df()
    struct_dict = spectr.structure_mapping(spectrum_df)
    # no need to retrieve antineo_drugs again
    selchem.flag_antineoplastic_compounds("both", antineo_drugs, 1, n_jobs=5)
    spectr.z_prime_statistics(value_var="all", method="median")
    spectr.save_data("stats_df")
    print(" ### DONE ###")
