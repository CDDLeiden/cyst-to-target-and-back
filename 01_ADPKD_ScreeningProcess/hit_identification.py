"""Module to identifiy hits on the ADPKD screening datasets. Though different methods are described, hit identification 
is performed on the NPI-normalized dataset, using the median of the per-well mean(cystic spheroid size)
distributions ± the mutual 1.5 * MAD (median absolute deviation) distance from the distributions

Treatment distribution (TD)
DMSO+FSK distribution (CScontrol)

Cyst swelling reducing compounds:
>>> median(TD) + 1.5 * MAD(TD) < median(CScontrol) - 1.5 * MAD(CScontrol)

Cyst swelling enhancing compounds:
>>> median(TD) - 1.5 * MAD(TD) > median(CScontrol) + 1.5 * MAD(CScontrol)
"""

import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from data_prepare import get_time_string
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy import stats
from scipy.stats import median_abs_deviation as mad
from screening_process import selleck_chem, spectrum, spectrum_validation
from tqdm import tqdm
from venn import venn


def smi_to_connectivity(smi):
    mol = Chem.MolFromSmiles(smi)
    connectivity = Chem.MolToInchiKey(mol).split("-")[0]
    return connectivity


def chiral_ECFP4(mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=4096, useChirality=True)


def ECFP4(mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=4096, useChirality=False)


def flatten_lists(list_of_lists):
    return [item for sublist in list_of_lists for item in sublist]


def plot_hits_per_plate(dataset_obj, distance=None):
    # plotting the number of hits per plate
    df = dataset_obj.norm_df.copy().sort_values(["PlateID", "PlateRow", "PlateColumn"]).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(8, 4))

    if dataset_obj.normalization_type in ["z-score", "npi"]:
        if dataset_obj.normalization_method == "median":
            criteria = f"treatment median < DMSO+FSK median - {distance} * MAD"
        else:
            criteria = f"treatment mean < DMSO+FSK mean - {distance} * SD"

    if dataset_obj.normalization_type == "z-score":
        score_col = "ZScore_hitflag != 'inactive'"
    elif dataset_obj.normalization_type == "npi":
        score_col = "NPIScore_hitflag != 'inactive'"
    elif dataset_obj.normalization_type == "b-score":
        score_col = "BScore_hitflag != 'inactive'"
        criteria = "lower 10% of B-Score distribution"
    (df.query(score_col).drop_duplicates(["Identifier"]).groupby("PlateID").count()["Identifier"].plot.bar())
    ax.set_title(
        f"{dataset_obj.name.capitalize()} dataset\n" "Number of hits per plate\n" f"Criteria: {criteria}"
    )
    fig.savefig(
        figures_root / f"{dataset_obj.name}_{dataset_obj.normalization_type}_hits_per_plate.png",
        dpi=150,
        facecolor="w",
        bbox_inches="tight",
    )
    plt.close(fig)
    return df


def add_hit_thresholds_toplot(
    reducer_threshold: float, enhancer_threshold: float, ax: plt.Axes, orientation: str
):
    """
    Function that takes as input a matplotlib axis object and adds
    the desired thresholds for the hit identification. Reducer threshold
    will be plotted red and enhancer threshold blue.

    Args:
        ax: matplotlib.pyplot axis object
        orientation: str, either "horizontal" or "vertical"
        reducer_threshold: float
        enhancer_threshold: float
    """
    # Setting horizontal and vertical lines
    ylims = ax.get_ylim()
    xlims = ax.get_xlim()
    if orientation == "vertical":
        ax.vlines(
            reducer_threshold,
            ylims[0],
            ylims[1],
            color="red",
            alpha=0.5,
            linestyle="--",
        )
        ax.vlines(
            enhancer_threshold,
            ylims[0],
            ylims[1],
            color="blue",
            alpha=0.5,
            linestyle="--",
        )
    elif orientation == "horizontal":
        ax.hlines(
            reducer_threshold,
            xlims[0],
            xlims[1],
            color="red",
            alpha=0.5,
            linestyle="--",
        )
        ax.hlines(
            enhancer_threshold,
            xlims[0],
            xlims[1],
            color="blue",
            alpha=0.5,
            linestyle="--",
        )
    ax.set_ylim(ylims)
    ax.set_xlim(xlims)
    return ax


def extract_identifiers_set(df, hitflag_col, flag):
    """
    Avoid code repetition with extracting the set of identifiers
    for the hit compounds.
    """
    return set(df.query(f"{hitflag_col} == '{flag}'").drop_duplicates(["Identifier"])["Identifier"].unique())


def comp_name_only(list_of_identifiers):
    """
    Function to extract only the compound name from the identifier.
    """
    return ["~".join(i.split("~")[:1]) for i in list_of_identifiers]


def add_screening_name(set_compnames, dataset_obj):
    """
    Function to add the screening name to the compound name.
    """
    return [hitname + f"~{dataset_obj.name}" for hitname in set_compnames]


def get_inf_cols(df: pd.DataFrame, verbose: bool = True):
    """
    Function to return column names with inf values.

    Args:
        df: pd.DataFrame to check for inf values.
        verbose: Print column names. Defaults to True.

    Returns:
        columns with inf values (list)
    """
    inf_cols = df.columns[df.isin([np.inf, -np.inf]).any()].tolist()
    if verbose:
        print(f"Columns with inf values: {inf_cols}")
    return inf_cols


def get_nan_cols(df: pd.DataFrame, threshold: int, verbose: bool = True):
    """
    Function to return column names with nan values.

    Args:
        df: pd.DataFrame to check for nan values.
        threshold: amount of nan values allowed in a column.
        verbose: prints column names. Defaults to True.

    Returns:
        columns with nan values (list)
    """
    nan_cols = df.columns[np.where(np.array([c for c in df.isnull().sum()]) > threshold)].tolist()
    if verbose:
        print(f"Columns with nan values: {nan_cols}")
    return nan_cols


selchem = selleck_chem()
spectr = spectrum()
specval = spectrum_validation()

datasets = [selchem, spectr, specval]

hit_dict = {}
root_dir = Path(__file__).parents[1]
dataset_root = root_dir / "data/adpkd_screening/screening_data/"
figures_root = root_dir / ("figures/hit_analysis")
if not figures_root.exists():
    figures_root.mkdir(parents=True)

# B-score hitflagging takes lowest % of the distribution
b_score_hit_rate = 0.1
# Distance from the median / mean of the distribution
zhit_distance = 1.5
npihit_distance = 1.5
# Whether use median & MAD or mean & SD
zhit_method = "median"
npihit_method = "median"
# Type of chemical structures to map
smiles_type = "as-is"
# If desired to compare the top % hits of each normalized dataset
top_compare = True
# If desired to save the final dataset with the identified hits
wanna_save = True

all_medians_datasets = []
all_hits_by_identifier = {
    "reducers": {
        "zscore": [],
        "npi": [],
        "bscore": [],
    },
    "enhancs": {
        "zscore": [],
        "npi": [],
        "bscore": [],
    },
}

# dataset_idx will be used to access z_scored_datasets
for dataset_idx, data in enumerate(datasets):
    data.drop_cols()
    # The only dataset requiring this function (TODO: change this)
    if data.name == "selleckchem":
        data.update_stimulation()
    data.update_id_attributes()
    data.update_features()
    data.flag_outliers()
    data.drop_bad_qc(which="df")

    # NPI normalization and hitflagging
    data.npi_normalize()
    data.npi_hitflagging(method=npihit_method, distance=npihit_distance)

    df = plot_hits_per_plate(data, distance=npihit_distance)
    if top_compare:
        data.top_percent_hitpicking(hit_rate=b_score_hit_rate, perplate=False)
        npi_reducers = extract_identifiers_set(data.norm_df, "Top_percent", "reducer")
        npi_enhancs = extract_identifiers_set(data.norm_df, "Top_percent", "enhanc")
    else:
        npi_reducers = extract_identifiers_set(df, "NPIScore_hitflag", "reducer")
        npi_enhancs = extract_identifiers_set(df, "NPIScore_hitflag", "enhanc")

    npi_series = (
        df[["Identifier", "NPIScore_hitflag", "obj.Mean(area)"]]
        .rename(columns={"obj.Mean(area)": "npiscored_obj_Mean_area"})
        .copy()
    )

    # Z-score normalization and hitflagging
    data.z_score_normalize(as_booij=False)
    data.zscore_hitflagging(method=zhit_method, distance=zhit_distance)

    if top_compare:
        data.top_percent_hitpicking(hit_rate=b_score_hit_rate, perplate=False)
        zscore_reducers = extract_identifiers_set(data.norm_df, "Top_percent", "reducer")
        zscore_enhancs = extract_identifiers_set(data.norm_df, "Top_percent", "enhanc")
    else:
        zscore_reducers = extract_identifiers_set(df, "ZScore_hitflag", "reducer")
        zscore_enhancs = extract_identifiers_set(df, "ZScore_hitflag", "enhanc")

    # Implementing chemical structures:
    chemstruct_df = data.get_chem_structs_df(verbose=False)
    data.structure_mapping(chemstruct_df, smiles_type=smiles_type, verbose=False)
    antineo_drugs = data.get_antineoplastic_from_chembl()
    data.flag_antineoplastic_compounds("both", antineo_drugs, 1, n_jobs=5)

    df = plot_hits_per_plate(data, distance=zhit_distance)
    assert "SMILES" in df.columns, "SMILES column not found in df."

    z_score_series = (
        df[["Identifier", "ZScore_hitflag", "obj.Mean(area)"]]
        .rename(columns={"obj.Mean(area)": "zscored_obj_Mean_area"})
        .copy()
    )

    # B-score hitflagging
    data.b_score_normalize(1, 0.01, which="df")
    data.bscore_hitflagging(hit_rate=b_score_hit_rate, perplate=False)
    data.norm_df = data.norm_df.sort_values(["PlateID", "PlateRow", "PlateColumn"]).reset_index(drop=True)

    plot_hits_per_plate(data)

    bscore_reducers = extract_identifiers_set(data.norm_df, "BScore_hitflag", "reducer")
    bscore_enhancs = extract_identifiers_set(data.norm_df, "BScore_hitflag", "enhanc")
    b_score_series = (
        data.norm_df[["Identifier", "BScore_hitflag", "obj.Mean(area)"]]
        .rename(columns={"obj.Mean(area)": "bscored_obj_Mean_area"})
        .copy()
    )

    assert all(df["Identifier"] == z_score_series["Identifier"])

    rename_dict = {
        "zscored_obj_Mean_area": f"zscored_obj_Mean_area_{zhit_method}",
        "npiscored_obj_Mean_area": f"npiscored_obj_Mean_area_{npihit_method}",
        "bscored_obj_Mean_area": "bscored_obj_Mean_area_median",
    }
    df = df.assign(
        ZScore_hitflag=z_score_series["ZScore_hitflag"],
        zscored_obj_Mean_area=z_score_series["zscored_obj_Mean_area"],
        NPIScore_hitflag=npi_series["NPIScore_hitflag"],
        npiscored_obj_Mean_area=npi_series["npiscored_obj_Mean_area"],
        BScore_hitflag=b_score_series["BScore_hitflag"],
        bscored_obj_Mean_area=b_score_series["bscored_obj_Mean_area"],
        Screening=lambda x: [data.name] * len(x),
    )
    all_medians_datasets.append(df)

    bscore_reduc_thresh = np.percentile(df["bscored_obj_Mean_area"], b_score_hit_rate * 100)
    bscore_enhanc_thresh = np.percentile(df["bscored_obj_Mean_area"], 100 - b_score_hit_rate * 100)

    zscore_reduc_allplates = comp_name_only(zscore_reducers)
    npi_reduc_allplates = comp_name_only(npi_reducers)
    bscore_reduc_allplates = comp_name_only(bscore_reducers)
    zscore_reducers = add_screening_name(zscore_reduc_allplates, data)
    npi_reducers = add_screening_name(npi_reduc_allplates, data)
    bscore_reducers = add_screening_name(bscore_reduc_allplates, data)

    zscore_enhanc_allplates = comp_name_only(zscore_enhancs)
    npi_enhanc_allplates = comp_name_only(npi_enhancs)
    bscore_enhanc_allplates = comp_name_only(bscore_enhancs)
    zscore_enhancs = add_screening_name(zscore_enhanc_allplates, data)
    npi_enhancs = add_screening_name(npi_enhanc_allplates, data)
    bscore_enhancs = add_screening_name(bscore_enhanc_allplates, data)

    all_hits_by_identifier["reducers"]["zscore"].extend(zscore_reducers)
    all_hits_by_identifier["reducers"]["npi"].extend(npi_reducers)
    all_hits_by_identifier["reducers"]["bscore"].extend(bscore_reducers)
    all_hits_by_identifier["enhancs"]["zscore"].extend(zscore_enhancs)
    all_hits_by_identifier["enhancs"]["npi"].extend(npi_enhancs)
    all_hits_by_identifier["enhancs"]["bscore"].extend(bscore_enhancs)

    zscore_hits_allplates = set(zscore_enhanc_allplates + zscore_reduc_allplates)
    npi_hits_allplates = set(npi_enhanc_allplates + npi_reduc_allplates)
    bscore_hits_allplates = set(bscore_enhanc_allplates + bscore_reduc_allplates)

    fig, ax = plt.subplots(figsize=(6, 6))
    tovenn_hits = {
        "Z-Score": set(zscore_hits_allplates),
        "B-Score": set(bscore_hits_allplates),
        "NPI": set(npi_hits_allplates),
    }
    venn(tovenn_hits, legend_loc="lower left", ax=ax)
    fig.savefig(
        figures_root / f"{data.name}_venn_hits.png",
        dpi=150,
        facecolor="w",
        bbox_inches="tight",
    )
    fig.suptitle(f"{data.name.capitalize()}: Venn diagram of hits")
    plt.close(fig)

    print("Generating figures per plate...")
    for plate in tqdm(df.PlateID.unique()):
        sub_df = df.query("PlateID == @plate").copy()
        poscontrol_median = sub_df.query("TreatmentType == 'pos_control'")["zscored_obj_Mean_area"].median()
        poscontrol_mad = mad(sub_df.query("TreatmentType == 'pos_control'")["zscored_obj_Mean_area"].values)
        # The threshold is more strict for the z-score hits
        zscore_enhanc_thresh = poscontrol_median + zhit_distance * poscontrol_mad * 2
        zscore_reduc_thresh = poscontrol_median - zhit_distance * poscontrol_mad
        newtreat_label = [
            "z-score_hit" if condition else label
            for label, condition, in zip(
                sub_df["TreatmentType"],
                (sub_df["ZScore_hitflag"] != "inactive") & (sub_df["TreatmentType"] == "treatment"),
            )
        ]
        # Scatter plot of the b-scores and the z-scores
        fig, ax = plt.subplots(figsize=(4, 4))
        ax = sns.scatterplot(
            data=(
                sub_df.assign(PlateColumn=lambda x: x["PlateColumn"].astype(str))
                .assign(TreatmentType=newtreat_label)
                .sort_values("TreatmentType")
            ),
            x="zscored_obj_Mean_area",
            y="bscored_obj_Mean_area",
            hue="TreatmentType",
            palette="tab10",
            ax=ax,
        )
        # Setting horizontal and vertical lines
        ax = add_hit_thresholds_toplot(
            zscore_reduc_thresh, zscore_enhanc_thresh, ax=ax, orientation="vertical"
        )
        ax = add_hit_thresholds_toplot(
            bscore_reduc_thresh, bscore_enhanc_thresh, ax=ax, orientation="horizontal"
        )
        ax.legend(bbox_to_anchor=(1.04, 1), borderaxespad=0)
        ax.set_ylabel("B-scored obj.Mean(area)")
        ax.set_xlabel("Z-scored obj.Mean(area)")
        ax.set_title(
            f"{data.name.capitalize()} dataset\n" f"Plate {plate}: B-scored vs Z-scored obj.Mean(area)"
        )
        fig.savefig(
            figures_root / f"{data.name}_Plate{plate}_scatterplot_B_and_Z-scores.png",
            dpi=150,
            facecolor="w",
            bbox_inches="tight",
        )
        plt.close(fig)

        # jointplot of the dataset
        grid = sns.jointplot(
            data=sub_df.query("TreatmentType == 'treatment'"),
            x="zscored_obj_Mean_area",
            y="bscored_obj_Mean_area",
            kind="reg",
            ax=ax,
            height=6,
        )

        grid.ax_joint = add_hit_thresholds_toplot(
            zscore_reduc_thresh,
            zscore_enhanc_thresh,
            ax=grid.ax_joint,
            orientation="vertical",
        )
        grid.ax_joint = add_hit_thresholds_toplot(
            bscore_reduc_thresh,
            bscore_enhanc_thresh,
            ax=grid.ax_joint,
            orientation="horizontal",
        )

        pearson_corr = stats.pearsonr(
            df.query(f"PlateID == {plate}")["zscored_obj_Mean_area"],
            df.query(f"PlateID == {plate}")["bscored_obj_Mean_area"],
        )

        plt.suptitle(
            f"{data.name.capitalize()} treatment conditions\n"
            f"Plate {plate}: person correlation = {pearson_corr[0]:.2f};"
            f" p-value={pearson_corr[1]:.4f}",
            y=1.05,
        )
        grid.figure.savefig(
            figures_root / f"{data.name}_Plate{plate}_jointplot_B_and_Z-scores.png",
            dpi=150,
            facecolor="w",
            bbox_inches="tight",
        )
        plt.close(grid.figure)

timestring = get_time_string()
dev_pattern = re.compile("\.SD")

all_medians_df = pd.concat(all_medians_datasets)
nan_cols = get_nan_cols(all_medians_df, 50, verbose=True)
inf_cols = get_inf_cols(all_medians_df, verbose=True)
undesired_cols = [c for c in all_medians_df.columns if dev_pattern.findall(c)]
all_medians_df = (
    all_medians_df.drop(columns=nan_cols + inf_cols + undesired_cols + ["Median_cystsize"])
    .assign(Identifier=lambda x: x["Identifier"] + "~" + x["Screening"])
    .rename(columns=(rename_dict))
    .dropna(axis=0)
    .reset_index(drop=True)
)

# Keep only unique identifiers
all_hits_by_identifier["enhancs"] = {k: list(set(v)) for k, v in all_hits_by_identifier["enhancs"].items()}
all_hits_by_identifier["reducers"] = {k: list(set(v)) for k, v in all_hits_by_identifier["reducers"].items()}
# json with all the hit compounds
json_filepath = figures_root / "alldatasets_hits.json"
with json_filepath.open("w") as f:
    json.dump(all_hits_by_identifier, f, indent=2)

# Saving the updated dataframe
all_medians_df.to_csv(
    figures_root / f"All_datasets_z-norm_hitflagged_{timestring}.csv.gz",
    index=False,
    compression="gzip",
)

# Summarizing bioactivity information for modelling:
unique_comps = (
    all_medians_df.copy()
    .query("Concentration == 1")
    .drop_duplicates("Compound")
    .sort_values("Compound")
    .set_index("Compound")
)
data_modelling = (
    all_medians_df.copy()
    .query("Concentration == 1")
    .groupby(["Compound"])["NPIScore_hitflag"]
    .value_counts()
    .unstack()
    .fillna(0)
    .assign(Activity=lambda x: x[["enhancer", "inactive", "reducer"]].idxmax(axis=1))
    .assign(
        Certainty=lambda x: x[["enhancer", "inactive", "reducer"]].max(axis=1)
        / x[["enhancer", "inactive", "reducer"]].sum(axis=1)
    )
    .drop(columns=["enhancer", "inactive", "reducer"])
    .sort_values("Compound")
    .assign(
        SMILES=unique_comps["SMILES"],
        Antineoplastic=unique_comps["Antineoplastic"],
        Screening=unique_comps["Screening"],
    )
)
assert all(unique_comps.index == data_modelling.index), "Index mismatch; Compounds are not the same"
data_modelling = data_modelling.reset_index().assign(
    Name_Activity=lambda x: x["Compound"] + "\n" + x["Activity"]
)
# data_modelling['Certainty'].plot.hist(ylim=(0,60), bins=30)
# data_modelling['Activity'].value_counts().plot.bar()
comps_to_rm = [
    "3-hydroxytyramine",
    "Modafinil",
    "Saracatinib",
    "Bez-235",
    "Carvedilol phosphate",
    "Decamethonium bromide",
    "Deoxysappanone b 7.3'-dimethyl ether acetate",
    "Deoxysappanone b 7.4'-dimethyl ether",
    "Desloratidine",
    "Osi-420",
    "Dovitinib Dilactic acid",
    "Edetate disodium",
    "Erythromycin stearate",
    "Solanesol",
    "Fomepizole hydrochloride",
    "Gossypol-acetic acid complex",
    "Imatinib Mesylate",
    "Theophylline",
    "Trihexyphenidyl hydrochloride",
    "R935788",
    "Tofacitinib citrate",
]

todrop_idxs = data_modelling.query("Compound in @comps_to_rm").index
print("dropped the following compounds:")
print(data_modelling.loc[todrop_idxs].Compound.values)
data_modelling = (
    data_modelling.drop(index=todrop_idxs)
    .drop(columns=["Name_Activity"])
    .assign(Connectivity=lambda x: x.SMILES.apply(smi_to_connectivity))
    .reset_index(drop=True)
)
if npihit_distance == 1.5:
    distance_metric = "Default"

hit_id_dir =  (dataset_root.parent / f"identified_hits").exists():
if not hit_id_dir.exists():
    hit_id_dir.mkdir()

if wanna_save:
    data_modelling.to_csv(
        hit_id_dir / f"pkd_HitCompounds_NPI-{npihit_method}"
        f"-{distance_metric}Distance-hitflag_{smiles_type}SMILES"
        f"_{timestring}.csv",
        index=False,
    )
