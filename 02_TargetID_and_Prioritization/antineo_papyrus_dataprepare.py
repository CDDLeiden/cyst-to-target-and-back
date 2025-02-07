import json
import timeit
from datetime import datetime
from pathlib import Path

import pandas as pd
from papyrus_scripts.preprocess import (
    consume_chunks,
    keep_contains,
    keep_organism,
    keep_quality,
    keep_similar,
)
from papyrus_scripts.reader import read_papyrus, read_protein_set
from papyrus_scripts.subsim_search import FPSubSim2
from rdkit import Chem, RDLogger


def get_time_string():
    return datetime.now().strftime("%Y%m%d-%H%M%S")


def rd_shut_the_hell_up():
    """Make the RDKit be a bit more quiet
    @return: None
    """
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)


def smi_to_connectivity(smi):
    mol = Chem.MolFromSmiles(smi)
    connectivity = Chem.MolToInchiKey(mol).split("-")[0]
    return connectivity


rd_shut_the_hell_up()


def connect_to_smi_mapping(
    papyrus_data: pd.DataFrame,
):
    """
    Function written to
    """
    df = papyrus_data.copy().drop_duplicates(subset=["SMILES", "connectivity"])
    id_series = (
        df.assign(ALLSMILES=lambda x: " " + x["SMILES"].astype(str))
        .groupby("connectivity")["ALLSMILES"]
        .sum()
    )
    comp_dict = id_series.str.split(" ")
    comp_dict = {k: list(set([v for v in value if v != ""])) for k, value in comp_dict.items()}
    return comp_dict


def get_available_fingerprint(fpsubsim2_file: str) -> list:
    """
    Function that loads the fpsubsim file and returns the registered fingerprint
    name inside the file. User as parameter for keep_similar(fingerprint=)
    """
    fpss2 = FPSubSim2()
    fpss2.load(fpsubsim2_file)
    return list(fpss2.available_fingerprints.keys())


def get_allantineoplastic_papyrus(
    savepath: str, method: str, simi_threshold: float = 1.0, fp_type: str = None
) -> pd.DataFrame:
    """
    Params:
    savepath -> path for saving the filtered papyrus dataset by the molecules present
                in the antineo_df.
    method -> either `connectivity`, `chembl_id` or `fp_similarity`, determining the `papyrus_scripts`
              filtering method.
    simi_threshold -> Float value between 0 and 1 used for the tanimoto similarity selection.
                      Defaults to `1.0`.
    fp_type -> file created (<filename.h5>) with fpsubsim2 (see `papyrus fpsubsim2 --fhelp`) with
               the pre-computed fingerprint. See path `/home/davidararipe/.data/papyrus/05.5/fingerprints`.
    """
    if isinstance(savepath, str):
        savepath = Path(savepath)
        assert savepath.exists(), f"{savepath} does not exist."

    available_methods = ["connectivity", "fp_similarity", "chembl_id"]
    if method not in available_methods:
        raise AttributeError(
            "method must be one of the following:" "`connectivity`, `fp_similarity`, `chembl_id`"
        )
    # Path to the antineoplastic compounds from DataPrepare.get_antineoplastic_from_chembl()
    antineo_dictpath = Path("data/known_antineoplastic.json")
    PAPYRUS_DATA_FOLDER = Path("/home/davidararipe/.data/")
    fpsubsim2_path = PAPYRUS_DATA_FOLDER / "papyrus/05.5/descriptors/"
    fingerprint_attribute = None

    # if method == "fp_similarity":
    #     # TODO: This is not implemented yet, I need to test it first
    #     # Check if the chosen fingerprint is available
    #     avail_fps = [path.split("/")[-1] for path in glob(f"{fpsubsim2_path}*.h5")]

    #     if fp_type not in avail_fps:
    #         raise AttributeError(f"fp_type not available. Try one from {avail_fps}")
    #     else:
    #         fpsubsim2_path = os.path.join(fpsubsim2_path, fp_type)

    #     # See what is the fingerprint name inside the h5 file.
    #     fingerprint_attribute = get_available_fingerprint(
    #         fpsubsim2_file=fpsubsim2_path
    #     )[0]

    time_string = get_time_string()

    with open(antineo_dictpath, "r") as json_file:
        antineo_dict = json.load(json_file)

    # Converting the dictionary to a pandas dataframe to use DataPrepare.get_connectivity_dict() later
    antineo_df = (
        pd.DataFrame.from_dict(antineo_dict, orient="index")
        .reset_index()
        .rename(columns={"index": "Compound"})
        .assign(Source=lambda x: ["ChEMBL"] * len(x))
    )
    protein_data = read_protein_set(source_path=PAPYRUS_DATA_FOLDER)
    query_connect = list()

    sample_data = read_papyrus(is3d=False, chunksize=100000, source_path=PAPYRUS_DATA_FOLDER)
    filter1 = keep_organism(
        data=sample_data,
        protein_data=protein_data,
        organism=["Homo sapiens", "Mus musculus", "Rattus norvegicus"],
        generic_regex=True,
    )
    filter2 = keep_quality(data=filter1, min_quality="medium")
    if method == "fp_similarity":
        smiles = antineo_df["iso_smiles"].values
        csv_fname = savepath / f"antineoplastic_papyrus_{method}_{fp_type.split('.')[0]}_{time_string}.csv"
        filter3 = keep_similar(
            data=filter2,
            molecule_smiles=smiles,
            fpsubsim2_file=fpsubsim2_path,
            fingerprint=fingerprint_attribute,
            threshold=simi_threshold,
            cuda=False,
        )
    elif method == "connectivity":
        if "Connectivity" not in antineo_df.columns:
            antineo_df = antineo_df.assign(
                Connectivity=lambda x: antineo_df["PapyrusSmiles"].apply(smi_to_connectivity)
            )
        query_connect = antineo_df["Connectivity"].unique()
        csv_fname = savepath / f"antineoplastic_papyrus_{method}_{time_string}.csv"
        pattern = "|".join(query_connect)
        filter3 = keep_contains(data=filter2, column="connectivity", value=pattern, case=True, regex=True)

    elif method == "chembl_id":
        chembl_ids = antineo_df["chembl_id"].unique()
        csv_fname = savepath / f"antineoplastic_papyrus_{method}_{time_string}.csv"
        pattern = "|".join(chembl_ids)
        filter3 = keep_contains(data=filter2, column="CID", value=pattern, case=True, regex=True)

    filtered_data = consume_chunks(filter3, progress=True, total=13)
    filtered_data.to_csv(csv_fname, index=False)

    number_unique_mols = len(filtered_data["connectivity"].unique())
    print("Size of the Papyrus dataframe:", len(filtered_data))
    print(f"method: {method} | fp_type: {fp_type}")
    print(f"unique connectivities: {number_unique_mols}")
    print(f"unique organisms: {filtered_data['Organism'].unique()}")
    print("done!")

    comp_dict = connect_to_smi_mapping(filtered_data)
    if all([fp_type is not None, method != "chembl_id"]):
        method = fp_type.split(".")[0]
        json_path = savepath / f"{fp_type}_connect_smi.json"
    elif method == "chembl_id":
        json_path = savepath / "chemblID_connect_smi.json"
    else:
        json_path = savepath / "connect_smi.json"

    with json_path.open("w") as json_file:
        json.dump(comp_dict, json_file, indent=4)
    return filtered_data


if __name__ == "__main__":
    start = timeit.default_timer()
    available_fingerprints = [
        "ecfp4_chiral.h5",
    ]
    get_allantineoplastic_papyrus(
        savepath="data/2_antineo_papyrus_data/",
        method="connectivity",
    )
    # get_allantineoplastic_papyrus(
    #     savepath="/home/davidararipe/databases/pkd_modelling/2_antineo_papyrus_data/",
    #     method="fp_similarity",
    #     simi_threshold=1.0,
    #     fp_type="rdkpatternfingerprint.h5",
    # )
    stop = timeit.default_timer()
    execution_time = stop - start
    print("Program Executed in " + str(execution_time))  # It returns time in seconds
