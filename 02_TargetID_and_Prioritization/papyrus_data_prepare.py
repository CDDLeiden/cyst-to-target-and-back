import json
import timeit
from pathlib import Path

import pandas as pd
from antineo_papyrus_dataprepare import (
    connect_to_smi_mapping,
    get_available_fingerprint,
    get_time_string,
    rd_shut_the_hell_up,
    smi_to_connectivity,
)
from papyrus_scripts import download_papyrus
from papyrus_scripts.preprocess import (
    consume_chunks,
    keep_contains,
    keep_organism,
    keep_quality,
    keep_similar,
)
from papyrus_scripts.reader import read_papyrus, read_protein_set


def papyrus_data_prepare(
    screening_df: pd.DataFrame,
    savepath: str,
    method: str,
    simi_threshold: float = 1.0,
    fp_type: str = None,
    is3d=False,
) -> pd.DataFrame:
    """
    Function to query the papyrus dataset with the molecules that I have
    on my pkd dataset. This is done by either using the connectivity
    of the molecules, or by fetching molecules within a certain tanimoto
    similarity threshold.

    Args:
        screening_df: PKD bioactivity dataset
            [file generated from hit_identification.py]
        savepath: path for saving the filtered papyrus dataset.
        method: either `connectivity` or `fp_similarity` -> to filer papyrus dataset.
        simi_threshold: When using the `fp_similarity`. Defaults to 1.0.
        fp_type: file created (<filename.h5>) with fpsubsim2
            (see `papyrus fpsubsim2 --fhelp`) with the pre-computed fingerprint.
            See path `papyrus/05.6/fingerprints`. Defaults to None.
        is3d: Whether to use stereochemistry. Defaults to False.

    Raises:
        AttributeError: when `method` is not `connectivity` or `fp_similarity`
        AttributeError: When `fp_type` not available.

    Returns:
        filtered_data -> Filtered Papyrus dataset.
    """
    if isinstance(savepath, str):
        savepath = Path(savepath)
        assert savepath.exists(), f"{savepath} does not exist."

    available_methods = ["connectivity", "fp_similarity"]
    if method not in available_methods:
        raise AttributeError("method must be either `connectivity` or `fp_similarity`")

    # Path to the FP similarity database file:
    PAPYRUS_DATA_FOLDER = Path(__file__).parents[1] / "data/papyrus_data"
    download_papyrus(outdir=str(PAPYRUS_DATA_FOLDER), version="05.5", nostereo=not is3d, progress=True)
    fpsubsim2_path = PAPYRUS_DATA_FOLDER / "papyrus/05.5/descriptors/"
    fingerprint_attribute = None
    time_string = get_time_string()
    protein_data = read_protein_set(source_path=PAPYRUS_DATA_FOLDER)

    if method == "connectivity":
        if "Connectivity" not in screening_df.columns:
            screening_df = screening_df.assign(
                Connectivity=lambda x: screening_df["standardised_smiles"].apply(smi_to_connectivity)
            )
        query_connect = screening_df["Connectivity"].unique()
        print(query_connect)

    elif method == "fp_similarity":
        avail_fps = sorted([str(fpath).split("/")[-1] for fpath in fpsubsim2_path.glob("*.h5")])
        if fp_type not in avail_fps:
            raise AttributeError(f"fp_type not available. Try one from {avail_fps}")
        else:
            fpsubsim2_path = fpsubsim2_path / fp_type

        fingerprint_attribute = get_available_fingerprint(fpsubsim2_file=fpsubsim2_path)[0]
        print(f"Using fingerprint: {fingerprint_attribute}")
        query_smiles = screening_df["SMILES"]

    sample_data = read_papyrus(is3d=is3d, chunksize=100000, source_path=PAPYRUS_DATA_FOLDER)
    filter1 = keep_organism(
        data=sample_data,
        protein_data=protein_data,
        organism=["Homo sapiens", "Mus musculus", "Rattus norvegicus"],
        generic_regex=True,
    )
    filter2 = keep_quality(data=filter1, min_quality="medium")
    if method == "fp_similarity":
        if is3d:
            csv_fname = (
                f"papyrus_ADPKDbioactivity_{method}_{fp_type.split('.')[0]}"
                f"_3D_simi{simi_threshold*100}_{time_string}.csv"
            )
        else:
            csv_fname = (
                f"papyrus_ADPKDbioactivity_{method}_{fp_type.split('.')[0]}"
                f"_simi{simi_threshold*100}_{time_string}.csv"
            )
        filter3 = keep_similar(
            data=filter2,
            molecule_smiles=query_smiles,
            fpsubsim2_file=fpsubsim2_path,
            fingerprint=fingerprint_attribute,
            threshold=simi_threshold,
            cuda=False,
        )
    elif method == "connectivity":
        csv_fname = f"papyrus_ADPKDbioactivity_{method}_{time_string}.csv"
        pattern = "|".join(query_connect)
        filter3 = keep_contains(data=filter2, column="connectivity", value=pattern, case=True, regex=True)

    filtered_data = consume_chunks(filter3, progress=True, total=13)
    number_unique_mols = len(filtered_data["connectivity"].unique())
    print("Size of the Papyrus dataframe:", len(filtered_data))
    print(f"method: {method} | fp_type: {fp_type}")
    print(f"unique connectivities: {number_unique_mols}")
    print(f"unique organisms: {filtered_data['Organism'].unique()}")
    print("done!")

    comp_dict = connect_to_smi_mapping(filtered_data)
    if fp_type is not None:
        method = fp_type.split(".")[0]
        json_path = savepath / f"{fp_type}_connect_smi.json"
    else:
        json_path = savepath / "connect_smi.json"

    with json_path.open("w") as json_file:
        json.dump(comp_dict, json_file, indent=4)
    filtered_data.to_csv(savepath / csv_fname, index=False)
    return filtered_data


root_dir = Path(__file__).parents[1]

if __name__ == "__main__":
    rd_shut_the_hell_up()
    start = timeit.default_timer()
    available_fingerprints = [
        "ecfp4_chiral.h5",
    ]
    savepath = root_dir / "data/papyrus_data"
    if not savepath.exists():
        savepath.mkdir(parents=True)

    hit_data_root = root_dir / "data/adpkd_screening/identified_hits"
    file_path = sorted(
        hit_data_root.glob(
            "*HitCompounds_NPI-median-DefaultDistance-hitflag_as-isSMILES*"
        )  # File output from hit_identification.py
    )[-1]
    screening_df = pd.read_csv(file_path)
    # df = papyrus_data_prepare(
    #     screening_df,
    #     savepath=savepath,
    #     method="fp_similarity",
    #     simi_threshold=1,
    #     fp_type="ecfp4_chiral.h5",
    #     is3d=True,
    # )
    df = papyrus_data_prepare(screening_df, savepath=savepath, method="connectivity")
    stop = timeit.default_timer()
    execution_time = stop - start
    print("Program Executed in " + str(execution_time))
