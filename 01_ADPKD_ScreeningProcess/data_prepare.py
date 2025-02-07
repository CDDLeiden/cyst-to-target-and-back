import codecs
import os
import re
from datetime import datetime
from functools import partial
from itertools import combinations, compress, product
from multiprocessing import Pool
from operator import itemgetter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import PIL
import pubchempy as pcp
import seaborn as sns
from chembl_structure_pipeline import standardizer
from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem, RDConfig, RDLogger
from rdkit.Chem import AllChem, DataStructs, Draw, PandasTools, rdMolDescriptors
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.SimDivFilters import rdSimDivPickers
from tqdm import tqdm

try:
    from chembl_webresource_client.new_client import new_client
except:
    pass


def get_time_string():
    return datetime.now().strftime("%Y%m%d-%H%M%S")


def rd_shut_the_hell_up():
    """Make the RDKit be a bit more quiet
    @return: None
    """
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)


def molname_clean(molname: str) -> str:
    """
    Takes a string and removes anything inside parenthesis. This function is used
    because there are duplicates with the same molecule name but with redundant
    information inside parenthesis. Some names are also defined on a dictionary,
    to retrieve them from pubchempy.

    Args:
        molname: name of the molecule.

    Returns:
        cleaned name of the molecule (without parenthesis).
    """
    pattern = re.compile(r" \(([^\)]+)\)| \[([^\)]+)\]")
    # Some molecule names that I need to fix for pubchempy on get_smiles_from_name()
    mol_dictionary = {
        "12-METHOXY-4,4-BISNOR-5beta-8,11,13-PODOCARPATRIEN-3-OL": "12-METHOXY-4,4-BISNOR-5alpha-8,11,13-PODOCARPATRIEN-3-OL",
        "DESLORATADINE HYDROCHLORIDE": "DESLORATADINE",
        "DMSO+FSK": "DMSO+FSK",
        "EPITHEAFLAVIC ACID": "Epitheaflavic acid 3'-gallate",
        "GENTAMICIN SULFATE": "GENTAMICIN",
        "PICEID METHYL ETHER": "PICEID",
        "Rapamycin-FSK": "Rapamycin",
        "Roscovitine-FSK": "Roscovitine",
    }
    new_string = pattern.sub("", molname)
    if molname in mol_dictionary.keys():
        new_string = mol_dictionary[molname]
    return new_string


def get_molecules_png(
    smiles: list,
    labels: list,
    n_cols: int = None,
    img_size: tuple = (300, 300),
    match_substructs: bool = False,
    font_path: str = None,
    font_size: int = 15,
    n_jobs: int = 1,
) -> PIL.Image.Image:
    """
    Args:
        smiles: list of smiles to display.
        labels: List of labels to be written on the mol images.
        n_cols: N columns with molecules. If None, n_cols = len(smiles).
            Defaults to None.
        match_substructs: Will compute the 2D coodinates so molecular structures
            match each other. Defaults to False.
        font_path: Path to a font to be used by the function. Defaults to None.
        font_size: Font size. Defaults to 15.
        n_jobs: Generate mol images in parallel. Defaults to 1.

    Returns:
        PIL.Image.Image
    """

    def list_split(list_a, chunk_size):
        """
        Returns a generator with a determined chunk size over list_a.
        Used to organize the figures in the correct `n_cols` format.
        """
        for i in range(0, len(list_a), chunk_size):
            yield list_a[i : i + chunk_size]

    images = list()
    if match_substructs:
        images = smilist_to_matching_img(smiles, size=img_size)
    else:
        with Pool(n_jobs) as pool:
            images = pool.map(partial(smi_to_img, size=img_size), smiles)

    if labels is not None:
        # Setting the configuration of the font
        if font_path is None:
            font_path = Path(
                "/mnt/c/Users/david/AppData/Local/Microsoft"
                "/Windows/Fonts/JetBrains Mono Light Italic Nerd Font Complete.ttf"
            )
        assert font_path.exists(), "Font path does not exist. "

        font = ImageFont.truetype(str(font_path), font_size)
        img_width, img_height = images[0].size

        # Writing the labels to each of the images
        for img, text in zip(images, labels):
            # Getting the size of the text to center the loation of the text
            font_width, font_height = font.getsize(text)

            # TODO: improve this so we can instead separate the text in two lines
            if font_width > img_width:
                while True:
                    font_size -= 1
                    font = ImageFont.truetype(str(font_path), font_size)
                    font_width, font_height = font.getsize(text)
                    if font_width < img_width:
                        break

            centered_w = (img_width - font_width) / 2
            centered_h = (img_height - font_height) / 99
            draw = ImageDraw.Draw(img)
            draw.text((centered_w, centered_h), text, fill="black", font=font)

    # Splitting the list of images into a list of lists with n_cols
    if n_cols is None:
        n_cols = len(smiles)
    images = list(list_split(images, n_cols))
    # Appending blank figures so we have the correct vector shapes
    while len(images[-1]) < len(images[0]):
        images[-1].append(Image.new("RGB", img_size, color=(255, 255, 255)))

    list_of_hstacked = list()
    # Creating list of horizontally stacked arrays
    for sublist in images:
        list_of_hstacked.append(np.hstack([np.asarray(img) for img in sublist]))
    # Vertically stacking horizontal arrays
    for item in list_of_hstacked:
        final_img = np.vstack([hstack for hstack in list_of_hstacked])
    # Creating and returning image from array
    final_img = Image.fromarray(final_img)
    return final_img


def neutralize_atoms(mol):
    """
    Author: Noel O`Boyle (Vincent Scalfani adapted code for RDKit)
    """
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol


def sanitize_smiles(smi: str) -> str:
    """
    Removes sulfurs, extermal molecules and salts, neutralizes charges
    using the function `neutralize_atoms` and returns the resulting smiles.

    Args:
        smi: single SMILES string to be sanitized.

    Returns:
        sanitized SMILES string.
    """

    salts = re.compile(r"\..?Cl|\..?Br|\..?Ca|\..?K|\..?Na|\..?Li|\..?Zn|/\..?Gd")
    s_acid_remover = re.compile(r"\.OS\(\=O\)\(\=O\)O")
    boron_pattern = re.compile(r"B")
    remover = SaltRemover(defnData="[Cl,Br,Ca,K,Na,Zn]")
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")

    mol = Chem.MolFromSmiles(smi)

    # Removing sulfuric acid (smiles = .OS(=O)(=O)O)
    if s_acid_remover.findall(smi):
        smi = re.sub(s_acid_remover, "", smi)
        try:
            Chem.MolFromSmiles(smi)
        except:
            print(f"{smi} could not be parsed after removing sulfuric acids!")
            return None

    # Removing external molecules by splitting on . and picking the largest smiles
    if "." in smi:
        smi = max(smi.split("."), key=len)
        try:
            mol = Chem.MolFromSmiles(smi)
        except:
            print(f"Compound, ({smi}) could not be parsed!!")
            return None

    # Trying to remove the salts
    if salts.findall(smi):
        res, deleted = remover.StripMolWithDeleted(mol)
        # avoid neutralizing smiles with boron atoms
        if all([res is not None, not boron_pattern.findall(smi)]):
            neutralize_atoms(res)
            # If it didn't remove, let's continue
            if salts.findall(Chem.MolToSmiles(res)):
                print(f"Unable to remove salts from compound {smi}")
                return None
            else:
                smi = Chem.MolToSmiles(res)
                mol = Chem.MolFromSmiles(smi)

    # Are the molecules charged according to the "pattern" variable?
    if mol.GetSubstructMatches(pattern):
        res, deleted = remover.StripMolWithDeleted(mol)
        # avoid neutralizing smiles with boron atoms
        if all([res is not None, not boron_pattern.findall(smi)]):
            neutralize_atoms(res)
        if salts.findall(Chem.MolToSmiles(res)):
            print(f"Unable to remove salts from compound {smi} after neutralizing")
            return None
        else:
            smi = Chem.MolToSmiles(res)

    return smi


def get_smiles_from_name(
    comp_name: str, smi_type: str, verbose: bool = False
) -> str or bool:
    """
    Gets the list of compounds fetched with PubChemPy for the compound
    name provided. Will return false if no compounds are retrieved. Smiles
    from the resulting list then go through sanitization, depending on the
    `smi_type` argument.

    Args:
        comp_name: name of the molecule to be searched.
        smi_type: type of smiles to return. Either `non_isomeric` or `parent`.
        verbose: turn on/off print statements. Defaults to False.

    Raises:
        AttributeError: if `smi_type` not in [`non_isomeric`, `parent`].

    Returns:
        `False` if no compounds are found, otherwise the obtained smiles.
    """

    if smi_type not in ["non_isomeric", "parent"]:
        raise AttributeError(
            "Invalid smi_type attribute. Check documentation \
            for the available types."
        )

    comp_smiles = []
    if comp_name == "DMSO+FSK":
        return False
    try:
        comp_list = pcp.get_compounds(comp_name, "name")
        comp_list[0]
    except:
        if verbose:
            print("No compounds in PubChem with this name", comp_name)
        return False

    for c in comp_list:
        if smi_type == "non_isomeric":
            comp_smiles.append(c.canonical_smiles)
        elif smi_type == "parent":
            comp_smiles.append(c.isomeric_smiles)

    # Making a set to have only unique smiles
    candidates = set()
    if smi_type == "non_isomeric":
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

    # From the sanitized molecules, let's take the smallest
    final_candidate = min(candidates, key=len)

    # Finally, check again for any charges...
    charge = re.compile(r"\+|\-!\d")
    if smi_type == "non_isomeric":
        if charge.findall(final_candidate):
            final_mol = Chem.MolFromSmiles(final_candidate)
            neutralize_atoms(final_mol)
            final_candidate = Chem.MolToSmiles(final_mol)
    return final_candidate


def chembl_smi_standardizer(smi: str) -> tuple:
    """
    Returns a tuple containing `parent smiles` and a `bool`,
    which defaults to if the standardizer failed and false when
    there were no errors.

    Args:
        smi: smiles string to be standardized with the `chembl_structure_pipeline`.

    Returns:
        Tuple containing the parent smiles and a bool indicating if the standardizer
        failed. If True -> standardizer failed..
    """
    mol = Chem.MolFromSmiles(smi)

    standard_mol = standardizer.standardize_mol(mol)
    result = standardizer.get_parent_mol(
        standard_mol
    )  # Tuple with molecule in #0 and Boolean in #1
    # Boolean states whether there was an exclusion flag. For more details, check:
    # https://github.com/chembl/ChEMBL_Structure_Pipeline/wiki/Exclusion-Flag

    parent_mol = result[0]
    parent_smi = Chem.MolToSmiles(
        parent_mol, kekuleSmiles=False, canonical=True, isomericSmiles=True
    )

    if result[1]:
        return parent_smi, True
    else:
        return parent_smi, False


def chembl_mol_standardizer(mol):
    """
    Function that takes RDKit mol objects as input
    and reutrns its parent smiles and a bool, which
    defaults to true in case the standardized failed
    and false when there were no errors.

    Args:
        mol: RDKit mol object.

    Returns:
        Stardardized parent smiles.
    """
    standard_mol = standardizer.standardize_mol(mol)
    result = standardizer.get_parent_mol(
        standard_mol
    )  # Tuple with molecule in #0 and Boolean in #1
    # Boolean states whether there was an exclusion flag. For more details, check:
    # https://github.com/chembl/ChEMBL_Structure_Pipeline/wiki/Exclusion-Flag

    parent_mol = result[0]
    parent_smi = Chem.MolToSmiles(parent_mol, kekuleSmiles=False, canonical=True)

    if result[1]:
        return parent_smi, True
    else:
        return parent_smi, False


def encode_check(filename, verbose=False) -> str:
    """
    Function source: https://stackoverflow.com/questions/436220/how-to-determine-the-encoding-of-text
    Checking the encoding of a file using Python. For different encodings, add it to the list...

    Args:
        filename: name of the file to check for the encoding

    Returns:
        Encoding type of the file (str).
    """
    encodings = ["utf-8", "windows-1250", "windows-1252"]
    for e in encodings:
        try:
            fh = codecs.open(filename, "r", encoding=e)
            fh.readlines()
            fh.seek(0)
        except UnicodeDecodeError:
            if verbose:
                print(f"got unicode error with {e}, trying different encoding")
            else:
                pass
        else:
            if verbose:
                print(f"opening the file with encoding:  {e}")
            else:
                pass
            return e


def get_morganFP(mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=4096, useChirality=False)


def smi_to_img(smi: str, size: tuple = (300, 300)) -> Image:
    """
    Wrapper function to convert smiles to images using RDKit.

    Args:
        smi: smiles string
        size: size of the image. Default: (300, 300)

    Returns:
        RDKit image object.
    """
    mol = Chem.MolFromSmiles(smi)
    img = Draw.MolToImage(mol, size=size)
    return img


def smilist_to_matching_img(smi_list: str, size: tuple = (300, 300)) -> Image:
    """
    Wrapper function to convert smiles to images using RDKit.

    Args:
        smi: smiles string
        size: size of the image. Default: (300, 300)

    Returns:
        RDKit image object.
    """
    mols = [Chem.MolFromSmiles(smi) for smi in smi_list]
    AllChem.Compute2DCoords(mols[0])
    for idx in range(1, len(smi_list)):
        AllChem.GenerateDepictionMatching2DStructure(mols[idx], mols[0])
    imgs = [Draw.MolToImage(mol, size=size) for mol in mols]
    return imgs


def smi_to_fp(smi: str, fp_name):
    """Function to convert smiles to fingerprints using RDKit.

    Args:
        smi: smiles string
        fp_name: name of the desired fingerprint. Available options:
        [`rdkpattern`, `atompair`, `ecfp4`]

    Returns:
        RDKit fingerprint object.
    """
    fp_dictionary = {
        "rdkpattern": Chem.rdmolops.RDKFingerprint,
        "atompair": rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect,
        "ecfp4": partial(
            AllChem.GetMorganFingerprintAsBitVect,
            radius=2,
            nBits=8192,
            useChirality=False,
        ),
    }
    fingerprint_function = fp_dictionary[fp_name]
    mol = Chem.MolFromSmiles(smi)
    smi_fp = fingerprint_function(mol)
    return smi_fp


def same_from_fps(fp1, fp2):
    """
    Calculates and returns the tanimoto similarity between the
    fingerprints (inputs) of two different molecules.

    Args:
        fp1: Fingerprint of molecule 1.
        fp2: Fingerprint of molecule 2.

    Returns:
        tanimoto similarity: float.
    """
    tani = DataStructs.FingerprintSimilarity(fp1, fp2)
    if tani > 0.99:
        return True
    else:
        return False


def smi_to_connectivity(smi):
    mol = Chem.MolFromSmiles(smi)
    connectivity = Chem.MolToInchiKey(mol).split("-")[0]
    return connectivity


def same_from_substructure(smi_pair: tuple) -> bool:
    mol1 = Chem.MolFromSmiles(smi_pair[0])
    mol2 = Chem.MolFromSmiles(smi_pair[1])
    if mol1.HasSubstructMatch(mol2) and mol2.HasSubstructMatch(mol1):
        return True
    else:
        return False


def check_repeated_mols(smilist: list, fp_func=get_morganFP, njobs: int = 1):
    """
    Takes a list of smiles and does a pairwise comparison across
    all components to search for identic molecules. Function can be run in
    parallel using the `njobs` argument.

    Args:
        smilist: list of smiles for which to check for duplicates.
        method: Fingerprint that will be used to check for duplicates
        using `same_from_fps`. Defaults to get_morganFP.
        njobs: Number of jobs to run in parallel. Defaults to 1.

    Returns:
        a nested list, corresponding to the index of the repeated molecules.
    """

    def flatten(x):
        result = []
        for el in x:
            if hasattr(el, "__iter__") and not isinstance(el, str):
                result.extend(flatten(el))
            else:
                result.append(el)
        return result

    molecules = [Chem.MolFromSmiles(smi) for smi in smilist]
    with Pool(njobs) as pool:
        # making this argument customizable to allow for different fingerprinting methods
        fps = pool.map(fp_func, molecules)
        idxs = list(range(len(molecules)))
        fps_tpl = list(combinations(fps, 2))
        idxs_tpl = list(combinations(idxs, 2))
        results = pool.starmap(same_from_fps, fps_tpl)
        duplicate_idxs = list(compress(idxs_tpl, results))

    allrepeats = list()

    # This part could be refactored but for now it works
    for number in np.unique(np.array(duplicate_idxs)):
        repeat = list()
        for x in np.where(np.array(duplicate_idxs) == number)[0]:
            repeat.append(duplicate_idxs[x])
        this_set = list(np.unique(np.array(repeat).flatten()))
        if this_set not in allrepeats:
            allrepeats.append(list(this_set))

    return allrepeats


def activity_to_categorical(mean_std: tuple, R_thresh: float, E_thresh: float):
    """
    Where mean_std corresponds to a tuple with mean at position [0] and
    standard deviation at position [1]. Function used by DataPrepare.define_hits()
    """
    if mean_std[0] + mean_std[1] < R_thresh:
        return "R"
    elif mean_std[0] - mean_std[1] > E_thresh:
        return "E"
    else:
        return "N"


class DataPrepare:
    """
    Callable class for preparing my bioactivity data.
    """

    def __init__(self):
        self.path = Path(__file__).parents[0]  # Path to the current file
        self.dataframe = None
        self.chemstructs_path = self.path / "chem_structurs"
        self.raw_selchemstruct_path = self.chemstructs_path / "sel_chem_structures.smi"
        self.raw_spectrumstruc_path = self.chemstructs_path / "SP-2400_utf8.sdf"
        self.screening_path = self.path / "data/1_bioactivity_data/screening_data"
        self.selchem_path = self.screening_path.glob("z-score*selleck*.csv")
        self.spectrum_path = self.screening_path.glob("z-score*selleck*.csv")
        self.dropped_mols = dict()
        self.failed_getparent = dict()
        # tcontrol stands for treatment controls
        self.selchem_tcontrol = None
        self.spectrum_tcontrol = None

    def prepare_selchem_chemical_structures(self):
        """
        Prepares the chemical structures the Screening.
        1) Reads a raw file containing the SMILES and aditional information
        2) Additional information is dropped from the file, except for a flagging on
        `antimitotic` or `antineoplastic` activity.
        3) Saves the StandardizedSmiles, the ParentSmiles and their isomeric versions,
        as well as the name used for each of the molecules.

        Returns:
            pd.DataFrame: Containing the StandardizedSmiles, the ParentSmiles (...)
        """
        e = encode_check(filename=self.raw_selchemstruct_path)

        sell_chem_df = pd.read_csv(
            self.raw_selchemstruct_path,
            sep="\t",
            lineterminator="\r",
            encoding=e,
        )

        smi_arr = sell_chem_df["Smiles"].values
        parent_smiles = list()
        idx_problem_mols = list()

        for idx, smi in enumerate(smi_arr):
            parent_smi, error = chembl_smi_standardizer(smi)
            if error:
                idx_problem_mols.append(idx)
            parent_smiles.append(parent_smi)

        stand_smiles = [sanitize_smiles(smi) for smi in parent_smiles]
        sell_chem_df["iso_ParentSmiles"] = parent_smiles
        sell_chem_df["iso_StandardizedSmiles"] = stand_smiles

        # Adding antineoplastic activity flag
        key_words = re.compile(r"(antimitotic|antineoplastic)")
        antineo_idx = list()

        for idx, c in enumerate(sell_chem_df["Brief Description"]):
            try:
                if key_words.findall(c):
                    antineo_idx.append(idx)
            except TypeError:
                continue

        flags = [True if idx in antineo_idx else False for idx in sell_chem_df.index]

        # Renaming and dropping some of the columns
        sell_chem_df = sell_chem_df.rename(columns={"Item Name": "MoleculeName"},).drop(
            columns=[
                "Catalog Number",
                "Concentration",
                "Plate Location",
                "Smiles",
                "Target",
                "Brief Description",
                "CAS Number",
            ],
        )

        # Making different columns for isomeric smiles and for smiles without isomeric information:
        parentsmiles = []
        for smi in sell_chem_df["iso_ParentSmiles"]:
            mol = Chem.MolFromSmiles(smi)
            notiso_smi = Chem.MolToSmiles(mol, isomericSmiles=False)
            parentsmiles.append(notiso_smi)

        standardized_smiles = []
        for smi in sell_chem_df["iso_StandardizedSmiles"]:
            mol = Chem.MolFromSmiles(smi)
            notiso_smi = Chem.MolToSmiles(mol, isomericSmiles=False)
            standardized_smiles.append(notiso_smi)

        updated_mol_names = [
            molname_clean(molname) for molname in sell_chem_df["MoleculeName"].tolist()
        ]

        sell_chem_df.assign(
            Screening=["SelleckChem"] * len(sell_chem_df),
            ParentSmiles=parentsmiles,
            StandardizedSmiles=standardized_smiles,
            MoleculeName=updated_mol_names,
            Antineoplastic=flags,
        )
        return sell_chem_df

    def prepare_spectrum_chemical_structures(self):
        # Reading the utf-8 file created
        spectrum_df = PandasTools.LoadSDF(str(self.raw_spectrumstruc_path))

        # Checking for molecules that weren't parsed:
        comp_names = list()
        with open(self.raw_spectrumstruc_path, "r") as source_file:
            for line in source_file:
                if line.startswith(">  <MOLENAME>"):
                    comp_names.append(next(source_file).rstrip("\n"))

        indexes = np.arange(2400)
        for idx in indexes:
            try:
                spectrum_df.loc[idx]
            except:
                # Adding to the dictionary self.dropped_mols
                self.dropped_mols[comp_names[idx]] = "Struct from DataFrame not parsed"

        spectrum_df.reset_index(inplace=True)

        # Adding antineoplastic activity flag
        key_words = re.compile(r"(antimitotic|antineoplastic)")
        antineo_idx = list()

        for idx, c in enumerate(spectrum_df["therap"]):
            try:
                if key_words.findall(c):
                    antineo_idx.append(idx)
            except TypeError:
                continue

        flags = [True if idx in antineo_idx else False for idx in spectrum_df.index]

        # Editing the dataframe:
        spectrum_df.drop(
            ["ID", "index", "ref", "cas#", "source", "therap", "tradename"],
            axis=1,
            inplace=True,
        )

        # Standardizing the compounds from the SPECTRUM library & saving failed:
        mol_array = np.array(spectrum_df["ROMol"])

        parent_smiles = list()
        names = spectrum_df["MOLENAME"]

        for idx, mol in enumerate(mol_array):
            parent_smi, error = chembl_mol_standardizer(mol)
            if error:
                self.failed_getparent[names[idx]] = parent_smi
                # I need to implement a function to visualize these compounds
            parent_smiles.append(parent_smi)

        # Formatting/dropping dataframe columns
        spectrum_df["iso_ParentSmiles"] = parent_smiles
        stand_smi = [sanitize_smiles(smi) for smi in parent_smiles]
        spectrum_df["iso_StandardizedSmiles"] = stand_smi

        # Commenting the collumns that I decided to drop
        spectrum_df.rename(
            columns={
                "MOLENAME": "MoleculeName",
                # "ref": "Reference",
                # "cas#": "CAS#",
                # "source": "Source",
                # "therap": "Bioactivity",
                # "tradename": "TradeName",
            },
            inplace=True,
        )
        spectrum_df.drop(columns=["status", "ROMol"], inplace=True)
        parentsmiles = []
        for smi in spectrum_df["iso_ParentSmiles"]:
            mol = Chem.MolFromSmiles(smi)
            notiso_smi = Chem.MolToSmiles(mol, isomericSmiles=False)
            parentsmiles.append(notiso_smi)

        standardized_smiles = []
        for smi in spectrum_df["iso_StandardizedSmiles"]:
            mol = Chem.MolFromSmiles(smi)
            notiso_smi = Chem.MolToSmiles(mol, isomericSmiles=False)
            standardized_smiles.append(notiso_smi)

        screening = ["Spectrum"] * len(spectrum_df)
        spectrum_df["Screening"] = screening
        spectrum_df["ParentSmiles"] = parentsmiles
        spectrum_df["StandardizedSmiles"] = standardized_smiles
        # Standardizing the molecule names!
        spectrum_df["MoleculeName"] = [
            molname_clean(molname) for molname in spectrum_df["MoleculeName"].tolist()
        ]
        spectrum_df["Antineoplastic"] = flags

        return spectrum_df

    def df_get_chemical_structures(
        self, bioactivity_df: pd.DataFrame, structures_df: pd.DataFrame
    ):
        """
        This function standardizes the molecule names of the two datasets and
        fetches the SMILES from the structures dataframe. Note that no object will
        be returned since the standardized smiles will be appended to the bioactivity_df.
        The molecule names that are not present in this dataset will then be used for
        fetching smiles with pubchempy. This pipeline can be easily modified for
        fetching isomeric smiles but this is not currently supported.

        REMINDER: avoid calling this function when there are repeated compound names.
        Pubchempy will take quite a while to fetch all the SMILES...

        Args:
            bioactivity_df: concatenated dataframe from "self.dataframes_for_modelling()"
            structures_df: output from "prepare_*_chemical_structures()"

        Raises:
            NameError: Whene the dataframe doesn't have the column "MoleculeName".
        """
        # Standard should be compound names under ['MoleculeName']
        if (
            "MoleculeName" not in bioactivity_df.columns
            or "MoleculeName" not in structures_df
        ):
            raise NameError(
                "column 'MoleculeName' not present in one of the dataframes"
            )

        antineo_flags = list()
        stand_smiles = list()
        parent_smiles = list()

        print("Fetching mol structs...")
        for molname in tqdm(bioactivity_df["MoleculeName"]):
            try:
                stand_smi = (
                    structures_df[structures_df["MoleculeName"] == molname]
                    .get("StandardizedSmiles")
                    .values[0]
                )
                parent_smi = (
                    structures_df[structures_df["MoleculeName"] == molname]
                    .get("iso_ParentSmiles")
                    .values[0]
                )
                flag = (
                    structures_df[structures_df["MoleculeName"] == molname]
                    .get("Antineoplastic")
                    .values[0]
                )
            except:
                stand_smi = get_smiles_from_name(molname, smi_type="non_isomeric")
                parent_smi = get_smiles_from_name(molname, smi_type="parent")
                # get_smiles_from_name() will return False in case it fails..
                if stand_smi:
                    stand_smiles.append(stand_smi)
                    antineo_flags.append(False)  # Don't know since came from pubchempy
                else:
                    stand_smiles.append(np.NaN)
                    antineo_flags.append(False)
                    self.dropped_mols[molname] = "Unknown structure"
                # The same for parent smiles...
                if parent_smi:
                    parent_smiles.append(parent_smi)
                else:
                    parent_smiles.append(np.NaN)
                continue
            stand_smiles.append(stand_smi)
            parent_smiles.append(parent_smi)
            antineo_flags.append(flag)
        bioactivity_df["StandardizedSmiles"] = stand_smiles
        bioactivity_df["Antineoplastic"] = antineo_flags
        bioactivity_df["ParentSmiles"] = parent_smiles

        try:
            todrop = bioactivity_df[
                (bioactivity_df["StandardizedSmiles"].isna())
                & (bioactivity_df["MoleculeName"] != "DMSO+FSK")
            ].index
            bioactivity_df.drop(index=todrop, inplace=True).reset_index(
                inplace=True, drop=True
            )
        except AttributeError:
            print("No NaN StandardizedSmiles detected")

    def define_hits(
        self, bioactivity_df: pd.DataFrame, R_threshold: float, E_threshold: float
    ) -> pd.DataFrame:
        """
        This function will return a the bioactivity_df with only Treatment
        molecules (column in pd.DataFrame). This new dataframe will have an
        activity column describing whether a certain molecule is a cyst size
        `enhancer (E)`, a `reducer (R)` or `inactive (N)`. A compound is an
        reducer when its area_Rhodamine_Mean + SD (standard deviation) are
        below the given `threshold` -> currently `DMSO+FSK mean - SD`.
        Compounds that go beyond DMSO+FSK mean + SD are encoded as enhancers.

        The input dataframe for this function should be the one created after
        applying dataframes_for_modelling().

        Args:
            bioactivity_df: pd.DataFrame output from DataPrepare.dataframes_for_modelling().
            R_threshold: bioactivity_df[DMSO+FSK] (MEAN) - STANDARD DEVIATION.
            E_threshold: bioactivity_df[DMSO+FSK] (MEAN) + STANDARD DEVIATION.

        Returns:
            treatments_df(pd.DataFrame) -> with hit molecules.
        """
        treatments_df = bioactivity_df[bioactivity_df["Control"] == "Treatment"].copy(
            deep=True
        )
        treatments_df.reset_index(inplace=True, drop=True)

        mean_SDs = list(
            zip(
                treatments_df["area_Rhodamine_Mean"].values,
                treatments_df["area_Rhodamine_SD"].values,
            )
        )

        # activities = ["A" if i in threshold_idx else "N" for i in treatments_df.index]
        with Pool(5) as pool:

            activities = pool.map(
                partial(
                    activity_to_categorical, R_thresh=R_threshold, E_thresh=E_threshold
                ),
                mean_SDs,
            )

        treatments_df["Activity"] = activities

        print(f"Total amount of compounds: {len(treatments_df)}")
        subset_df = treatments_df[treatments_df["Activity"] == "E"]
        print(f"Number of cyst-swelling enhancing compounds: {len(subset_df)}")

        subset_df = treatments_df[treatments_df["Activity"] == "R"]
        print(f"Number of cyst-swelling reducing compounds: {len(subset_df)}")

        subset_df = treatments_df[treatments_df["Activity"] == "N"]
        print(f"Number of inactive compounds: {len(subset_df)}")
        return treatments_df

    def drop_repeats_from_dataset(
        self, bioactivity_df: pd.DataFrame, smiles_col: pd.DataFrame, njobs: int = 1
    ):
        """
        WARNING: THIS FUNCTION IS BROKEN, I STILL HAVE TO FIX IT

        Takes as input the normalized dataframe from `get_merged_datasets()`,
        the total dataframe with all chemical structures 'prepare_*_chemical_structures()'
        and the output of the function 'get_repeats_dataset(). This will drop the
        selected molecules from the dataframe according to the following rules:

        1) Bioactivity - If both repeats' obj.Mean(area) standard deviation is < 1, then ;

        Args:
            bioactivity_df: The merged bioactivity dataframe.
            chemstructs_df: The merged chemical structures dataframe.
            njobs: Checking for duplicate molecules in parallel. Defaults to 1.

        Returns:
            Two dataframes, where [0] contains the duplicates and [1] contains
            the ones that will be dropped.
        """
        repeats = check_repeated_mols(bioactivity_df["StandardizedSmiles"], njobs=njobs)

        repeats_df = pd.DataFrame(columns=bioactivity_df.columns)
        # Making a new dataframe with a "Repeats" index for the duplicates
        for idx, repeat in enumerate(repeats):
            subset_df = bioactivity_df.loc[repeat]
            subset_df["Repeats"] = [idx] * len(subset_df)

            ihave = list()
            names = list(subset_df["MoleculeName"])
            for molname in names:
                # Checking if the molecule is not from pubchem
                if molname in list(chemstructs_df["MoleculeName"]):
                    ihave.append(True)
                else:
                    ihave.append(False)

            subset_df["Ihave"] = ihave
            repeats_df = pd.concat([repeats_df, subset_df], ignore_index=True)

        # picked is a list of selected molecules based on rules 1, 2, 3.
        picked = []
        for r in repeats_df["Repeats"].unique():
            subset_df = repeats_df[repeats_df["Repeats"] == r]
            activeonly = subset_df[
                subset_df["Activity"] != "N"
            ]  ## modified (do I get errors?)
            haveonly = subset_df[subset_df["Ihave"] == True]
            if len(activeonly) == 1:
                picked.append(activeonly.iloc[0]["MoleculeName"])
                continue
            if len(haveonly) == 1:
                picked.append(haveonly.iloc[0]["MoleculeName"])
                continue
            smallerSD_idx = subset_df["area_Rhodamine_SD"].idxmin()
            picked.append(subset_df.loc[smallerSD_idx]["MoleculeName"])

        # NOT IN
        todrop_df = repeats_df[~repeats_df["MoleculeName"].isin(picked)]
        # Adding information on the self.dropped_mols
        for molname in todrop_df["MoleculeName"]:
            self.dropped_mols[molname] = "Repeated chemical structure"

        print(f"{len(todrop_df)} repeated molecules will be removed from the dataset")
        return repeats_df, todrop_df

    def get_antineoplastic_from_chembl(self):
        """
        Returns a dictionary with the names of antineoplastic compounds with
        the keys [iso_smiles] and [noniso_smiles].
        """

        def get_noniso_smiles(smi):
            mol = Chem.MolFromSmiles(smi)
            noniso_smi = Chem.MolToSmiles(
                mol, kekuleSmiles=False, canonical=True, isomericSmiles=False
            )
            return noniso_smi

        molecule = new_client.molecule
        approved_drugs = molecule.filter(max_phase=4).order_by(
            "molecule_properties__mw_freebase"
        )
        approved_drugs.set_format("json")

        # L01 referes to the atc antineoplastic class, such as described in:
        # https://www.whocc.no/atc_ddd_index/?code=L01&showdescription=yes
        antineoplastic = re.compile("(L01)")
        drug_dict = dict()

        for drug in approved_drugs:

            name = drug["pref_name"]

            try:
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
                atc_class = " ".join(drug_dict[name]["atc_class"])
                if antineoplastic.findall(atc_class):
                    drug_dict[name]["antineoplastic"] = True

        antineo_drugs = dict()

        for key in drug_dict:
            if drug_dict[key]["antineoplastic"]:
                antineo_drugs[key] = dict()
                antineo_drugs[key]["chembl_id"] = drug_dict[key]["chembl_id"]
                antineo_drugs[key]["iso_smiles"] = drug_dict[key]["smiles"]
                antineo_drugs[key]["noniso_smiles"] = get_noniso_smiles(
                    drug_dict[key]["smiles"]
                )
        return antineo_drugs

    def flag_antineoplastic_compounds(
        self,
        screening_df: pd.DataFrame,
        antineo_drugs: dict,
        simi_threshold: float,
        fingerprint="ecfp4",
    ):
        """
        Adds screening_df['Antineoplastic'] with True based on the tanimoto similarity.

        Params:
        screening_df -> merged bioactivity dataset with chemical structures.
        antineo_drugs -> output from get_antineoplastic_from_chembl().
        simi_threshold -> 0 < float value < 1.
        fingerprint -> fingerprint name to be used on the smi_to_fp() function.
        """

        known_flags = screening_df["Antineoplastic"].values
        # Antineoplastic compounds from ChEMBL
        antineo_smiles = [
            antineo_drugs[key]["noniso_smiles"] for key in antineo_drugs.keys()
        ]
        # Antineoplastic compounds from dataset
        for idx, flag in enumerate(known_flags):
            if flag:
                antineo_smiles.append(screening_df.iloc[idx]["StandardizedSmiles"])

        with Pool(5) as pool:
            antineo_fp_arr = pool.map(
                partial(smi_to_fp, fp_name=fingerprint), antineo_smiles
            )
            data_fp_arr = pool.map(
                partial(smi_to_fp, fp_name=fingerprint),
                screening_df["StandardizedSmiles"],
            )

        antineo_flag = list()

        for antineo_fp in tqdm(antineo_fp_arr):
            for idx, fp in enumerate(data_fp_arr):
                # When we find a True bioactivity flag, we don't need to assess again
                if known_flags[idx]:
                    antineo_flag.append(idx)
                    continue
                tani = DataStructs.FingerprintSimilarity(antineo_fp, fp)
                if tani >= simi_threshold:
                    antineo_flag.append(idx)

        antineo_flag = list(set(antineo_flag))
        print(f"{len(antineo_flag)} antineoplastic compounds identified in the dataset")
        antineo_flag = [
            True if idx in antineo_flag else False for idx in screening_df.index
        ]
        screening_df["Antineoplastic"] = antineo_flag

    @staticmethod
    def get_connectivity_dict(
        screening_df: pd.DataFrame,
        filtered_data: pd.DataFrame,
        method: str,
        fingerprint: str = "ecfp4",
    ) -> dict:
        """
        Params:
            screening_df -> Concatenated bioactivity dataframes.
            filtered_data -> Output from papyrus_data_prepare.py.
            method -> How to find the same molecule: either `connectivity`,
                      `fp_similarity` or `same_substructure`.
            Note: `same_substructure` will lead to faulty results in some cases.
        """

        if method not in ["connectivity", "fp_similarity", "same_substructure"]:
            raise AttributeError(
                "method should be either `connectivity`, `fp_similarity`, or `same_substructure`"
            )

        smi_to_connect_dict = dict()

        # Dictionary with the connectivity for unique SMILES
        for smi in filtered_data["SMILES"].unique():
            subset_df = filtered_data[filtered_data["SMILES"] == smi]
            connect = subset_df["connectivity"].unique()
            # Verify if there's more than 1 connectivity per SMILES
            if len(connect) > 1:
                print("watch out for ", smi)
            else:
                smi_to_connect_dict[smi] = connect[0]

        data_smis = list(smi_to_connect_dict.keys())
        names = screening_df["MoleculeName"].values

        if method == "connectivity":
            my_smis = screening_df["ParentSmiles"]
        else:
            my_smis = screening_df["StandardizedSmiles"].values

        connectivity_dict = dict()

        if method == "fp_similarity":
            with Pool(5) as pool:
                my_fps = pool.map(partial(smi_to_fp, fp_name=fingerprint), my_smis)
                toquery_fps = pool.map(
                    partial(smi_to_fp, fp_name=fingerprint), data_smis
                )
                # Compare fingerprints from my dataset to Papyrus
                for idx, fp in enumerate(tqdm(my_fps)):
                    for i, f in enumerate(toquery_fps):
                        if same_from_fps(fp, f):
                            connectivity_dict.update(
                                {names[idx]: smi_to_connect_dict[data_smis[i]]}
                            )  # The name of the compound : unique connectivity

        elif method == "connectivity":
            # Reverse keys and values from smi_to_connect_dict
            connect_to_smi_dict = {v: k for k, v in smi_to_connect_dict.items()}
            with Pool(5) as pool:
                my_connects = pool.map(smi_to_connectivity, my_smis)
                query_connects = subset_df["connectivity"].unique()
                # Compare connectivities from my dataset to Papyrus
                for idx, my_c in enumerate(tqdm(my_connects)):
                    for i, c in enumerate(query_connects):
                        if my_c == c:
                            connectivity_dict.update(
                                {names[idx]: connect_to_smi_dict[query_connects[i]]}
                            )  # The name of the compound : unique connectivity

        elif method == "same_substructure":
            with Pool(8) as pool:
                my_smiles_to_name_dict = {
                    smi: name for smi, name in zip(my_smis, names)
                }
                combination_smis = list(product(my_smis, data_smis))
                # Will return a list of bool
                is_same = tqdm(
                    pool.imap(
                        partial(same_from_substructure),
                        combination_smis,
                    ),
                    total=len(combination_smis),
                )
                # Only the smiles that represent the same molecule
                same_molecules = compress(combination_smis, is_same)
                for smi_pair in same_molecules:
                    connectivity_dict.update(
                        {
                            my_smiles_to_name_dict[smi_pair[0]]: smi_to_connect_dict[
                                smi_pair[1]
                            ]
                        }
                    )

        screening_df["Connectivity"] = [
            connectivity_dict[name] if name in connectivity_dict.keys() else np.NaN
            for name in screening_df["MoleculeName"]
        ]
        return connectivity_dict

    @staticmethod
    def get_papyrus_protein_info(input_df: pd.DataFrame):
        """
        Params:
        input_df -> papyrus bioactivity dataframe

        output: the dataframe with the proteins
        """
        papyrus_path = (
            Path(__file__).home()
            / ".data/papyrus/05.5/05.5_combined_set_protein_targets.tsv.xz"
        )
        papyrus_proteins = pd.read_csv(
            papyrus_path,
            sep="\t",
            keep_default_na=False,
        )
        result_subset = papyrus_proteins[
            papyrus_proteins["target_id"].isin(input_df["target_id"])
        ]
        return result_subset

    @staticmethod
    def update_papyrus_protclass(papyrus_protein_df: pd.DataFrame):
        """
        Params:
        papyrus_protein_df -> The papyrus protein dataset
        """

        def get_class0(string):
            return string.split("->")[0]

        def get_class1(string):
            cla_list = string.split("->")
            try:
                string = cla_list[1]
            except IndexError:
                string = "Unknown"
            return string

        papyrus_protein_df["Class[0]"] = papyrus_protein_df["Classification"].apply(
            partial(get_class0)
        )
        papyrus_protein_df["Class[1]"] = papyrus_protein_df["Classification"].apply(
            partial(get_class1)
        )

        return papyrus_protein_df

    @staticmethod
    def get_data_from_papyrus(input_df: pd.DataFrame, papyrus_df: pd.DataFrame):
        """
        Params:
        input_df -> subset of the output from main.py with columns `MoleculeName`, `Connectivity`.
        papyrus_df -> dataframe obtained from the papyrus_data_prepare.py module.

        returns the dataset with extra information such as [Activity, Class[0] and Class[1]].
        """

        connect_dict = dict()
        query_connectivities = input_df[~input_df["Connectivity"].isnull()][
            "Connectivity"
        ].unique()
        # Removing nan values (yield c==c as False):
        # query_connectivities = [c for c in query_connectivities if c == c]

        # Dictionary connecting molecule names and connectivities
        for connect in query_connectivities:
            try:
                connect_dict[connect] = "~".join(
                    input_df[input_df["Connectivity"] == connect]["Compound"].values
                )
            except:
                print(f"Exception reached for {connect}")
                continue

        mymols_papyrus_df = papyrus_df[
            papyrus_df["connectivity"].isin(query_connectivities)
        ].copy(deep=True)

        def apply_activity_class(value, threshold):
            if value > threshold:
                return "A"
            else:
                return "N"

        threshold = 6.0
        while True:
            if threshold > 9:
                break
            mymols_papyrus_df[f"Activity_{threshold:.1f}"] = mymols_papyrus_df[
                "pchembl_value_Median"
            ].apply(partial(apply_activity_class, threshold=threshold))

            threshold += 0.5

        my_proteins_df = DataPrepare.get_papyrus_protein_info(mymols_papyrus_df)
        my_proteins_df = DataPrepare.update_papyrus_protclass(my_proteins_df)

        class_0 = list()
        class_1 = list()
        print("Updating target information")
        for t in tqdm(my_proteins_df["target_id"]):
            description0 = my_proteins_df[my_proteins_df["target_id"] == t]["Class[0]"]
            description1 = my_proteins_df[my_proteins_df["target_id"] == t]["Class[1]"]
            subset_df = mymols_papyrus_df[mymols_papyrus_df["target_id"] == t]
            class_0 += list(description0) * len(subset_df)
            class_1 += list(description1) * len(subset_df)

        # Updating additional information to the final dataframe
        mymols_papyrus_df["Class[0]"] = class_0
        mymols_papyrus_df["Class[1]"] = class_1
        mymols_papyrus_df["Compond"] = [
            connect_dict[c] for c in mymols_papyrus_df["connectivity"]
        ]
        mymols_papyrus_df.reset_index(inplace=True, drop=True)
        return mymols_papyrus_df


class DataVisualize:
    """
    Collection of static methods for plotting the data
    generated from the modules main.py and papyrus_data_prepare.py
    """

    def __init__(self) -> None:
        pass

    @staticmethod
    def compare_controls(normaliz_df: pd.DataFrame, name: str, save_path: str):
        """
        Function for visualizing the data distribution of the positive and
        negative controls and their standard deviation.

        Params:
        normaliz_df = Screening dataframe with normalized values.
                      (Output from DataPrepare.normalized_*_bioactivity())
        name = Name of the screening for the title of the plot.
        save_path = Path for saving the plot. If save == '', it won't be saved.
        """

        subset_df = normaliz_df[normaliz_df["Control"].isin(["neg", "pos"])]

        fig, ax = plt.subplots(figsize=(8, 4))

        neg_control_df = subset_df[subset_df["Control"] == "neg"]
        pos_control_df = subset_df[subset_df["Control"] == "pos"]

        negmean = float(neg_control_df.groupby("MoleculeName").mean()["area_Rhodamine"])
        negSD = float(neg_control_df.groupby("MoleculeName").std()["area_Rhodamine"])
        posmean = float(pos_control_df.groupby("MoleculeName").mean()["area_Rhodamine"])
        posSD = float(pos_control_df.groupby("MoleculeName").std()["area_Rhodamine"])

        neg_sd1 = negmean - negSD
        neg_sd2 = negmean + negSD
        pos_sd1 = posmean - posSD
        pos_sd2 = posmean + posSD

        sns.set_style("whitegrid")
        sns.histplot(
            data=subset_df,
            x="area_Rhodamine",
            bins=15,
            alpha=0.3,
            kde=True,
            hue="Control",
        )

        ####### Plotting the means #######
        ax.vlines(
            x=negmean,
            color="blue",
            linestyle=(0, (5, 10)),
            ymin=0,
            ymax=100,
            label="DMSO_mean",
        )
        ax.vlines(
            x=posmean,
            color="red",
            linestyle=(0, (5, 10)),
            ymin=0,
            ymax=60,
            label="DMSO+FSK_mean",
        )

        # ####### Plotting the SDs #######
        ax.vlines(
            x=[neg_sd1, neg_sd2],
            color="blue",
            linestyle=(0, (1, 10)),
            label="DMSO_mean ± SD",
            ymin=0,
            ymax=80,
        )
        ax.vlines(
            x=[pos_sd1, pos_sd2],
            color="red",
            linestyle=(0, (1, 10)),
            label="DMSO+FSK_mean ± SD",
            ymin=0,
            ymax=40,
        )

        ax.legend(
            title="Legend",
            labels=[
                "DMSO+FSK",
                "DMSO",
                "DMSO mean",
                "DMSO+FSK mean",
                "DMSO mean ± SD",
                "DMSO+FSK mean ± SD",
            ],
            bbox_to_anchor=(1.04, 0),
            loc="lower left",
            borderaxespad=0,
        )
        ax.set_xlabel("Rhodamine Area (Z-score normalized per plate)")
        ax.set_title(f"{name} Screen: positive & negative controls")
        if save_path == "":
            return
        if os.path.isdir(save_path) == False:
            print("Figure not saved. Must enter a path to a directory")
            return
        else:
            fig.savefig(
                os.path.join(save_path, f"{name}_controls.png"),
                format="png",
                dpi=1200,
                bbox_inches="tight",
            )
        return fig

    @staticmethod
    def plot_selchem_hits(
        normaliz_df: pd.DataFrame,
        screening_df: pd.DataFrame,
        antineo: bool,
        public_data: bool,
        allhits: bool = False,
        tophits: bool = False,
    ):
        """
        Params:
        normaliz_df -> selleck_normal_df generated from the function data_prepare.normalized_selchem_bioactivity()
        screening_df -> final dataset generated from the main.py module
        antineo -> whether to plot the antineoplastic compounds or not
        public_data -> whether to plot the compounds with public data available or not
        """
        actives_df = screening_df[screening_df["Activity"] == "A"]
        if allhits:
            activ_names = actives_df["MoleculeName"].unique()
            if tophits:
                fig, ax = plt.subplots(figsize=(6, 4))
            else:
                fig, ax = plt.subplots(figsize=(6, 20))
        else:
            if public_data:
                if antineo:
                    fig, ax = plt.subplots(figsize=(6, 6))
                    activ_names = actives_df[
                        ~(actives_df["Connectivity"].isnull())
                        & (actives_df["Antineoplastic"] == True)
                    ]["MoleculeName"].unique()
                else:
                    ig, ax = plt.subplots(figsize=(6, 18))
                    activ_names = actives_df[
                        ~(actives_df["Connectivity"].isnull())
                        & (actives_df["Antineoplastic"] == False)
                    ]["MoleculeName"].unique()
            else:
                if antineo:
                    fig, ax = plt.subplots(figsize=(6, 3))
                    activ_names = actives_df[
                        (actives_df["Connectivity"].isnull())
                        & (actives_df["Antineoplastic"] == True)
                    ]["MoleculeName"].unique()
                else:
                    fig, ax = plt.subplots(figsize=(6, 6))
                    activ_names = actives_df[
                        (actives_df["Connectivity"].isnull())
                        & (actives_df["Antineoplastic"] == False)
                    ]["MoleculeName"].unique()

        sns.set_style("ticks")
        norma_subset = normaliz_df[
            (normaliz_df["combined.string"] == "DMSO+FSK_pos_2.5")
            | (normaliz_df["MoleculeName"].isin(activ_names))
        ]

        bioactiv_sorted_average = (
            norma_subset.groupby("MoleculeName")
            .median()
            .sort_values("area_Rhodamine", ascending=False)
            .reset_index()
        )
        bioactiv_sorted_average = list(bioactiv_sorted_average["MoleculeName"].values)
        if tophits:
            bioactiv_sorted_average = itemgetter(0, -5, -4, -3, -2, -1)(
                bioactiv_sorted_average
            )

        toplot_df = pd.DataFrame(columns=norma_subset.columns)
        for molname in bioactiv_sorted_average:
            subset_df = norma_subset[norma_subset["MoleculeName"] == molname]
            toplot_df = pd.concat((toplot_df, subset_df), axis=0)

        sns.boxplot(
            data=norma_subset,
            y="MoleculeName",
            x="area_Rhodamine",
            palette="vlag",
            width=0.6,
            whis=[0, 100],
            order=bioactiv_sorted_average,
        )

        sns.stripplot(
            x="area_Rhodamine",
            y="MoleculeName",
            data=norma_subset,
            size=4,
            color=".3",
            linewidth=0,
            order=bioactiv_sorted_average,
        )

        sns.despine(left=True, top=True)
        ax.set_xlim(-5.5, 6)
        ax.xaxis.grid(True)
        ax.set(ylabel="")

        thresh = norma_subset[norma_subset["MoleculeName"] == "DMSO+FSK"]
        thresh = thresh["area_Rhodamine"].mean() - thresh["area_Rhodamine"].std()
        plt.axvline(
            x=thresh,
            color="tab:red",
            linestyle="--",
            label="DMSO+FSK Mean - SD",
            alpha=0.75,
        )
        plt.legend(
            bbox_to_anchor=(0.6, 0.005, 1, 0.2),
            loc="lower left",
            borderaxespad=0,
            ncol=3,
        )
        if tophits:
            plt.title("SelleckChem Screen - positive control & top 5 hit compounds")
            plt.xlabel("Rhodamine Area (Z-score normalized per plate)")
            fig.savefig(
                "Top5_SelChem_HitComps.png", format="png", dpi=1200, bbox_inches="tight"
            )
        if allhits:
            plt.title("SelleckChem Screen\n Positive control & hit compounds")
            plt.xlabel("Rhodamine Area (Z-score normalized per plate)")
            fig.savefig(
                "all_sellchem_hits.png",
                format="png",
                dpi=1200,
                bbox_inches="tight",
            )
        return fig

    @staticmethod
    def plot_spectrum_hits(
        normaliz_df: pd.DataFrame,
        screening_df: pd.DataFrame,
        antineo: bool,
        public_data: bool,
        allhits: bool = False,
        tophits: bool = False,
    ):
        """
        Params:
        normaliz_df -> selleck_normal_df generated from the function data_prepare.normalized_selchem_bioactivity().
        screening_df -> final dataset generated from the main.py module.
        antineo -> True or False whether you want to plot antineoplastic active compounds or not.
        publicdata -> Whether you plot the compounds with public data available or not.
        """
        actives_df = screening_df[screening_df["Activity"] == "A"]

        if allhits:
            activ_names = actives_df["MoleculeName"].unique()
            if tophits:
                fig, ax = plt.subplots(figsize=(6, 4))
            else:
                fig, ax = plt.subplots(figsize=(6, 20))
        else:
            if public_data:
                if antineo:
                    fig, ax = plt.subplots(figsize=(6, 8))
                    activ_names = actives_df[
                        ~(actives_df["Connectivity"].isnull())
                        & (actives_df["Antineoplastic"] == True)
                    ]["MoleculeName"].unique()
                else:
                    fig, ax = plt.subplots(figsize=(6, 9))
                    activ_names = actives_df[
                        ~(actives_df["Connectivity"].isnull())
                        & (actives_df["Antineoplastic"] == False)
                    ]["MoleculeName"].unique()
            else:
                if antineo:
                    fig, ax = plt.subplots(figsize=(6, 8))
                    activ_names = actives_df[(actives_df["Antineoplastic"] == True)][
                        "MoleculeName"
                    ].unique()
                else:
                    fig, ax = plt.subplots(figsize=(6, 12))
                    activ_names = actives_df[(actives_df["Antineoplastic"] == False)][
                        "MoleculeName"
                    ].unique()

        sns.set_style("ticks")
        norma_subset = normaliz_df[
            (normaliz_df["Control"] == "pos")
            | (normaliz_df["MoleculeName"].isin(activ_names))
        ]

        bioactiv_sorted_average = (
            norma_subset.groupby("MoleculeName")
            .median()
            .sort_values("area_Rhodamine", ascending=False)
            .reset_index()
        )
        bioactiv_sorted_average = list(bioactiv_sorted_average["MoleculeName"].values)
        if tophits:
            bioactiv_sorted_average = itemgetter(0, -5, -4, -3, -2, -1)(
                bioactiv_sorted_average
            )

        def capitalize_comps(comp_name):
            if "DMSO" not in comp_name:
                return comp_name.capitalize()
            else:
                return comp_name

        toplot_df = pd.DataFrame(columns=norma_subset.columns)
        for molname in bioactiv_sorted_average:
            subset_df = norma_subset[norma_subset["MoleculeName"] == molname]
            toplot_df = pd.concat((toplot_df, subset_df), axis=0)
        toplot_df["MoleculeName"] = toplot_df["MoleculeName"].apply(capitalize_comps)

        sns.boxplot(
            data=toplot_df,
            y="MoleculeName",
            x="area_Rhodamine",
            palette="vlag",
            width=0.6,
            whis=[0, 100],
        )

        sns.stripplot(
            x="area_Rhodamine",
            y="MoleculeName",
            data=toplot_df,
            size=4,
            color=".3",
            linewidth=0,
        )

        sns.despine(left=True, top=True)
        ax.set_xlim(-5.5, 6)
        ax.xaxis.grid(True)
        ax.set(ylabel="")

        thresh = norma_subset[norma_subset["MoleculeName"] == "DMSO+FSK"]
        thresh = thresh["area_Rhodamine"].mean() - thresh["area_Rhodamine"].std()
        plt.axvline(
            x=thresh,
            color="tab:red",
            linestyle="--",
            label="DMSO+FSK Mean - SD",
            alpha=0.75,
        )
        plt.legend(
            bbox_to_anchor=(0.6, 0.005, 1, 0.2),
            loc="lower left",
            borderaxespad=0,
            ncol=3,
        )
        if tophits:
            plt.title("Spectrum Screen\n Positive control & top 5 hit compounds")
            plt.xlabel("Rhodamine Area (Z-score normalized per plate)")
            fig.savefig(
                "Top5_Spectrum_HitComps.png",
                format="png",
                dpi=1200,
                bbox_inches="tight",
            )
        if allhits:
            plt.title("Spectrum Screen\n Positive control & hit compounds")
            plt.xlabel("Rhodamine Area (Z-score normalized per plate)")
            fig.savefig(
                "all_spectrum_hits.png",
                format="png",
                dpi=1200,
                bbox_inches="tight",
            )
        return fig

    @staticmethod
    def plot_papyrus_activ_percent(
        my_papyrus_df: pd.DataFrame, plt_title: str, save: bool = False
    ):
        """
        Params:
        my_papyrus_df -> output from papyrus_data_prepare.py.
        plt_title -> Title to be displayed on the plot
        """
        median_pchembl = np.array(my_papyrus_df["pchembl_value_Median"])

        fig, ax = plt.subplots(figsize=(8, 5))
        pchembl_plot = sns.histplot(median_pchembl, bins=30, alpha=0.6, kde=True)

        actives_over8 = len(np.where(median_pchembl > 8)[0]) / len(median_pchembl)
        print(f"percentage of targets with bioactivity above 6.0: {actives_over8:.3%}")

        actives_over6 = len(np.where(median_pchembl > 6)[0]) / len(median_pchembl)
        print(f"percentage of targets with bioactivity above 6.0: {actives_over6:.3%}")

        actives_over6_5 = len(np.where(median_pchembl > 6.5)[0]) / len(median_pchembl)
        print(
            f"percentage of targets with bioactivity above 6.5: {actives_over6_5:.3%}"
        )

        plt.axvline(
            x=8,
            color="tab:red",
            linestyle="--",
            label=f"Activity above 8.0: {actives_over8:.3%}",
        )

        # plt.axvline(
        #     x=6,
        #     color="tab:red",
        #     linestyle="--",
        #     label=f"Activity above 6.0: {actives_over6:.3%}",
        # )
        # plt.axvline(
        #     x=6.5,
        #     color="tab:orange",
        #     linestyle="--",
        #     label=f"Activity above 6.5: {actives_over6_5:.3%}",
        # )

        plt.title(f"{plt_title}")
        plt.xlabel("pChEMBL Value (Median)")
        plt.legend()
        if save:
            fig.savefig(
                "Papyrus_data.png",
                format="png",
                dpi=1200,
                bbox_inches="tight",
            )
        return fig

    @staticmethod
    def plot_selchem_treatmentcontrols(meansd_selchem_df, selleck_normal_df):
        """
        Params:
        meansd_selchem_df -> output from DataPrepare.dataframes_for_modelling().
        selleck_normal_df -> output from DataPrepare.normalized_selchem_bioactivity().
        """
        sns.set_palette(sns.color_palette("tab10"))

        pos_control_df = meansd_selchem_df[meansd_selchem_df["Control"] == "pos"]
        poscontrol_threshold = float(pos_control_df["area_Rhodamine_Mean"]) - float(
            pos_control_df["area_Rhodamine_SD"]
        )

        neg_control_df = meansd_selchem_df[meansd_selchem_df["Control"] == "neg"]
        negcontrol_threshold = float(neg_control_df["area_Rhodamine_Mean"]) - float(
            neg_control_df["area_Rhodamine_SD"]
        )
        print(negcontrol_threshold, poscontrol_threshold)

        selchem_controls_df = selleck_normal_df[
            selleck_normal_df["Control"].isin(["ControlTreatment", "pos", "neg"])
        ]

        fig, ax = plt.subplots(figsize=(8, 3))
        plt.axhline(
            y=negcontrol_threshold, color="blue", linestyle="--", label="DMSO mean - SD"
        )
        plt.axhline(
            y=poscontrol_threshold,
            color="red",
            linestyle="--",
            label="DMSO+FSK mean - SD",
        )
        ax = sns.violinplot(
            x="MoleculeName", y="area_Rhodamine", data=selchem_controls_df
        )
        ax.set_ylabel("Area Rhodamine (Z-score)", size=12)
        ax.set_title("SelleckChem screen control treatments", size=14)
        ax.set_xlabel("Compound Name", size=12)
        ax.set_ylim(-6, 6)
        plt.legend(bbox_to_anchor=(1.04, 0), loc="lower left", borderaxespad=0)

        for tick in ax.get_xticklabels():
            tick.set_rotation(45)
            tick.set_fontsize(12)
            tick.set_horizontalalignment("right")

        plt.show()
        return fig

    @staticmethod
    def plot_spectrum_treatmentcontrols(meansd_spectrum_df, spect_normal_df):
        """
        Params:
        meansd_spectrum_df -> output from DataPrepare.dataframes_for_modelling().
        spect_normal_df -> output from DataPrepare.normalized_spectrum_bioactivity().
        """
        pos_control_df = meansd_spectrum_df[meansd_spectrum_df["Control"] == "pos"]
        poscontrol_threshold = float(pos_control_df["area_Rhodamine_Mean"]) - float(
            pos_control_df["area_Rhodamine_SD"]
        )

        neg_control_df = meansd_spectrum_df[meansd_spectrum_df["Control"] == "neg"]
        negcontrol_threshold = float(neg_control_df["area_Rhodamine_Mean"]) - float(
            neg_control_df["area_Rhodamine_SD"]
        )
        print(negcontrol_threshold, poscontrol_threshold)

        specontrol_df = spect_normal_df[
            spect_normal_df["Control"].isin(["ControlTreatment", "pos", "neg"])
        ]

        fig, ax = plt.subplots(figsize=(8, 3))
        plt.axhline(
            y=poscontrol_threshold, color="red", linestyle="--", label="DMSO mean - SD"
        )
        plt.axhline(
            y=negcontrol_threshold,
            color="blue",
            linestyle="--",
            label="DMSO+FSK mean - SD",
        )
        ax = sns.violinplot(x="MoleculeName", y="area_Rhodamine", data=specontrol_df)
        ax.set_ylabel("Area Rhodamine (Z-score)", size=12)
        ax.set_title("Spectrum screen control treatments", size=14)
        ax.set_xlabel("Compound Name", size=12)
        ax.set_ylim(-6, 6)
        plt.legend(bbox_to_anchor=(1.04, 0), loc="lower left", borderaxespad=0)

        for tick in ax.get_xticklabels():
            tick.set_rotation(45)
            tick.set_fontsize(12)
            tick.set_horizontalalignment("right")

        return fig


if __name__ == "__main__":
    pass
