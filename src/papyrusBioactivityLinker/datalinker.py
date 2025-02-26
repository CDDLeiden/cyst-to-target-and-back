import json
import re
from pathlib import Path
from typing import List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from chemFilters.chem.standardizers import ChemStandardizer, InchiHandling
from matplotlib.ticker import MaxNLocator
from papyrus_scripts.preprocess import (
    consume_chunks,
    keep_contains,
    keep_match,
    keep_organism,
    keep_quality,
)
from tqdm import tqdm
from UniProtMapper import ProtMapper
from venn import venn

from .interfaces import PapyrusLinker

DEFAULT_ORGANISM_LIST = ["Homo sapiens", "Mus musculus", "Rattus norvegicus"]


class PapyrusCompoundLinker(PapyrusLinker):
    """Class to link a bioactivity dataset to the Papyrus dataset.

    To download the dataset, check download_papyrus from papyrus_scripts.
    """

    def __init__(
        self,
        df: pd.DataFrame,
        smiles_col: str,
        bioactivity_col: str,
        papyrus_path: Optional[str] = None,
        papyrus_kwargs: dict = {},
        njobs: int = 8,
    ) -> None:
        """Initialize the PapyrusDataLinker class.

        Args:
            df: a pandas dataframe with the bioactivity data to be linked.
            smiles_col: the name of the column containg the SMILES strings.
            bioactivity_col: name of the column containing the bioactivity labels.
            papyrus_path: path to the papyrus dataset if a subset is already saved.
            papyrus_kwargs: keyword arguments to be passed to
                `papyrus_scripts.read_papyrus`.

        Note: if the papyrus datset is saved under a custom `source_path`, this needs
        to be provided within papyrus_kwargs.
        """
        super().__init__(df, papyrus_path, papyrus_kwargs, njobs)
        self.df = df
        self.name = None  # Attribute used by save and fromFile
        self.smiles_col = smiles_col
        self.bioactivity_col = bioactivity_col
        self.linkerConfig = {
            "papyrusThreshold": None,
            "targetSets": dict(),
            "papyrus_kwargs": papyrus_kwargs,
            "target_cols": list(),
        }
        self.connectivity_name_dict = None
        self.all_targets = None
        self.bioactivity_dict = None

    def _convertTargCols(self, convert_to: str):
        """Function to toggle the target columns as either str or list
        representation.

        Args:
            convert_to: ["list", "string"]
        """
        if convert_to == "list":
            for col in self.linkerConfig["target_cols"]:
                self.df[col] = (
                    self.df[col]
                    .fillna("")
                    .apply(lambda x: x.split(";"))
                    .replace({"": np.nan})
                )
        if convert_to == "string":
            for col in self.linkerConfig["target_cols"]:
                self.df[col] = (
                    self.df[col]
                    .fillna("")
                    .apply(lambda x: ";".join(x))
                    .replace({"": np.nan})
                )

    def applyPapyrusFilters(
        self,
        ids: Union[List[str], str],
        organisms: List[str] = DEFAULT_ORGANISM_LIST,
        min_quality: str = "medium",
    ):
        """Applies several filters to the papyrus dataset. The filters are:
        [keep_organism, keep_quality, keep_type, keep_contains].

        Args:
            ids: a list 'connectivity' or 'accession' ids to filter the papyrus dataset.
            organisms: Organisms to keep. Defaults to ["Homo sapiens", "Mus musculus",
                "Rattus norvegicus"].
            min_quality: minimum quality of the dataset. For mor information check
                (Béquignon et al., 2023). Defaults to "medium".

        Returns:
            pd.DataFrame: a filtered papyrus dataset.
        """
        sample_data = self.papyrusDataset
        filter1 = keep_organism(
            data=sample_data,
            protein_data=self.papyrusProteins,
            organism=organisms,
            generic_regex=True,
        )
        filter2 = keep_quality(data=filter1, min_quality=min_quality)
        # Removed because takes too long to run
        # filter3 = keep_type(data=filter2, activity_types=activity_types)
        filter3 = keep_match(data=filter2, column="connectivity", values=ids)
        filtered_data = consume_chunks(filter3, progress=True, total=60)
        self.papyrusDataset = filtered_data
        print("Size of the Papyrus dataframe:", len(filtered_data))
        print(f"{len(filtered_data['connectivity'].unique())}unique connectivities.")
        print(f"{len(filtered_data['accession'].unique())}unique accessions.")
        print(filtered_data.columns)
        return filtered_data

    def prepareDataset(
        self,
        progress: bool = True,
        verbose: bool = False,
    ):
        """A function to prepare the user's dataset for linking with papyrus.

        Args:
            standardize_smiles: whether to standardize the smiles usingthe papyrus
                structure pipeline. Defaults to True.
        """
        standardizer = ChemStandardizer(
            "papyrus",
            n_jobs=self.njobs,
            progress=progress,
            from_smi=True,
            rdkit_loglevel=("critical" if not verbose else "warning")
        )
        connect_calculator = InchiHandling(
            convert_to="connectivity",
            n_jobs=self.njobs,
            progress=progress,
            from_smi=True,
            rdkit_loglevel=("critical" if not verbose else "warning")
        )
        print("Standardizing SMILES...")
        std_smiles = standardizer(self.df[self.smiles_col].tolist())
        print("Adding SMILES connectivities...")
        self.df = self.df.assign(
            papyrus_smiles=std_smiles, connectivity=connect_calculator(std_smiles)
        ).replace({None: np.NaN})
        print(
            f"{std_smiles.count(None)} failed to be standardized and will be dropped."
        )
        todrop_idxs = self.df.query("papyrus_smiles.isna()").index.tolist()
        self.df = self.df.drop(index=todrop_idxs).reset_index(drop=True)
        return self.df

    def linkData(self, string_type: str = "target_only", drop: bool = True):
        """Link the user's dataset to papyrus. The result is a dataframe with the target
        information for several pChEMBL concentration thresholds

        Args:
            string_type: how the target information is displayed. It's also possible to
                return the median values by setting this to "with_affinity".
                Defaults to "target_only".
            drop: drop the columns where the connectivities are not found in papyrus.

        Raises:
            AttributeError: for invalid `string_type` parameter.

        Returns:
            A dataframe with the target information for several pChEMBL concentrations.
        """
        string_type = string_type.lower()
        if string_type not in ["target_only", "with_affinity"]:
            raise AttributeError(
                "`string_type` parameter not available."
                "Use 'target_only' or 'with_affinity'"
            )

        papyrus_df = self.papyrusDataset
        user_df = self.df.query(
            'connectivity.isin(@papyrus_df["connectivity"].unique())'
        )
        th_range = np.arange(6.0, 9.5, 0.5)
        colnames = [f"Targets_{th:.1f}" for th in th_range]
        self.linkerConfig["target_cols"] = colnames
        results_df = pd.DataFrame(columns=["connectivity"]).assign(
            connectivity=papyrus_df["connectivity"].unique()
        )
        for th, col in tqdm(zip(th_range, colnames), total=len(th_range)):
            subset = papyrus_df.query("pchembl_value_Median > @th")
            len_connects = len(subset["connectivity"].unique())
            if string_type == "target_only":
                targets = (
                    subset.groupby("connectivity")["target_id"]
                    .unique()
                    .reset_index()
                    .rename(columns={"target_id": col})
                )
                assert len(targets) == len_connects
            elif string_type == "with_affinity":
                targets = (
                    subset.assign(
                        id_activ=lambda x: list(
                            x["target_id"] + "~" + x["pchembl_value_Median"].astype(str)
                        )
                    )
                    .groupby("connectivity")["targ_affinity"]
                    .unique()
                    .reset_index()
                )
                assert len(targets) == len_connects
            results_df = pd.merge(results_df, targets, on="connectivity", how="outer")
        # Assign bioactivity class from original dataset
        c_map = dict(zip(user_df["connectivity"], user_df[self.bioactivity_col]))
        results_df = results_df.assign(
            bioactivity=lambda x: x["connectivity"].map(c_map)
        )
        self.df = self.df.merge(results_df, on="connectivity", how="outer")
        self._convertTargCols(convert_to="string")
        # dropping the ones not found in the papyrus dataset
        not_found = np.where(self.df["bioactivity"].isna())[0]
        print(f"{len(not_found)} compounds not found in papyrus dataset.")
        if drop:
            print("Dropping compounds...")
            self.df = (
                self.df.reset_index(drop=True)
                .drop(index=not_found)
                .reset_index(drop=True)
            )
        return results_df

    def save(self, name: str, save_path: str = None, force=False):
        """Save the dataset to the save_path attribute. The directory
        will contain the following files:
        - The original dataset - Optional additional column with standardised smiles;
        - A subset of the papyrus dataset with the desired bioactivity;
        - A json file with the configuration necessary to load the dataset using the
            @staticmethod `fromFile` method.
        """
        self.name = name
        targcol_regex = re.compile(r"Targets_\d+\.\d+$")
        if not self.df.columns.str.contains(targcol_regex).any():
            raise AttributeError(
                "No column with targets found. Please run `linkData()`"
            )
        if save_path is None:
            save_path = Path(name).resolve()
        else:
            save_path = (Path(save_path) / name).resolve()
        papyrus_path = save_path / "papyrus_subset.csv"
        user_data_path = save_path / f"{name}.csv"
        json_path = save_path / "papyrus_linker_config.json"
        # Make target sets json serializable
        targ_sets = self.linkerConfig["targetSets"]
        targ_sets = {k: list(v) for k, v in targ_sets.items()}
        self.linkerConfig["targetSets"] = targ_sets

        if not save_path.exists():
            print(f"Creating directory {save_path}...")
        elif not force:
            print(f"Directory {save_path} already exists.")
            print("To overwrite the files, set force=True.")
            return
        save_path.mkdir(parents=True, exist_ok=True)
        self.linkerConfig.update(
            {
                "bioactivity_col": self.bioactivity_col,
                "smiles_col": self.smiles_col,
                "papyrus_path": str(papyrus_path),  # enable json serialization
                "user_data_path": str(user_data_path),  # enable json serialization
                "njobs": self.njobs,
            }
        )
        self.df.to_csv(user_data_path, index=False)
        self.papyrusDataset.to_csv(papyrus_path, index=False)
        with open(json_path, "w") as f:
            json.dump(self.linkerConfig, f, indent=4)
        return

    @staticmethod
    def fromFile(config_path: str):
        """Instantiate a PapyrusLinker object from a json configuration file.

        Args:
            config_path: either be the direct path to the json configuration
        file or to the directory where the configuration file is located.

        Raises:
            AttributeError: if the configuration file is not found or if there
                are more than one configuration file in the directory.

        Returns:
            a `papyrusCompoundLinker` object.
        """
        config_path = Path(config_path)
        if config_path.is_dir():
            configs = sorted(list(config_path.glob("*_config.json")))
            if len(configs) > 1:
                raise AttributeError(
                    "More than one configuration file found. Please specify the file."
                )
            else:
                config_path = configs[0]
        if config_path.exists():
            with open(config_path, "r") as f:
                papyrus_config = json.load(f)
            linker = PapyrusCompoundLinker(
                df=pd.read_csv(papyrus_config["user_data_path"]),
                smiles_col=papyrus_config["smiles_col"],
                bioactivity_col=papyrus_config["bioactivity_col"],
                papyrus_path=papyrus_config["papyrus_path"],
                papyrus_kwargs=papyrus_config["papyrus_kwargs"],
                njobs=papyrus_config["njobs"],
            )
            linker.linkerConfig = papyrus_config
            linker.getTargetSets(targ_col=papyrus_config["papyrusThreshold"])
            return linker
        else:
            raise AttributeError("Configuration file not found.")

    def getTargetSets(
        self, targ_col: str = "Targets_6.5", rm_mutation: bool = True
    ) -> dict:
        """Take the dataframe returned by linkData() and creates sets
        that are stored as attributes. This function is a prerequisite for
        plotting the venn diagram of the bioactivities.

        Args:
            targ_col: column with discrete active targets considering a bioactivity
                threshold. Defaults to "Targets_6.5".
            rm_mutation: to remove mutation info from `target_id`. Defaults to True.

        Raises:
            AttributeError: if the dataframe doesn't contain targets following the
            regex for targ_col.

        Returns:
            A dictionary with the sets for plotting.
        """
        self._convertTargCols(convert_to="list")
        self.linkerConfig.update({"papyrusThreshold": targ_col})
        targ_col = targ_col
        targ_col_regex = re.compile(r"Targets_\d+\.\d+$")
        if not self.df.columns.str.contains(targ_col_regex).any():
            raise AttributeError(
                "No column with targets found. Please run `linkData()`"
            )
        df = self.df.copy().assign(bioactivity=lambda x: x[self.bioactivity_col])
        bioactivities = df["bioactivity"].unique()
        bioactiv_dict = dict()

        for bio in bioactivities:
            subset = df[(df["bioactivity"] == bio) & (~df[targ_col].isnull())]
            if rm_mutation:
                targets = set(
                    subset[targ_col].explode().apply(lambda x: x.split("_")[0])
                )
            else:
                targets = set(subset[targ_col].explode())
            if "" in targets:
                targets.remove("")
            bioactiv_dict.update({bio: targets})

        all_targets = set.union(*bioactiv_dict.values())
        self.all_targets = all_targets
        self.bioactivity_dict = bioactiv_dict
        self.linkerConfig["targetSets"].update(bioactiv_dict)
        self._convertTargCols(convert_to="string")
        return bioactiv_dict

    def addProteinInfo(
        self,
        pooling_interval: float = 3,
        total_retries: int = 5,
        backoff_factor: float = 0.25,
    ) -> Tuple[pd.DataFrame, List[str]]:
        """Retrieve protein names, organism, tissue specificity, GO terms and a few
        other UniProt crossreferences and merges that to the papyrus dataframe.

        Args:
            pooling_interval: for UniProtRetriever. Defaults to 3.
            total_retries: for UniProtRetriever. Defaults to 5.
            backoff_factor: for UniProtRetriever. Defaults to 0.25.

        Returns:
            Tuple: Result from the UniProtRetriever and a list of failed accessions.
        """
        mapper = ProtMapper(pooling_interval, total_retries, backoff_factor)
        result_df, failed = mapper(
            ids=self.papyrusDataset["accession"].unique(),
            fields=[
                "protein_name",
                "organism_name",
                "cc_tissue_specificity",
                "go_p",
                "go_c",
                "go_f",
                "xref_pdb",
                "xref_chembl",
            ],
        )
        if failed:
            print(f"Failed to retrieve accessions {failed}")
        # p for parentheses
        allafter_p = re.compile(r"\s\(.*")
        inside_p = re.compile(r"\((.*?)\)")

        result_df = result_df.assign(
            mainName=result_df["Protein names"].str.replace(allafter_p, "", regex=True),
            orgName=result_df["Organism"].str.findall(inside_p).str[0],
            nameOrg=lambda x: x["mainName"] + " (" + x["orgName"] + ")",
        ).rename(columns={"From": "accession"})

        self.papyrusDataset = self.papyrusDataset.merge(
            result_df, on="accession", how="outer"
        )
        return result_df, failed

    def addConnectivityNames(self, name_col) -> dict:
        """Adds a connectivity to name dictionary to enable a more readable
        plot when plotting bioactivities for a single connectivity.

        Args:
            name_col: column with compound name on `self.df`.
        """
        self.connectivity_name_dict = dict(
            zip(self.df["connectivity"], self.df[name_col])
        )
        return self.connectivity_name_dict

    def plotCompounsByTarget(self):
        pass

    def plotTargetsByCompound(
        self,
        connectivity: str,
        bioactivity_col: str = "pchembl_value_Median",
        label_col: str = "nameOrg",
        plot_top=10,
        ax=None,
    ):
        """Plots the bioactivities for a single molecule based on its connectivity.

        Args:
            connectivity: connectivy of the desired compound.
            bioactivity_col: column from self.papyrusData with the bioactivity to be
                plotted. Defaults to "pchembl_value_Median".
            label_col: column from self.papyrusDataset containing the desired label.
                Defaults to "nameOrg".
            plot_top: if desired to plot only `n` top-activity compounds.
                Defaults to 10.
            ax: plt.Axes object to plot on. Defaults to None.

        Returns:
            ax: plt.Axes object with the plot.
        """
        connectivity = connectivity

        def set_ax_params(plot_df, comp_name, ax):
            ax.set_title(
                f"{comp_name}'s targets and respective affinities\n"
                "Median pChEMBL value - source: Papyrus"
            )
            ax.grid(axis="x", alpha=0.5)
            ax.set_axisbelow(True)
            maxval = np.ceil(plot_df.max())
            if maxval > 8:
                ax.set_xlim(0, int(maxval))
            else:
                ax.set_xlim(0, 8)
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            ax.set_ylabel("")
            return ax

        toplot_df = self.papyrusDataset.query(
            "connectivity == @connectivity"
        ).set_index(label_col)[bioactivity_col]

        if plot_top is not None:
            if len(toplot_df) > plot_top:
                print(
                    f"More than {plot_top} affinities to plot."
                    f" Will plot only the top {plot_top}"
                )
                # To take the top n values we sort it with the highest first
                toplot_df = toplot_df.sort_values(ascending=False)[:plot_top]
        toplot_df = toplot_df.sort_values(ascending=True)
        if ax is None:
            ax = plt.gca()
        ax = toplot_df.plot.barh(edgecolor="black", color="lightblue", ax=ax)
        # If user already provides the id_mapping, no need to query uniprot
        if self.connectivity_name_dict is not None:
            comp_name = self.connectivity_name_dict.get(connectivity, connectivity)
        else:
            try:
                self.addConnectivityNames("Compound")
            except KeyError:
                print(
                    "No connectivity to name mapping found. "
                    "Run self.addConnectivityNames() first."
                )
                return
        th = float(self.linkerConfig["papyrusThreshold"].split("_")[1])
        ax.axvline(th, color="red", linestyle="--")

        # set ax parameters preference
        ax = set_ax_params(plot_df=toplot_df, comp_name=comp_name, ax=ax)
        return ax

    def bioactivityVenn(
        self,
        toplot_dict: dict = None,
        ax: Optional[plt.Axes] = None,
        cmap="plasma",
        **kwargs,
    ):
        """Plot the output from getTargetSets() as a venn diagram.

        Args:
            toplot_dict: if left as none will use self.bioactivity_dict.
                Defaults to None.
            ax: custom axis to place the plot. Defaults to None.
            kwargs: keyword arguments to pass to the venn function.

        Returns:
            plt.Axes object
        """
        if toplot_dict is None:
            toplot_dict = self.bioactivity_dict
        if ax is None:
            ax = plt.gca()
        ax = venn(
            toplot_dict,
            cmap=cmap,
            alpha=0.6,
            legend_loc="lower right",
            ax=ax,
            **kwargs,
        )
        return ax


class PapyrusProteinLinker(PapyrusLinker):
    """Class to link a bioactivity dataset to the Papyrus dataset.

    To download the dataset, check download_papyrus from papyrus_scripts.
    """

    def __init__(
        self,
        df: pd.DataFrame,
        protein_col: str,
        papyrus_path: Optional[str] = None,
        papyrus_kwargs: dict = {},
        njobs: int = 8,
    ) -> None:
        """Initialize the PapyrusDataLinker class.

        Args:
            df: a pandas dataframe with the bioactivity data to be linked.
            protein_col: the column containing the Uniprot accessions.
            papyrus_path: path to the papyrus dataset if a subset is already saved.
            papyrus_kwargs: keyword arguments to be passed to
                `papyrus_scripts.read_papyrus`.

        Note: if the papyrus datset is saved under a custom `source_path`, this needs
        to be provided within papyrus_kwargs.
        """
        super().__init__(df, papyrus_path, papyrus_kwargs, njobs)
        self.protein_col = protein_col

    def applyPapyrusFilters(
        self,
        ids: Union[List[str], str],
        organisms: List[str] = DEFAULT_ORGANISM_LIST,
        min_quality: str = "medium",
    ):
        """Applies several filters to the papyrus dataset. The filters are:
        [keep_organism, keep_quality, keep_contains]. The keep_contains filter is
        applied to filter the dataset based on uniprot IDs.

        Args:
            ids: a list 'connectivity' or 'accession' ids to filter the papyrus dataset.
            organisms: Organisms to keep. Defaults to ["Homo sapiens", "Mus musculus",
                "Rattus norvegicus"].
            min_quality: minimum quality of the dataset. For mor information check
                (Béquignon et al., 2023). Defaults to "medium".

        Returns:
            pd.DataFrame: a filtered papyrus dataset.
        """
        pattern = "|".join(ids)
        sample_data = self.papyrusDataset
        filter1 = keep_organism(
            data=sample_data,
            protein_data=self.papyrusProteins,
            organism=organisms,
            generic_regex=True,
        )
        filter2 = keep_quality(data=filter1, min_quality=min_quality)
        # Removed because takes too long to run
        # filter3 = keep_type(data=filter2, activity_types=activity_types)
        filter3 = keep_contains(
            data=filter2, column="accession", value=pattern, case=True, regex=True
        )
        filtered_data = consume_chunks(filter3, progress=True, total=60)
        self.papyrusDataset = filtered_data
        print("Size of the Papyrus dataframe:", len(filtered_data))
        print(f"unique connectivities: {filtered_data['connectivity']}")
        print(f"unique accessions: {filtered_data['accession'].unique()}")
        print(filtered_data.columns)
        return filtered_data

    def linkData(self, string_type: str = "target_only"):
        """Link the user's dataset to papyrus. The result is a dataframe with the target
        information for several pChEMBL concentration thresholds

        Args:
            string_type: how the target information is displayed. It's also possible to
                return the median values by setting this to "with_affinity".
                Defaults to "target_only".

        Raises:
            AttributeError: for invalid `string_type` parameter.

        Returns:
            A dataframe with the target information for several pChEMBL concentrations.
        """
        protIDs = self.protein_col
        string_type = string_type.lower()
        if string_type not in ["target_only", "with_affinity"]:
            raise AttributeError(
                "`string_type` parameter not available."
                "Use 'target_only' or 'with_affinity'"
            )

        papyrus_df = self.papyrusDataset
        # user_df = self.df.query(f'{protIDs}.isin(@papyrus_df["accession"].unique())')
        # Will I really link this?
        th_range = np.arange(6.0, 9.5, 0.5)
        colnames = [f"Targets_{th:.1f}" for th in th_range]
        results_df = pd.DataFrame().assign(accession=papyrus_df["accession"].unique())
        for th, col in tqdm(zip(th_range, colnames), total=len(th_range)):
            subset = papyrus_df.query("pchembl_value_Median > @th")
            len_connects = len(subset["accession"].unique())
            if string_type == "target_only":
                targets = (
                    subset.groupby("accession")["target_id"]
                    .unique()
                    .reset_index()
                    .rename(columns={"target_id": col})
                )
                assert len(targets) == len_connects
            elif string_type == "with_affinity":
                targets = (
                    subset.assign(
                        id_activ=lambda x: list(
                            x["target_id"] + "~" + x["pchembl_value_Median"].astype(str)
                        )
                    )
                    .groupby("accession")["targ_affinity"]
                    .unique()
                    .reset_index()
                )
                assert len(targets) == len_connects
            results_df = pd.merge(results_df, targets, on="accession", how="outer")
        # Assign bioactivity class from original dataset
        # will this be of any use??
        # c_map = dict(zip(user_df["accession"], user_df[self.bioactivity_col]))
        # results_df = results_df.assign(
        #     bioactivity=lambda x: x["connectivity"].map(c_map)
        # )
        self.df = self.df.assign(accession=lambda x: x[protIDs]).merge(
            results_df, on="accession", how="outer"
        )
        return results_df
