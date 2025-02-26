from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from papyrus_scripts import read_papyrus, read_protein_set
from scipy.stats import median_abs_deviation


class PapyrusLinkerABC(ABC):
    """Class to link a bioactivity dataset to the Papyrus dataset.

    To download the dataset, check download_papyrus from papyrus_scripts.
    """

    @abstractmethod
    def _loadPapyrus(self):
        pass

    @property
    @abstractmethod
    def papyrusProteins(self):
        pass

    @abstractmethod
    def applyPapyrusFilters(self):
        pass

    @abstractmethod
    def separateBioactivityTypes(self):
        pass

    @abstractmethod
    def linkData(self):
        pass

    @abstractmethod
    def save(self):
        pass

    @staticmethod
    @abstractmethod
    def fromFile():
        pass


class PapyrusLinker(PapyrusLinkerABC):
    """Class to link a bioactivity dataset to the Papyrus dataset.

    To download the dataset, check download_papyrus from papyrus_scripts.
    """

    def __init__(
        self,
        df: pd.DataFrame,
        papyrus_path: Optional[str] = None,
        papyrus_kwargs: dict = {},
        njobs: int = 8,
    ) -> None:
        """Initialize the PapyrusDataLinker class.

        Args:
            df: a pandas dataframe with the bioactivity data to be linked.
            papyrus_path: path to the papyrus dataset if a subset is already saved.
            papyrus_kwargs: keyword arguments to be passed to
                `papyrus_scripts.read_papyrus`.

        Note: if the papyrus datset is saved under a custom `source_path`, this needs
        to be provided within papyrus_kwargs.
        """
        self.df = df
        self.papyrus_path = papyrus_path
        self.papyrus_kwargs = {
            "is3d": False,
            "version": "05.6",
            "plusplus": False,
            "chunksize": 1000000,
            **papyrus_kwargs,
        }
        self.papyrusDataset = self._loadPapyrus()
        self.njobs = njobs

    def _loadPapyrus(self):
        if self.papyrus_path is None:
            return read_papyrus(**self.papyrus_kwargs)
        else:
            return pd.read_csv(self.papyrus_path)

    @property
    def papyrusProteins(self):
        if "source_path" in self.papyrus_kwargs.keys():
            source_path = self.papyrus_kwargs["source_path"]
        else:
            source_path = Path(__file__).home() / ".data"
        return read_protein_set(
            source_path=source_path, version=self.papyrus_kwargs["version"]
        )

    def applyPapyrusFilters(self):
        """This method will be different for compound and target linking"""
        return

    def separateBioactivityTypes(self):
        """Separate the bioactivity types found in Papyrus into the following columns:
        Mean_<activity> Median_<activity> SD_<activity> MAD_<activity>.

        Args:
            df: filtered papyrus dataframe with all the bioactivity types.

        Returns:
            pd.DataFrame: filtered papyrus dataframe with separated bioactivity types.
        """
        df = self.papyrusDataset.copy()

        def series_to_arrays(pd_series, data_type: str = "int"):
            if data_type == "int":
                pd_series = pd_series.str.split(";").apply(
                    lambda x: np.array(x).astype(int)
                )
            if data_type == "float":
                pd_series = pd_series.str.split(";").apply(
                    lambda x: np.array(x).astype(float)
                )
            return pd_series

        def compress_bioactivities(pd_dataframe, bioactiv_type):
            """
            After calling `series_to_arrays`, this function will compress the
            values for pchembl_value according to the binary bioactiv_type.
            """
            boolean = pd_dataframe[bioactiv_type]
            values = pd_dataframe["pchembl_value"]
            results = []

            for b, v in zip(boolean, values):
                # If the values are not same length, no need to compress
                if len(b) != len(v):
                    if b[0] == 1:
                        results.append(v)
                    else:
                        results.append(np.NaN)
                else:
                    compressed = np.compress(b, v)
                    # If we have no bioactivity like that, return NaN
                    if compressed.shape[0] == 0:
                        results.append(np.NaN)
                    else:
                        results.append(compressed)
            return results

        bioactiv_types = ["type_IC50", "type_EC50", "type_KD", "type_Ki", "type_other"]

        df = (
            df.assign(type_IC50=series_to_arrays(df["type_IC50"]))
            .assign(type_EC50=series_to_arrays(df["type_EC50"]))
            .assign(type_KD=series_to_arrays(df["type_KD"]))
            .assign(type_Ki=series_to_arrays(df["type_Ki"]))
            .assign(type_other=series_to_arrays(df["type_other"]))
            .assign(
                pchembl_value=series_to_arrays(df["pchembl_value"], data_type="float")
            )
        )
        df = (
            df.assign(type_IC50=compress_bioactivities(df, bioactiv_type="type_IC50"))
            .assign(type_EC50=compress_bioactivities(df, bioactiv_type="type_EC50"))
            .assign(type_KD=compress_bioactivities(df, bioactiv_type="type_KD"))
            .assign(type_Ki=compress_bioactivities(df, bioactiv_type="type_Ki"))
            .assign(type_other=compress_bioactivities(df, bioactiv_type="type_other"))
        )

        func_dict = {
            "mean": lambda x: np.mean(x) if isinstance(x, np.ndarray) else np.NaN,
            "median": lambda x: np.median(x) if isinstance(x, np.ndarray) else np.NaN,
            "std": lambda x: np.std(x) if isinstance(x, np.ndarray) else np.NaN,
            "mad": lambda x: median_abs_deviation(x)
            if isinstance(x, np.ndarray)
            else np.NaN,
        }

        for bioactiv in bioactiv_types:
            name = bioactiv.split("_")[1]
            df[f"Mean_{name}"] = df[bioactiv].apply(func_dict["mean"])
            df[f"Median_{name}"] = df[bioactiv].apply(func_dict["median"])
            df[f"SD_{name}"] = df[bioactiv].apply(func_dict["std"])
            df[f"MAD_{name}"] = df[bioactiv].apply(func_dict["mad"])

        self.papyrusDataset = df
        return df

    def linkData(self):
        """This method will be different for compound and target linking"""
        return

    def save(self):
        """This method will be different for compound and target linking"""
        return

    def fromFile(self):
        """This method will be different for compound and target linking"""
        return

