"""Optional QSPRpred split adapters retained from model development.

The committed model-training notebook uses the scikit-learn ``ShuffleSplit``
re-exported by this module. MOOD is therefore not part of the reproducibility
environment; the two MOOD adapters below import it lazily only if explicitly
used.
"""

from typing import Iterable

import numpy as np
import pandas as pd
from qsprpred.data import QSPRDataset
from qsprpred.data.sampling.splits import DataSplit
from qsprpred.logs import logger
from qsprpred.utils.interfaces.randomized import Randomized
from sklearn.model_selection import ShuffleSplit


class QSPRMoodPerimeterSplit(DataSplit, Randomized):
    """Splits dataset in random train and test subsets.

    Attributes:
        testFraction (float):
            fraction of total dataset to testset
        seed (int):
            Random state to use for shuffling and other random operations.
    """

    def __init__(
        self,
        test_fraction=None,
        train_fraction=None,
        n_folds: int = 5,
        dataset: QSPRDataset | None = None,
        seed: int | None = None,
        split_kwargs: dict | None = None,
    ) -> None:
        self.testFraction = test_fraction
        self.trainFraction = train_fraction
        super().setSeed(seed or (dataset.randomState if self.hasDataSet else None))
        self.setDataSet(dataset=dataset)
        self.nFolds = n_folds
        self.split_kwargs = split_kwargs if split_kwargs else {}

    def setDataSet(self, dataset: "MoleculeDataTable"):
        super().setDataSet(dataset)

    def split(
        self, X: np.ndarray | pd.DataFrame, y: np.ndarray | pd.DataFrame | pd.Series
    ) -> Iterable[tuple[list[int], list[int]]]:
        from mood.splitter import PerimeterSplit

        logger.debug(f"Split arguments: {self.split_kwargs}")
        splitter = PerimeterSplit(
            n_clusters=25,
            n_splits=self.nFolds,
            test_size=self.testFraction,
            train_size=self.trainFraction,
            random_state=self.seed,
        )
        for train_indices, test_indices in splitter.split(X.values, y.values):
            yield train_indices, test_indices

# Retained as commented reference because it may be useful in future work. It is
# not used by this study and would otherwise add datamol solely for one method.

# class QSPRMoodScaffoldSplit(DataSplit, Randomized):
#     """Splits dataset in random train and test subsets.

#     Attributes:
#         testFraction (float):
#             fraction of total dataset to testset
#         seed (int):
#             Random state to use for shuffling and other random operations.
#     """

#     def __init__(
#         self,
#         test_fraction=None,
#         train_fraction=None,
#         n_folds: int = 5,
#         dataset: QSPRDataset | None = None,
#         seed: int | None = None,
#         split_kwargs: dict | None = None,
#     ) -> None:
#         self.testFraction = test_fraction
#         self.trainFraction = train_fraction
#         super().setSeed(seed or (dataset.randomState if self.hasDataSet else None))
#         self.setDataSet(dataset=dataset)
#         self.nFolds = n_folds
#         self.split_kwargs = split_kwargs if split_kwargs else {}
#         self.setScaffolds()

#     def setScaffolds(self):
#         # logger.info("Computing scaffolds")
#         ds = self.dataSet
#         if ds.X_ind.empty:  # calculate independent test set
#             smiles = ds.df[ds.smilesCol]
#         else:
#             smiles = ds.df.loc[ds.X.index, ds.smilesCol]
#         n_jobs = ds.nJobs
#         with Parallel(n_jobs=n_jobs) as parallel:
#             scaffolds = parallel(delayed(self._scaffold_from_smi)(smi) for smi in smiles)
#         self._scaffolds = scaffolds

#     @staticmethod
#     def _scaffold_from_smi(smi: str) -> str:
#         return dm.to_smiles(dm.to_scaffold_murcko(dm.to_mol(smi)))

#     def setDataSet(self, dataset: "MoleculeDataTable"):
#         super().setDataSet(dataset)

#     def split(
#         self, X: np.ndarray | pd.DataFrame, y: np.ndarray | pd.DataFrame | pd.Series
#     ) -> Iterable[tuple[list[int], list[int]]]:
#         logger.debug(f"Split arguments: {self.split_kwargs}")
#         splitter = PredefinedGroupShuffleSplit(
#             n_splits=self.nFolds,
#             groups=self._scaffolds,
#             test_size=self.testFraction,
#             train_size=self.trainFraction,
#             random_state=self.seed,
#         )
#         for train_indices, test_indices in splitter.split(X.values, y.values):
#             yield train_indices, test_indices


class QSPRMoodMaxDissimilaritySplit(DataSplit, Randomized):
    """Splits dataset in random train and test subsets.

    Attributes:
        testFraction (float):
            fraction of total dataset to testset
        seed (int):
            Random state to use for shuffling and other random operations.
    """

    def __init__(
        self,
        test_fraction=None,
        train_fraction=None,
        n_folds: int = 5,
        dataset: QSPRDataset | None = None,
        seed: int | None = None,
        split_kwargs: dict | None = None,
    ) -> None:
        self.testFraction = test_fraction
        self.trainFraction = train_fraction
        super().setSeed(seed or (dataset.randomState if self.hasDataSet else None))
        self.setDataSet(dataset=dataset)
        self.nFolds = n_folds
        self.split_kwargs = split_kwargs if split_kwargs else {}

    def setDataSet(self, dataset: "MoleculeDataTable"):
        super().setDataSet(dataset)

    def split(
        self, X: np.ndarray | pd.DataFrame, y: np.ndarray | pd.DataFrame | pd.Series
    ) -> Iterable[tuple[list[int], list[int]]]:
        from mood.splitter import MaxDissimilaritySplit

        logger.debug(f"Split arguments: {self.split_kwargs}")
        splitter = MaxDissimilaritySplit(
            n_clusters=25,
            n_splits=self.nFolds,
            test_size=self.testFraction,
            train_size=self.trainFraction,
            random_state=self.seed,
        )
        for train_indices, test_indices in splitter.split(X.values, y.values):
            yield train_indices, test_indices


class QSPRShuffleSplit(DataSplit, Randomized):
    """Splits dataset in random train and test subsets.

    Attributes:
        testFraction (float):
            fraction of total dataset to testset
        seed (int):
            Random state to use for shuffling and other random operations.
    """

    def __init__(
        self,
        test_fraction=None,
        n_folds: int = 5,
        dataset: QSPRDataset | None = None,
        seed: int | None = None,
        split_kwargs: dict | None = None,
    ) -> None:
        self.testFraction = test_fraction
        super().setSeed(seed or (dataset.randomState if self.hasDataSet else None))
        self.setDataSet(dataset=dataset)
        self.nFolds = n_folds
        self.split_kwargs = split_kwargs if split_kwargs else {}

    def setDataSet(self, dataset: "MoleculeDataTable"):
        super().setDataSet(dataset)

    def split(
        self, X: np.ndarray | pd.DataFrame, y: np.ndarray | pd.DataFrame | pd.Series
    ) -> Iterable[tuple[list[int], list[int]]]:
        logger.debug(f"Split arguments: {self.split_kwargs}")
        splitter = ShuffleSplit(
            n_splits=self.nFolds,
            test_size=self.testFraction,
            random_state=self.seed,
        )
        for train_indices, test_indices in splitter.split(X.values, y.values):
            yield train_indices, test_indices
