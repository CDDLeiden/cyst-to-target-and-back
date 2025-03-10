import os
from pickle import dump
from typing import Any

import numpy as np
import pandas as pd
from qsprpred.logs import logger
from qsprpred.models import SklearnModel
from qsprpred.tasks import ModelTasks
from sklearn.exceptions import NotFittedError
from sklearn.utils.validation import check_is_fitted


class XGBoostModel(SklearnModel):
    def __init__(
        self,
        base_dir: str,
        alg=None,
        name: str | None = None,
        parameters: dict | None = None,
        autoload: bool = True,
        random_state: int | None = None,
    ):
        """Initialize SklearnModel model.

        Args:
            base_dir (str): base directory for model
            alg (Type): sklearn model class
            name (str): customized model name
            parameters (dict): model parameters
            autoload (bool): load model from file
            random_state (int): seed for the random state
        """
        super().__init__(base_dir, alg, name, parameters, autoload, random_state)
        self.fit_verbose = False
        # Initialize models with defined parameters:
        try:
            # check if alg can be initialized with parameters
            if self.parameters is not None:
                self.alg(**self.parameters)
            else:
                self.alg()
        except:
            logger.error(f"Cannot initialize alg {self.alg} with parameters {self.parameters}.")
            raise
        # set parameters if defined
        if (self.parameters not in [None, {}]) and hasattr(self, "estimator") and self.estimator is not None:
            try:
                check_is_fitted(self.estimator)
            except NotFittedError:
                self.estimator.set_params(**self.parameters)
        # log some things
        logger.info("parameters: %s" % self.parameters)
        logger.debug(f'Model "{self.name}" initialized in: "{self.baseDir}"')
        self.eval_set = None  # for early stopping

    @property
    def supportsEarlyStopping(self) -> bool:
        """Whether the model supports early stopping or not."""
        return True

    def set_validation_set(
        self,
        X_test: pd.DataFrame | np.ndarray = None,
        y_test: pd.DataFrame | np.ndarray = None,
        fit_verbose: bool = False,
    ):
        """Set the validation set for the model to use in early stopping mode.

        Args:
            X_test: values with the validation features. Defaults to None.
            y_test: values with the validation target property. Defaults to None.
            fit_verbose: whether to fit with verbose or not. This will affect the callback
                passed to the model, for example. Defaults to False.
        """
        self.fit_verbose = fit_verbose
        self.eval_set = [(X_test, y_test)]

    def _check_model_is_supported(self) -> None:
        """Check if the model is supported."""
        if self.model.__class__.__name__ not in ["XGBClassifier", "XGBRegressor"]:
            raise ValueError(f"Model {self.model.__class__.__name__} not supported")

    def fit(
        self,
        X: pd.DataFrame | np.ndarray,
        y: pd.DataFrame | np.ndarray,
        estimator: Any = None,
        mode: Any = None,
        monitor: None = None,
        **kwargs,
    ):
        if any(["callbacks" in kwargs, "early_stopping_rounds" in kwargs]):
            logger.warning(
                "When setting early stopping, please set the validation set for early "
                "stopping using the `set_validation_set` method. Otherwise, early stoppoing "
                "won't work."
            )
        # check for incompatible tasks
        try:
            if self.task == ModelTasks.MULTITASK_MIXED:
                raise ValueError(
                    "MultiTask with a mix of classification and regression tasks "
                    "is not supported for sklearn models."
                )
            if self.task == ModelTasks.MULTITASK_MULTICLASS:
                raise NotImplementedError(
                    "At the moment there are no supported metrics "
                    "for multi-task multi-class/mix multi-and-single class classification."
                )
        except TypeError:
            pass
        estimator = self.estimator if estimator is None else estimator
        X, y = self.convertToNumpy(X, y)
        # sklearn models expect 1d arrays
        # for single target regression and classification
        if not self.task.isMultiTask():
            y = y.ravel()
        return estimator.fit(X, y, eval_set=self.eval_set, verbose=self.fit_verbose)

    def save(self, save_estimator=False):
        """Save model to file.

        Args:
            save_estimator (bool):
                Explicitly save the estimator to file, if `True`.
                Note that some models may save the estimator by default
                even if this argument is `False`.

        Returns:
            str:
                absolute path to the metafile of the saved model
            str:
                absolute path to the saved estimator, if `include_estimator` is `True`
        """
        os.makedirs(self.outDir, exist_ok=True)
        try:
            meta_path = self.toFile(self.metaFile)
        except:
            logger.error(f"Could not save model to {self.metaFile}")
            pass
        if save_estimator:
            try:
                est_path = self.saveEstimator()
            except TypeError:
                # save it as a pickle file
                est_path = self.toFile(f"{self.outPrefix}.pkl")
                with open(est_path, "wb") as f:
                    dump(self.estimator, f)
