from datetime import datetime
from typing import Callable, Iterable

import numpy as np
from qsprpred.data import QSPRDataset
from qsprpred.logs import logger
from qsprpred.models import QSPRModel
from qsprpred.models.assessment.methods import ModelAssessor
from qsprpred.models.hyperparam_optimization import OptunaOptimization
from qsprpred.models.monitors import HyperparameterOptimizationMonitor


class NewOptunaOptimization(OptunaOptimization):
    def __init__(
        self,
        param_grid: dict,
        model_assessor: ModelAssessor,
        score_aggregation: Callable[[Iterable], float] = np.mean,
        monitor: HyperparameterOptimizationMonitor | None = None,
        n_trials: int = 100,
        n_jobs: int = 1,
    ):
        """Initialize the class for hyperparameter optimization
        of QSPRModels using Optuna.

        Args:
            param_grid (dict):
                search space for bayesian optimization, keys are the parameter names,
                values are lists with first element the type of the parameter and the
                following elements the parameter bounds or values.
            model_assessor (ModelAssessor):
                assessment method to use for the optimization
                (default: CrossValAssessor)
            score_aggregation (Callable):
                function to aggregate the scores of different folds if the assessment
                method returns multiple predictions
            monitor (HyperparameterOptimizationMonitor):
                monitor for the optimization, if None, a BaseMonitor is used
            n_trials (int):
                number of trials for bayes optimization
            n_jobs (int):
                number of jobs to run in parallel.
                At the moment only n_jobs=1 is supported.
        """
        super().__init__(param_grid, model_assessor, score_aggregation, monitor, n_trials, n_jobs)

    def optimize(
        self,
        model: QSPRModel,
        ds: QSPRDataset,
        save_params: bool = True,
        refit_optimal: bool = False,
        study_name: str = None,
        storage_url: str = None,
        continue_existing: bool = False,
    ) -> dict:
        """Bayesian optimization of hyperparameters using optuna.

        Args:
            model (QSPRModel): the model to optimize
            ds (QSPRDataset): dataset to use for the optimization
            save_params (bool):
                whether to set and save the best parameters to the model
                after optimization
            refit_optimal (bool):
                Whether to refit the model with the optimal parameters on the
                entire training set after optimization. This implies 'save_params=True'.
            study_name (str): name for the optuna study. Defaults to None.
            storage_url (str):
                url for the optuna study storage. If None, and study_name is defined,
                a sqlite database with the study_name is created in the current directory.
                Defaults to None.
            continue_existing (bool):
                whether to continue an existing study or start a new one. Defaults to False.


        Returns:
            dict: best parameters found during optimization
        """
        import optuna

        if storage_url is None and study_name is not None:
            storage_url = f"sqlite:///{study_name}.db"

        self.monitor.onOptimizationStart(model, ds, self.config, self.__class__.__name__)

        logger.info("Bayesian optimization can take a while " "for some hyperparameter combinations")
        # create optuna study
        sampler = optuna.samplers.TPESampler(seed=model.randomState)
        if continue_existing:
            try:
                study = optuna.load_study(study_name=study_name, storage=storage_url, sampler=sampler)
                logger.info(f"Loaded existing study '{study_name}'")
            except optuna.exceptions.DuplicatedStudyError:
                study = optuna.create_study(
                    study_name=study_name, storage=storage_url, direction="maximize", sampler=sampler
                )
                logger.info(f"Created new study '{study_name}'")
        else:
            study = optuna.create_study(
                study_name=study_name, storage=storage_url, direction="maximize", sampler=sampler
            )
            logger.info(f"Created new study '{study_name}'")
        logger.info("Bayesian optimization started: %s" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        study.optimize(lambda t: self.objective(t, model, ds), n_trials=self.nTrials, n_jobs=self.nJobs)
        logger.info("Bayesian optimization ended: %s" % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        # save the best study
        trial = study.best_trial
        # log the best study
        logger.info("Bayesian optimization best params: %s" % trial.params)
        # save the best score and parameters, return the best parameters
        self.bestScore = trial.value
        self.bestParams = trial.params

        self.monitor.onOptimizationEnd(self.bestScore, self.bestParams)
        # save the best parameters to the model if requested
        self.saveResults(model, ds, save_params, refit_optimal)
        return study, self.bestParams
