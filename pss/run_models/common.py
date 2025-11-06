from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, Sequence
import os
from datetime import datetime
import pandas as pd
import numpy as np
import optuna


def _model_string(model_type: str, is_stratified: bool) -> str:
    return f"{'stratified-' if is_stratified else ''}{model_type}"

def _default_hp_space(model_type: str) -> dict:
    """Default Optuna search spaces per ML/DL model."""
    if model_type == "coxnet":
        return {
            "alphas": {"type": "categorical", "choices": [[2.0**i] for i in np.arange(-6, 3, 2)]},
            "l1_ratio": {"type": "float", "low": 0.1, "high": 1.0},
        }
    if model_type == "ssvm":
        return {
            "alpha": {"type": "float", "low": 2.0**-5, "high": 2.0**5, "log": True},
            "kernel": {"type": "categorical", "choices": ["linear", "poly", "rbf"]},
            "rank_ratio": {"type": "float", "low": 0.0, "high": 1.0},
        }
    if model_type == "rsf":
        return {
            "n_estimators": {"type": "categorical", "choices": [50, 100, 500]},
            "min_samples_split": {"type": "int", "low": 3, "high": 15},
        }
    if model_type == "sgb":
        return {
            "learning_rate": {"type": "float", "low": 0.001, "high": 1},
            "n_estimators": {"type": "categorical", "choices": [50, 100, 500]},
            "min_samples_split": {"type": "int", "low": 3, "high": 15},
        }
    if model_type == "deepsurv-torch":
        return {
            "num_nodes": {"type": "categorical",
                          "choices": [  
                                [128], [64], [32],
                                [128, 64], [64, 32], [32, 16],
                                [64, 64, 32], [32, 32, 16]
                            ]},
            "dropout": {"type": "float", "low": 0.1, "high": 0.5},
            "learning_rate": { "type": "float", "low": 0.0001, "high": 0.005, "log": True},
            "weight_decay": {"type": "float",  "low": 0.0001,"high": 0.05, "log": True},
            "batch_size": { "type": "categorical",  "cahoices": [32, 64, 128]}
        }
    raise ValueError(f"Unknown `model_type`: {model_type!r}. Specify one of 'deepsurv-torch', 'coxnet', 'rsf', 'ssvm', 'sgb'.")


class TunablePipelineBase(ABC):
    """
    A reusable base class for ML/DL training pipelines that share:
      - Hyperparameter search space
      - Optuna study management (create/load, reuse best params)
      - trial threshold logic to skip re-tuning
      - Other common optuna interface: tune_hyperparameters(), set_best_params(), get_best_params()
    Subclasses must implement:
      - _objective(trial, **kwargs) -> float  (maximize)
    Optionally override:
      - _study_name(n_samples) -> str
      - _direction() -> "maximize" | "minimize"
    """

    def __init__(
        self,
        batchNormType: Optional[str] = None,
        dataName: Optional[str] = None,
        modelString: Optional[str] = None,
        hyperparameters: Optional[Dict[str, Any]] = None,
        storage_url: str = "sqlite:///hp-log.db"
    ):
        self.hyperparameters = hyperparameters or {}
        self.storage_url = storage_url
        self._best_params: Dict[str, Any] = {}

        # Scenario identifiers
        self.batchNormType = batchNormType
        self.dataName = dataName
        self.modelString = modelString#_model_string(model_typee, is_stratified)  # set in subclass if computed from flags
        
    # ----------------------------------------------------------------
    # Wrapper for _suggest and _wrapped_objective functions for tuning
    # ----------------------------------------------------------------
    def tune_hyperparameters(
        self,
        df: pd.DataFrame,
        n_samples: int,
        params: Dict[str, any] | None = None,
        n_trials: int = 30,
        n_splits: int = 5,
        n_jobs: int = 1,
        trial_threshold: int = 30,
        timeout: Optional[int] = 2400,
        extra_name_parts: Optional[Sequence[str]] = None,
        **objective_kwargs,
    ) -> optuna.study.Study:
        """
        Perform hyperparameter tuning within the k-fold CV framework using the Optuna library.
        Automatically stores tuning results or retrieves from existing studies to avoid duplicate optimization processes. 
        Defualt objective() is to to maximize the average validation concordance index (C-index) scores.  
        
        Args:
            df: Training dataset.
            n_samples: default = None
            params: Search grid for hyperparameter tuning.
            n_splits: Number of splits for cross-validation.
            n_trials: Number of trials to run for tuning.
            n_jobs: Number of parallel jobs for tuning.
            trial_threshold: Minimum number of trials required to skip tuning.
            timeout: Maximum run time for a tuning study.
        """
        params = self.hyperparameters if params is None else params
        if not params:
            raise ValueError("Hyperparameters search space is empty.")

        study_name = self._study_name(n_samples, extra_name_parts)
        try:          
            study = optuna.load_study(study_name=study_name, storage=self.storage_url)
            self._best_params = study.best_params
            print(f"✅Retrieved best hyperparameters from existing study '{study_name}': {study.best_params}")
            return study
        
        except Exception:
            print(f"⚠️No Optuna study '{study_name}' found. Start hyperparameter tuning...")
            study = optuna.create_study(
                direction=self._direction(),
                storage=self.storage_url,
                study_name=study_name,
                load_if_exists=True
            )
        # Reuse already-completed trials to avoid re-tuning
        successful = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]
        fixed_params = study.best_params if len(successful) >= trial_threshold else {}
        new_space = {k: v for k, v in self.hyperparameters.items() if k not in fixed_params}

        if len(new_space) == 0:
            self._best_params = fixed_params
            print(f"All hyperparameters already tuned: {self._best_params}\nSkipping optimization...")
            return study

        def _suggest(trial: optuna.trial.Trial, space: Dict[str, Any]) -> Dict[str, Any]:
            """ Convert a dictionary of hyperparameter ranges into Optuna-compatible format.
            
            Args:
            - trial (optuna.trial.Trial): The Optuna trial object to sample parameters for.
            - space (Dict): A dictionary with hyperparameter names as keys. Values should be a dictionary specifying the type of parameter ('int', 'float', or 'categorical')
            and the range or choices.
            
            Returns:
            - params (Dict): A dictionary of sampled hyperparameters from the given ranges (continuous) / choices (discrete).
            """
            params = {}
            for name, info in space.items():
                t = info["type"]
                if t == "int":
                    params[name] = trial.suggest_int(name, info["low"], info["high"])
                elif t == "float":
                    params[name] = trial.suggest_float(name, info["low"], info["high"], log=info.get("log", False))
                elif t == "categorical":
                    params[name] = trial.suggest_categorical(name, info["choices"])
                else:
                    raise ValueError(f"Unknown param type: {t}")
            return params

        def _wrapped_objective(trial: optuna.trial.Trial) -> float:
            sampled = _suggest(trial, new_space)
            full = {**fixed_params, **sampled} if fixed_params else sampled
            return self._objective(
                trial=trial,
                df=df,
                params=full,
                n_splits=n_splits,
                **objective_kwargs)

        study.optimize(_wrapped_objective, 
                        n_trials=n_trials, 
                        n_jobs=n_jobs, 
                        timeout=timeout)

        if study.best_params:
            self._best_params = {**fixed_params, **study.best_params}
        else:
            self._best_params = fixed_params

        return study

    def set_best_params(self, params: Dict[str, Any]) -> None:
        self._best_params = dict(params)

    def get_best_params(self) -> Dict[str, Any]:
        # if self._best_params:
        #     return self._best_params
        
        # # default fallback (first choice / lower bound)
        # best = {}
        # for k, v in self.hyperparameters.items():
        #     if isinstance(v, dict) and "choices" in v:
        #         best[k] = v["choices"][0]
        #     elif isinstance(v, dict) and "low" in v:
        #         best[k] = v["low"]
        # self._best_params = best
        return self._best_params


    # ----- To be implemented/overridden by subclasses -----
    @abstractmethod
    def _objective(
        self,
        # trial: optuna.Trial,
        df: pd.DataFrame,
        n_splits: int,
        params: Dict[str, Any],
        **kwargs,
    ) -> float:
        """Return a metric to maximize (e.g., mean C-index across folds)."""
        raise NotImplementedError

    def _study_name(self, n_samples: int, parts: str = []) -> str:
        suffix = "-".join([str(n_samples)] + list(parts or []))
        return f"{self.batchNormType}-{self.dataName}-{self.modelString}-{suffix}"

    def _direction(self) -> str:
        return "maximize"
    


class ResultsWriterMixin:
    """
    Mixin that adds .write(df, ...) to save model results to:
        models/<batchNormType>/<dataName>/<modelString>/model_results_MMDDYY.csv
    """
    def write(self, 
              model_results: pd.DataFrame, *,
              root = "models",
              file_path: str | None = None,
              file_name: str | None = None,
              append: bool = True
    ):
        parts = [
            getattr(self, "batchNormType", None),
            getattr(self, "dataName", None),
            getattr(self, "modelString", None),
        ]
        parts = [p for p in parts if p]
        out_dir = os.path.join(root, *parts) if file_path is None else file_path
        os.makedirs(out_dir, exist_ok=True)

        today = datetime.now().strftime("%m%d%y")
        file_name = f'model_results_{today}.csv' if file_name is None else file_name
        path = os.path.join(out_dir, file_name)

        exists = os.path.exists(path)
        if path.endswith(".txt"):
            model_results.to_csv(path, mode=("a" if (append and exists) else "w"), sep="\t",
                      header=not exists, index=False)
        else:
            model_results.to_csv(path, mode=("a" if (append and exists) else "w"),
                      header=not exists, index=False)