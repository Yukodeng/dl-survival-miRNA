from __future__ import annotations
import time
import numpy as np
import pandas as pd
from typing import Dict, Any

from sklearn.model_selection import StratifiedKFold
from sksurv.metrics import integrated_brier_score
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.svm import FastKernelSurvivalSVM
from sksurv.ensemble import RandomSurvivalForest, GradientBoostingSurvivalAnalysis

from .common import TunablePipelineBase, ResultsWriterMixin, _model_string, _default_hp_space
from ..utils import dataframe_to_scikitsurv_ds


MODEL_MAP = {
    "coxnet": CoxnetSurvivalAnalysis,
    "rsf": RandomSurvivalForest,
    "ssvm": FastKernelSurvivalSVM,
    "sgb": GradientBoostingSurvivalAnalysis
}


class MLSurvivalPipeline(TunablePipelineBase, ResultsWriterMixin):
    """
    ML suvival risk prediction models with built-in hyperparameter tuning (optuna). 

    Args:
        model_type (str): Type (``coxnet``, ``rsf``,  ``ssvm``, ``sgb``)
        sim_type (str): Simulation condition code (e.g., ``BE00Asoo00_TC_linear-moderate``)
        train_df (pd.DataFrame): Training data
        test_df (pd.DataFrame): Testing data
        time_col (str): Default ``time``
        status_col (str): Default ``status``
        batch_col (str): Default ``batch.id``
        hyperparameters (Dict[str), Any] | None = None): Optuna-compatible hyperparameter search space
        storage_url (str): Default ``sqlite:///survmodels-hp-log.db``
        is_stratified (bool): Whether to apply 
    """
    
    def __init__(
        self,
        model_type: str,
        train_df: pd.DataFrame,
        test_df: pd.DataFrame,
        *,
        batchNormType: str | None = None,
        dataName: str | None = None,
        is_stratified: bool = False,
        time_col: str = "time",
        status_col: str = "status",
        batch_col: str = "batch.id",
        hyperparameters: Dict[str, Any] | None = None,
        storage_url: str = "sqlite:///survmodels-hp-log.db"
    ):
        modelString = _model_string(model_type, is_stratified)
        hyperparameters = _default_hp_space(model_type) if hyperparameters is None else hyperparameters
        super().__init__(
            batchNormType=batchNormType,
            dataName=dataName, 
            modelString=modelString,
            hyperparameters=hyperparameters,
            storage_url=storage_url
            )
        self.model_type = model_type
        self.train_df = train_df
        self.test_df = test_df
        self.time_col = time_col
        self.status_col = status_col
        self.batch_col = batch_col
        self.is_stratified = is_stratified
        self._best_params = {}

    # ------------------------------------------------------
    #  helper to conditionally remove/encode batch as features
    # ------------------------------------------------------
    def _with_batch_features(self, df: pd.DataFrame) -> pd.DataFrame:
        if self.batch_col not in df.columns:
            return df
        
        # ----- non-stratified (remove batch id ) ------
        if not self.is_stratified:
            return df.drop(columns=[self.batch_col])
        
        # ---- stratified (one-hot encode batch id) ----
        ## NOTE: implement naive approach of batch stratification by including batch.id as predictor
        return pd.get_dummies(
            df,
            columns=[self.batch_col],
            prefix="batch",
            drop_first=False,
            dtype=float)
        
    # --------------------------------------------------
    #  Optuna objective (k-fold mean validation C-index) 
    # --------------------------------------------------
    def _objective(self, df, n_splits, params, fixed_params=None, **_):
        
        full_params = {**fixed_params, **params} if fixed_params else params
        model = MODEL_MAP[self.model_type](**full_params)
        kf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)

        # stratify by "status_batch" for balanced assignment
        stratify_labels = df[[self.status_col, self.batch_col]].astype(str).agg("_".join, axis=1).values
        
        df = df._with_batch_features(df)
        X, y = dataframe_to_scikitsurv_ds(df, time_col=self.time_col, status_col=self.status_col)

        scores = []
        for train_idx, val_idx in kf.split(df, stratify_labels):
            model.fit(X[train_idx], y[train_idx])
            scores.append(model.score(X[val_idx], y[val_idx]))
            
        return np.mean(scores)
    
    # ---------------------------------------------
    # Train/Eval once using self.get_best_params()
    # ---------------------------------------------
    def train(
        self,
        params: Dict[str, Any] | None = None,
        *,
        train_df: pd.DataFrame | None = None,
        test_df: pd.DataFrame | None = None) -> Dict[str, float]:
        
        df_tr = self.train_df.copy() if train_df is None else train_df
        df_te = self.test_df.copy() if test_df is None else test_df
        n_train = df_tr.shape[0]
        
        Xtr, ytr = dataframe_to_scikitsurv_ds(df_tr._with_batch_features(df_tr))
        Xte, yte = dataframe_to_scikitsurv_ds(df_te._with_batch_features(df_te))

        # =================== Train Model ====================
        params = params or self.get_best_params()
        
        # Hard Override of n_estimators for RSF for large subsets
        if "rsf" in self.modelString and n_train > 5000:
            params["n_estimators"] = 50
            print(f"[INFO] Overriding {self.modelString} n_estimators={params['n_estimators']} for n={n_train}")
                     
        model = MODEL_MAP[self.model_type](**params)
        start = time.time()
        model.fit(Xtr, ytr)
        duration = float(round(time.time() - start, 2))
        
        # ==================== Evaluation ====================
        # C-index -------------------------
        tr_c = float(model.score(Xtr, ytr))
        te_c = float(model.score(Xte, yte))

        # IBS -----------------------------
        try:
            tmin = np.ceil(max(ytr[self.time_col].min(), yte[self.time_col].min()))
            tmax = np.floor(min(ytr[self.time_col].max(), yte[self.time_col].max()))
            times = np.linspace(tmin, tmax, 20)
            
            s_tr = np.asarray([[fn(t) for t in times] for fn in model.predict_survival_function(Xtr)])
            s_te = np.asarray([[fn(t) for t in times] for fn in model.predict_survival_function(Xte)])
            
            tr_brier = float(integrated_brier_score(ytr, ytr, s_tr, times))
            te_brier = float(integrated_brier_score(ytr, yte, s_te, times))
        except Exception:
            tr_brier = te_brier = float("nan")

        return duration, tr_brier, te_brier, tr_c, te_c
