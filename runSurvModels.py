import os
import time
import logging
import numpy as np
import pandas as pd
from datetime import datetime

import optuna
from utils import dataframe_to_scikitsurv_ds
from sklearn.model_selection import StratifiedKFold, train_test_split
from sksurv.metrics import integrated_brier_score
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.svm import FastKernelSurvivalSVM
from sksurv.ensemble import RandomSurvivalForest, GradientBoostingSurvivalAnalysis

class SurvivalModelPipeline():
    def __init__(self,
                train_df, 
                test_df,
                batchNormType,
                dataName,
                hyperparameters=None,
                test_size=1000,
                time_col='time', status_col='status', batch_col="batch.id",
                storage_url="sqlite:///survmodels-hp-log.db"
                ):
        self.train_df = train_df
        self.test_df = test_df
        self.test_size = test_size
        self.time_col = time_col
        self.status_col = status_col
        self.batch_col = batch_col
        self.dataName = dataName
        self.batchNormType = batchNormType
        self.hyperparameters = hyperparameters
        self.storage_url = storage_url
        self.best_params = {}
    
    
    def _create_model(self, **kwargs):
        raise NotImplementedError
    
    def write(self, model_results, out_dir=None, fileName=None):
        '''Save model performance results as a txt or csv file.
        '''
        out_dir = os.path.join('models', self.batchNormType, self.dataName, self.modelString) if out_dir is None else out_dir
        os.makedirs(out_dir, exist_ok=True)
        today = datetime.now().strftime("%Y%m%d")
        fileName = f'model_results_{today}.csv' if fileName is None else fileName
        
        if 'txt' in fileName:
            model_results.to_csv(os.path.join(out_dir,fileName), sep='\t')
        elif 'csv' in fileName:
            model_results.to_csv(os.path.join(out_dir,fileName), index=False)
        else:
            print('Please specify a file name with either a txt or csv extension.')
    
    
    def parse_hyperparameters(self, trial, config):
        params = {}
        for name, info in config.items():
            if info["type"] == "int":
                params[name] = trial.suggest_int(name, info["low"], info["high"])
            elif info["type"] == "float":
                params[name] = trial.suggest_float(name, info["low"], info["high"], log=info.get("log", False))
            elif info["type"] == "categorical":
                params[name] = trial.suggest_categorical(name, info["choices"])
        return params
    
    
    def _objective(self, trial, train_df,  n_splits, config=None, fixed_params=None):
        
        # Prepare model hyperparameters setting
        config = self.hyperparameters if config is None else config
        tunable_params = self.parse_hyperparameters(trial, config=config)
        full_params = {**fixed_params, **tunable_params} if fixed_params else tunable_params
        model = self._create_model(**full_params)
        
        kf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
        
        stratify_labels = train_df[[self.status_col, self.batch_col]].astype(str).agg('_'.join, axis=1).values
        
        train_df = train_df.drop(columns=self.batch_col)
        X, y = dataframe_to_scikitsurv_ds(train_df)
        
        val_scores = []
        for train_idx, val_idx in kf.split(train_df, stratify_labels):
            
            X_train, X_val = X[train_idx], X[val_idx]
            y_train, y_val = y[train_idx], y[val_idx]
            
            model.fit(X_train, y_train)
            score = model.score(X_val, y_val)
            val_scores.append(score)
            
        return np.mean(val_scores)
    
    def tune_hyperparameters(self, 
                            train_df, #has batch id column
                            n_samples,
                            n_splits,
                            n_trials,
                            trial_threshold,
                            study_name=None,
                            config=None,
                            n_jobs=1,
                            timeout=7200):
        
        # Create an Optuna study for hyperparameter optimization
        study_name = f"{self.modelString}-{self.batchNormType}-{self.dataName}-{n_samples}" if study_name is None else study_name
        study = optuna.create_study(direction="maximize",
                                    storage=self.storage_url,
                                    study_name=study_name, 
                                    load_if_exists=True)
        
        # Define search space
        config = self.hyperparameters if config is None else config
        if config is None:
            raise ValueError("No hyperparameter search space provided. Set `config` or `self.hyperparameters`.")
        
        # Extract already optimized hyperparameters and prepare new hyperparameter search grid
        successful_trials = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]
        fixed_params = study.best_params if len(successful_trials) >= trial_threshold else {}
        new_config = {k: v for k, v in config.items() if k not in fixed_params}
        if len(new_config) == 0:
            self.best_params = fixed_params
            print(f"All hyperparameters already tuned: {self.best_params}\nSkipping optimization...")
            return study
        
        # Hyperparameter optimization using the objective function
        study.optimize(
            lambda trial: self._objective(trial, train_df,
                                        n_splits=n_splits,
                                        config=new_config,
                                        fixed_params=fixed_params),
            n_trials=n_trials,
            n_jobs=n_jobs,
            timeout=timeout
        )
        
        # Save the best hyperparameters
        self.best_params = {**fixed_params, **study.best_params}
        print(f"Found best hyperparameters: {self.best_params}")
        
        
    def train(self, X_train, y_train, X_test, y_test):
        # =================== Train Model ====================
        model = self._create_model(**self.best_params)
        start = time.time()
        model.fit(X_train, y_train)
        duration = round(time.time() - start, 2)
        
        # ==================== Evaluation ====================
        # Concordance index ----------------
        train_c = model.score(X_train, y_train)
        test_c = model.score(X_test, y_test)
        
        # Integrated Brier score -----------
        try:
            min_t = np.ceil(max(y_train[self.time_col].min(), y_test[self.time_col].min()))
            max_t = np.floor(min(y_train[self.time_col].max(), y_test[self.time_col].max()))
            times = np.linspace(min_t, max_t, 20)

            preds_train = np.asarray([[fn(t) for t in times] for fn in model.predict_survival_function(X_train)])
            preds_test = np.asarray([[fn(t) for t in times] for fn in model.predict_survival_function(X_test)])

            train_brier = integrated_brier_score(y_train, y_train, preds_train, times)
            test_brier = integrated_brier_score(y_train, y_test, preds_test, times)
        except (AttributeError, ValueError):
            train_brier = test_brier = np.nan

        return duration, train_c, test_c, train_brier, test_brier


    def train_over_subsets(self,
                           subset_sizes=[100,500,1000,2000,5000,10000],
                           runs_per_size=[20,20,20,20,20,20],
                           splits_per_size=[3,5,5,10,10,10],
                           trials_per_size=[20,20,20,20,20,20],
                           trial_threshold=20,
                           n_jobs=1,
                           is_tune=False,
                           is_save=False,
                           fileName="model_results.csv"):
        
        model_results = []
        for n, runs, splits, n_trials in zip(subset_sizes, runs_per_size, splits_per_size, trials_per_size):
            
            print(f"Running for training size N={n}...")
            
            run_time, run_train_cind, run_test_cind, run_train_brier, run_test_brier = [],[],[],[],[]
            for run in range(runs):
                train_sub, _ = train_test_split(self.train_df, train_size=n, 
                                                stratify=self.train_df[[self.status_col, self.batch_col]],
                                                shuffle=True, random_state=run)
                test_sub, _  = train_test_split(self.test_df,
                                                train_size=self.test_size,
                                                stratify=self.test_df[[self.status_col, self.batch_col]],
                                                shuffle=True, random_state=run)
                
                if run == 0:
                    if is_tune:     
                        self.tune_hyperparameters(
                            train_sub, n_samples=n, n_splits=splits, n_trials=n_trials, 
                            trial_threshold=trial_threshold, n_jobs=n_jobs
                        )
                    else:
                        try:
                            study_name = f"{self.modelString}-{self.batchNormType}-{self.dataName}-{n}"
                            study = optuna.load_study(study_name=study_name, storage=self.storage_url)
                            best_params = study.best_params
                            print(f"Retrieved best hyperparameters from {study_name}: {best_params}")
                        except Exception:
                            print(f"⚠️No Optuna study found for '{study_name}'. Fall back to hyperparameters from config.")
                            best_params = {}
                            for key, val in self.hyperparameters.items():
                                if isinstance(val, dict) and "choices" in val:
                                    best_params[key] = val["choices"][0]  # Pick first default choice
                                elif isinstance(val, dict) and "low" in val and "high" in val:
                                    best_params[key] = val["low"]  # Pick lower bound as default
                                else:
                                    best_params[key] = val  # directly use the value
                                    
                    # Apply override if defined in self
                    if hasattr(self, "param_override") and self.param_override:
                        best_params.update(self.param_override)
                        print(f"Final parameter set after override: %s", best_params)
                
                ### NOTE: conversion to scikitsurv format assumes no batch.id column
                X_train, y_train = dataframe_to_scikitsurv_ds(train_sub.drop(columns=self.batch_col))
                X_test, y_test = dataframe_to_scikitsurv_ds(test_sub.drop(columns=self.batch_col))
                
                duration, tr_cind, te_cind, tr_brier, te_brier = self.train(X_train, y_train, X_test, y_test)
                model_results.append({
                    "n train": n,
                    "train time": duration,
                    "train C": tr_cind,
                    "test C": te_cind,
                    "train brier": tr_brier,
                    "test brier": te_brier
                })
                run_train_cind.append(tr_cind)
                run_test_cind.append(te_cind)
                run_train_brier.append(tr_brier)
                run_test_brier.append(te_brier)
                run_time.append(duration)
            
            print(
                f"(Avg. runtime: {np.mean(run_time):.2f}s)   |\
                (C-index)  Train: {round(np.mean(run_train_cind),3)}, Test: {round(np.mean(run_test_cind),3)}   |\
                (Brier)  Train: {round(np.nanmean(run_train_brier),3)}, Test: {round(np.nanmean(run_test_brier),3)} (Mean)\n"
            )                
        
        model_results = pd.DataFrame(model_results)  
        if is_save:
            self.write(model_results=model_results, fileName=fileName)
                
        return model_results


class CoxPHElasticNetModel(SurvivalModelPipeline):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.modelString = 'coxnet'
        if self.hyperparameters is None:
            self.hyperparameters = {
                # 'alphas': {'type': 'categorical', 'choices': [[2.0**i] for i in np.arange(-6, 3, 2)]},
                'l1_ratio': {'type': 'float', 'low': 0.1, 'high': 1.0}
            }
    def _create_model(self, **kwargs):
        return CoxnetSurvivalAnalysis(fit_baseline_model=True, **kwargs)


class SVMSurvivalModel(SurvivalModelPipeline):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.modelString = 'svm'
        if self.hyperparameters is None:
            self.hyperparameters = {
                'alpha': {'type': 'float', 'low': 2.0**-5, 'high': 2.0**5, 'log': True},
                'kernel': {'type': 'categorical', 'choices': ['linear', 'poly', 'rbf']},
                'rank_ratio': {'type': 'float', 'low': 0.0, 'high': 1.0}
            }
    def _create_model(self, **kwargs):
        return FastKernelSurvivalSVM(**kwargs)


class RandomSurvivalForestModel(SurvivalModelPipeline):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.modelString = 'rsf'
        if self.hyperparameters is None:
            self.hyperparameters = {
                'n_estimators': {'type': 'categorical', 'choices': [100, 500, 1000]},
                # 'max_depth': {'type': 'int', 'low': 3, 'high': 10},
                'min_samples_split': {'type': 'int', 'low': 3, 'high': 15}
            }
    def _create_model(self, **kwargs):
        return RandomSurvivalForest(**kwargs)

class GradientBoostingSurvivalModel(SurvivalModelPipeline):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.modelString = 'gb'
        if self.hyperparameters is None:
            self.hyperparameters = {
                'learning_rate': {'type': 'float', 'low': 0.001, 'high': 1},
                'n_estimators': {'type': 'categorical', 'choices': [100, 500, 1000]},
                # 'max_depth': {'type': 'int', 'low': 3, 'high': 10},
                'min_samples_split': {'type': 'int', 'low': 3, 'high': 150}
            }
    def _create_model(self, **kwargs):
        return GradientBoostingSurvivalAnalysis(**kwargs)
