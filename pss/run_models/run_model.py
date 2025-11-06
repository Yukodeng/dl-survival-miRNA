from __future__ import annotations
from typing import Sequence
import numpy as np
import pandas as pd
import optuna
from sklearn.model_selection import train_test_split


def train_over_subsets(
    pipeline,
    *,
    test_size: int = 1000,
    subset_sizes: Sequence[int]=[100,500,1000,2000,5000,10000],
    runs_per_size: Sequence[int]=[20,20,20,20,20,20],
    splits_per_size: Sequence[int]=[3,5,5,10,10,10],
    trials_per_size: Sequence[int]=[30,30,30,30,30,30],
    run_seeds: Sequence[int] | None = None, ## update 7/28 (YD): specify seeds for iterations
    trial_threshold: int = 30,
    n_jobs: int = 1,
    is_tune: bool = False,
    is_save: bool = False,
    file_path: str | None = None,
    file_name: str | None = None
    ) -> pd.DataFrame:
    """
    Standalone simulation executer for both MLSurvivalPipeline and DeepSurvPipeline that tunes hyperparameters once per N, then uses the best params for training iterations.

    Args:
        pipeline: an instance of MLSurvivalPipeline or DeepSurvPipeline adapter.
        test_size (int): default ``1000``.
        subset_sizes (List): [100,500,1000,2000,5000,10000].
        runs_per_size (List): [20,20,20,20,20,20].
        splits_per_size (List): [3,5,5,10,10,10].
        trials_per_size (List): [30,30,30,30,30,30].
        run_seeds (List): list of random seeds (int) for the iteration (length should match with the rest List arguments).
        trial_threshold: Number of existing trials in a study to bypass tuning process, default ``30``.
        n_jobs (int): default ``1``.
        is_tune (bool): defualt ``False``.
        is_save (bool) defualt ``False``.
        file_path (str): Folder where the simulaiton results will be saved, default ``None``.
            If ``is_save`` is ``True``, results by default will be saved to `models/<batchNormType>/<dataName>/<modelString>/`
        file_name (str): Name of the reslut file, default ``None``.
            If ``is_save`` is ``True``, results by default will be saved as ``model_results_MMDDYY.csv``
    
    Returns:
        pd.DataFrame: Survival modeling results across conditions.
            Columns:
                Name: 'n train', dtype: int
                Name: 'train time', dtype: float
                Name: 'train C', dtype: float
                Name: 'test C', dtype: float
                Name: 'train brier', dtype float
                Name: 'test brier', dtype: float
            
    Examples
    --------
    """
    train_df, test_df = pipeline.train_df, pipeline.test_df

    model_results = []
    for n, n_runs, n_splits, n_trials in zip(subset_sizes, runs_per_size, splits_per_size, trials_per_size):
        
        print(f"Running for training size N={n}...")
            
        # Apply flexible early stopping based on training sample size
        if hasattr(pipeline, "early_stop_per_size"):
            stop_params = pipeline.early_stop_per_size.get(str(n), {})
            pipeline.patience = stop_params.get("patience", 30)
            pipeline.min_delta = stop_params.get("min_delta", 1e-3)
        
        # ======== Run multiple iterations with custom random seeds ========
        seeds = run_seeds if run_seeds is not None else list(range(n_runs))
        if len(seeds) != n_runs:
            print("Length of run_seeds must match number of runs. Use default seeds.")
            seeds = list(range(n_runs))
        
        run_time, run_train_cind, run_test_cind, run_train_brier, run_test_brier = [],[],[],[],[]
        for run in seeds:
            if n < train_df.shape[0]:
                train_sub, _ = train_test_split(
                    train_df, train_size=n, 
                    stratify=train_df[[pipeline.status_col, pipeline.batch_col]],
                    shuffle=True, 
                    random_state=run
                )
            else:
                train_sub = train_df
           
            test_sub, _  = train_test_split(
                test_df, train_size=test_size,
                stratify=test_df[[pipeline.status_col, pipeline.batch_col]],
                shuffle=True, 
                random_state=run
            )
            
            # ================== Tune hyperparameters ====================
            if is_tune:
                _ = pipeline.tune_hyperparameters(
                    train_sub, n_samples=n, n_splits=n_splits, n_trials=n_trials, 
                    trial_threshold=trial_threshold, n_jobs=n_jobs
                )
            else:
                try:
                    study_name = pipeline._study_name(n)
                    study = optuna.load_study(study_name=study_name, storage=pipeline.storage_url)
                    pipeline._best_params = study.best_params
                    print(f"Found existing Optuna study '{study_name}'. Retrieved hyperparameters: {study.best_params}")
                    
                except Exception:
                    print("No existing Optuna study found. Default hyperparameters will be used.") 
   
            # ============== Train with best hyperparameters =============
            duration,tr_brier, te_brier, tr_cind, te_cind = pipeline.train(train_sub, test_sub)
            
            model_results.append({
                "n train": int(n),
                "train time": float(duration),
                "train C": float(tr_cind),
                "test C": float(te_cind),
                "train brier": float(tr_brier),
                "test brier": float(te_brier)
            })
            run_train_cind.append(tr_cind)
            run_test_cind.append(te_cind)
            run_train_brier.append(tr_brier)
            run_test_brier.append(te_brier)
            run_time.append(duration)
            
            print(f"(runtime: {duration:.2f}s)   |\
                    (C-index)  Train: {tr_cind:.3f}, Test: {te_cind:.3f}   |\
                    (Brier)  Train: {tr_brier:.3f}, Test: {te_brier:.3f}\n")  
            
        print(
            f"(Avg. runtime: {np.mean(run_time):.2f}s)   |\
            (C-index)  Train: {round(np.mean(run_train_cind),3)}, Test: {round(np.mean(run_test_cind),3)}   |\
            (Brier)  Train: {round(np.nanmean(run_train_brier),3)}, Test: {round(np.nanmean(run_test_brier),3)} (Mean)\n"
        )                
    
    model_results = pd.DataFrame(model_results)  
    if is_save:
        pipeline.write(model_results=model_results, file_path=file_path, file_name=file_name)

    return model_results
