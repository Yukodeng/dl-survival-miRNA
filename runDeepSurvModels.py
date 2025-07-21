import os
import time
from datetime import datetime
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn_pandas import DataFrameMapper

import torch
import torchtuples as tt
from utils import *
from pycox.models.cox import CoxPH, CoxPHStratified, StratifiedDataset
from pycox.evaluation.eval_surv import EvalSurv
import optuna

np.random.seed(123)
torch.manual_seed(123)

def parse_hyperparameters(trial, config=None):
    """
    Convert a dictionary of hyperparameter ranges into Optuna-compatible format.
    
    Args:
    - trial (optuna.trial.Trial): The Optuna trial object to sample parameters from.
    - config (dict): A dictionary with hyperparameter names as keys. Values should be a dictionary specifying the type of parameter ('int', 'float', or 'categorical')
    and the range or choices.
    
    Returns:
    - dict: A dictionary of sampled hyperparameters compatible with the given ranges.
    """
    params = {}
    for param_name, param_info in config.items():
        param_type = param_info.get("type")
        if param_type == "int":
            params[param_name] = trial.suggest_int(param_name, param_info["low"], param_info["high"])
        elif param_type == "float":
            params[param_name] = trial.suggest_float(param_name, param_info["low"], param_info["high"], log=param_info.get("log", False))
        elif param_type == "categorical":
            params[param_name] = trial.suggest_categorical(param_name, param_info["choices"])
        else:
            raise ValueError(f"Unsupported hyperparameter type: {param_type}")
    return params


def filter_hyperparams_by_size(config, n_samples, exclusion_dict):
    exclude_keys = exclusion_dict.get(n_samples, [])
    return {k: v for k, v in config.items() if k not in exclude_keys}


class DeepSurvPipeline():
    """_summary_
    """
    def __init__(
        self,
        train_df, 
        test_df, 
        test_size=1000,
        time_col='time', status_col='status', batch_col="batch.id",
        batchNormType=None,
        dataName=None, 
        hyperparameters=None,
        is_stratified=False, 
        storage_url="sqlite:///deepsurv-torch-hp-log.db",
        random_state=42
    ):
        self.train_df = train_df
        self.test_df = test_df
        self.n = train_df.shape[0]
        self.test_size = test_size
        self.time_col = time_col
        self.status_col = status_col
        self.batch_col = batch_col
        self.batchNormType = batchNormType
        self.dataName = dataName
        self.is_stratified = is_stratified
        self.hyperparameters = hyperparameters
        self.storage_url = storage_url
        self.random_state = random_state
        self.best_params = None
        self.model = None
        self.modelString = "%sdeepsurv-torch" % (['','stratified-'][self.is_stratified])
        self.patience = 30
        self.min_delta = 1e-3
    
    
    def _preprocess_data(self, df, mapper=None, fit_scaler=True):
        """
        Applies StandardScaler() to transform the input miRNAseq data frame.
        
        Parameters:
        df (pd.DataFrame): Input data.
        mapper (DataFrameMapper or None): If provided, uses this mapper to transform data.
        fit_scaler (bool): Whether to fit a new mapper on the current data.
        
        Returns:
        x (np.ndarray): Scaled features.
        y (tuple): (time, event) tuple for survival modeling.
        mapper (DataFrameMapper): The fitted mapper used for transformation.
        """
        survival_cols = [self.time_col, self.status_col]
        covariate_cols = [col for col in df.columns if col not in survival_cols]
        
        if fit_scaler or mapper is None:
            standardize = [([col], StandardScaler()) for col in covariate_cols]
            mapper = DataFrameMapper(standardize)
            # Transform features (miRNA expression)
            x = mapper.fit_transform(df[covariate_cols]).astype('float32')
        else:
            x = mapper.transform(df[covariate_cols]).astype('float32')
        
        # Prepare labels (survival data)
        y = (df[self.time_col].values, df[self.status_col].values)
        
        return x, y, mapper
    
    def write(self, model_results, out_dir=None, fileName=None):
        
        out_dir=os.path.join('models', self.batchNormType, self.dataName, self.modelString) if out_dir is None else out_dir
        os.makedirs(out_dir, exist_ok=True)
        
        # save results as txt or csv file
        today = datetime.now().strftime("%Y%m%d")
        fileName = f'model_results_{today}.csv' if fileName is None else fileName
        
        if 'txt' in fileName:
            model_results.to_csv(os.path.join(out_dir,fileName), sep='\t')
        elif 'csv' in fileName:
            model_results.to_csv(os.path.join(out_dir,fileName), index=False)
        else:
            fileName = fileName + '.csv'
            model_results.to_csv(os.path.join(out_dir,fileName), index=False)
            # print('Please specify a file name with either a txt or csv extension.')
    
    def train_mlp(self, 
                train_df=None,
                val_df=None, 
                params=None,
                patience=None, 
                min_delta=None, 
                verbose=False, 
                is_tuning=False,
                calculate_brier=True):
        """Training function that works for both tuning and main training.
        
        Args:
            x_train (pandas.DataFrame): _description_
            y_train (pandas.DataFrame): _description_
            params (dict): neural network training parameters
                - batch_size (int): default = 32
                    Number of samples used for calculating loss and updating weights.  
                - epochs (int): default = 500 
                    Number of times the entire dataset is being processed.
                - num_nodes (list): default = [32,32]
                    Number of hidden layers and nodes in each layer. 
                - learning_rate (float): default = 0.001
                    The rate at which weights are updated.
                - weight_decay (float): default = 0.0001
                    Regularization parameter applied to the loss function. 
                - dropout (float): range 0-1, default = 0.1 
                    Drop out rate.
                - activation (str): {'ReLU', 'LeakyReLU', 'SELU'}, default = "ReLU" 
                    Activation function applied the hidden layers.
                - batch_norm (bool): default = True
                    ...
                - output_bias (bool): default = True
                    ...
            patience (int): The number of epochs allowed for no improvement in loss before stopping the network training.
            min_delta (float): Minimum improvement for early stopping.
            verbose (bool): Whether to print training logs.
            is_tuning (bool): Whether to skip saving model.
            calculate_brier (bool): Whether to compute Brier score.
            
        Returns:
            duration: float 
                Training time in seconds.
            tr_c_index, val_c_index: float
                Harrell's concordance index scores.
            tr_brier, val_brier: float 
                Integrated Brier scores calculated at 50 time points
                selected evenly throughout min and max survival time in the train and test data.
        """
        # ================== Prepare data ====================
        train_df = train_df if train_df is not None else self.train_df
        val_df  = val_df if val_df is not None else self.test_df
        
        batch_ids_train = train_df[self.batch_col].to_numpy().reshape(-1)
        batch_ids_val = val_df[[self.batch_col]].to_numpy().reshape(-1)
        
        train_df = train_df.drop(columns=self.batch_col)
        val_df = val_df.drop(columns=self.batch_col)
        x_train, y_train, mapper = self._preprocess_data(train_df, fit_scaler=True)
        x_val, y_val, _ = self._preprocess_data(val_df, mapper=mapper, fit_scaler=False)
        
        # ================== GPU integration =================
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        
        # Convert numpy arrays to torch tensors and move to GPU
        x_train = torch.from_numpy(x_train).to(device)
        x_val = torch.from_numpy(x_val).to(device)
        
        durations_train = torch.from_numpy(y_train[0]).float().to(device)
        events_train = torch.from_numpy(y_train[1]).float().to(device)
        durations_val = torch.from_numpy(y_val[0]).float().to(device)
        events_val = torch.from_numpy(y_val[1]).float().to(device)
        batch_ids_train = torch.from_numpy(batch_ids_train).long().to(device)
        batch_ids_val = torch.from_numpy(batch_ids_val).long().to(device)
        
        y_train = (durations_train, events_train)
        y_val = (durations_val, events_val)
        
        # =============== Set hyperparameters ================
        if params is None and self.best_params is None:
            print("No hyperparameters provided. Training with defaults.")
            params = {}
        elif params is None:
            params = self.best_params
            
        input_size = x_train.shape[1]
        output_size = 1
        num_nodes = params.get("num_nodes", [32,32])            # Default num of layers & nodes
        dropout = params.get("dropout", 0.1)                    # Default dropout rate
        learning_rate = params.get("learning_rate", 1e-4)       # Default learning rate
        batch_size = params.get("batch_size", 64)               # Default batch size
        epochs = params.get("epochs", 500)                      # Default number of epochs
        batch_norm = params.get("batch_norm", True)             # Default batch normalization
        output_bias = params.get("output_bias", True)           # Default output bias
        weight_decay = params.get("weight_decay", 1e-4)         # Default weight decay
        activation_map = {                                      
            "ReLU": torch.nn.ReLU,
            "LeakyReLU": torch.nn.LeakyReLU,
            "SELU": torch.nn.SELU
        } 
        activation = activation_map.get(params.get("activation", "ReLU")) # Activation function
        
        # =================== Build Neural Net ===================
        # Define network
        net = tt.practical.MLPVanilla(
            in_features=input_size,
            out_features=output_size,
            num_nodes=num_nodes,
            dropout=dropout, 
            batch_norm=batch_norm,
            activation=activation,
            output_bias=output_bias
        )
        net = net.to(device) # send to GPU
        
        # Define optimizer 
        optimizer = tt.optim.Adam(weight_decay=weight_decay, lr=learning_rate)
        
        # Get default early stopping settings if not defined 
        patience = self.patience if patience is None else patience
        min_delta = self.min_delta if min_delta is None else min_delta
        callbacks = [tt.callbacks.EarlyStopping(patience=patience, min_delta=min_delta)]
        
        # =================== Train Model ====================
        start = time.time() # Record iteration start time
        if self.is_stratified:
            train_dataset = torch.utils.data.TensorDataset(x_train, durations_train, events_train, batch_ids_train)
            val_dataset = torch.utils.data.TensorDataset(x_val, durations_val, events_val, batch_ids_val)
            train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
            val_loader = torch.utils.data.DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
            
            model = CoxPHStratified(net, optimizer)
            log = model.fit_dataloader(
                train_loader,
                epochs=epochs,
                callbacks=callbacks,
                verbose=verbose,
                val_dataloader=val_loader
            )
        else:
            model = CoxPH(net, optimizer)
            log = model.fit(
                x_train, y_train,
                batch_size=batch_size,
                epochs=epochs,
                callbacks=callbacks, 
                verbose=verbose,
                val_data=(x_val, y_val),
                val_batch_size=batch_size
            )
        stop = time.time() # Record time when training finished
        duration = round(stop - start, 2)
        
        # ==================== Evaluation ====================
        _ = model.compute_baseline_hazards(input=x_train, target=(durations_train, events_train))   
        
        # Convert torch tensors back to numpy objects for evaluation
        x_train_np = x_train.detach().cpu().numpy()
        x_val_np = x_val.detach().cpu().numpy()
        
        durations_train_np = durations_train.detach().cpu().numpy()
        durations_val_np   = durations_val.detach().cpu().numpy()
        events_train_np    = events_train.detach().cpu().numpy()
        events_val_np      = events_val.detach().cpu().numpy()
        
        # Initialize EvalSurv objects 
        tr_surv  = model.predict_surv_df(x_train_np)
        val_surv = model.predict_surv_df(x_val_np)
        tr_ev = EvalSurv(tr_surv, durations_train_np, events_train_np, censor_surv='km')
        val_ev = EvalSurv(val_surv, durations_val_np, events_val_np, censor_surv='km')
        
        # Concordance index ----------------
        if self.is_stratified:
            tr_c_index  = tr_ev.stratified_concordance_td(batch_indices=batch_ids_train) 
            val_c_index = val_ev.stratified_concordance_td(batch_indices=batch_ids_val) 
        else:
            tr_c_index  = tr_ev.concordance_td() 
            val_c_index = val_ev.concordance_td() 
        
        # Integrated Brier score -----------
        min_surv = np.ceil(max(np.min(durations_train_np), np.min(durations_val_np)))
        max_surv = np.floor(min(np.max(durations_train_np), np.max(durations_val_np)))
        times = np.linspace(min_surv, max_surv, 20)
        
        if calculate_brier:
            if self.is_stratified:
                tr_brier  = tr_ev.stratified_integrated_brier_score(time_grid=times, batch_indices=batch_ids_train) 
                val_brier = val_ev.stratified_integrated_brier_score(time_grid=times, batch_indices=batch_ids_val)
            else:
                tr_brier  = tr_ev.integrated_brier_score(time_grid=times) 
                val_brier = val_ev.integrated_brier_score(time_grid=times) 
        else:
            tr_brier = val_brier = np.nan
        
        if not is_tuning:
            self.model = model
            self.optimizer = optimizer
            
        # print(
        #     f"""
        #     N={self.n} Training time: ({duration}s)
        #         C index:  Train ({round(tr_c_index, 3)})  |  Test ({round(val_c_index, 3)})
        #         Brier score:  Train ({round(tr_brier, 3)})  |  Test ({round(val_brier, 3)})
        #     """ 
        # )  
        return duration, tr_brier, val_brier, tr_c_index, val_c_index
    
    
    def _objective(self, trial, train_df, config, n_splits, fixed_params=None):
        """Perform K-Fold CV and return the average validation loss across folds.""" 
        
        tunable_params = parse_hyperparameters(trial, config=config)
        
        full_params = {**fixed_params, **tunable_params} if fixed_params else tunable_params
        
        kf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=self.random_state)
        stratify_labels = train_df[self.status_col].astype(str) #+ "_" + train_df[self.batch_col].astype(str)
        
        val_c_indexes = []
        for train_index, val_index in kf.split(train_df, stratify_labels):
            tr_df, val_df = train_df.iloc[train_index,:], train_df.iloc[val_index,:]

            # Train and evaluate on the current fold -- validation score
            *_, val_c_index = self.train_mlp(tr_df, val_df,
                                            params=full_params,
                                            patience=self.patience, 
                                            min_delta=self.min_delta,
                                            verbose=False,
                                            is_tuning=True)
            val_c_indexes.append(val_c_index)
        
        return np.mean(val_c_indexes)
    
    
    def tune_hyperparameters(self, 
                            train_df, 
                            n_samples=None,
                            config=None,
                            storage_url=None, 
                            study_name=None,
                            n_splits=10,
                            trial_threshold=30,
                            n_trials=30,
                            n_jobs=1, 
                            timeout=2400):
        """
        Perform hyperparameter tuning within the K-fold CV framework using the Optuna library.
        Provides option to store tuning results and retrieve from existing studies to avoid duplicate optimization processes. 
        Studies are designed to maximize the average validation Harrell's concordance index (C-index) scores.  
        
        Args:
            train_df: Training dataset.
            n_samples: default = None
            config: Search grid for hyperparameter tuning.
            storage_url: URL for Optuna storage.
            study_name: Name of the study.
            n_splits: Number of splits for cross-validation.
            n_trials: Number of trials to run for tuning.
            n_jobs: Number of parallel jobs for tuning.
            trial_threshold: Minimum number of trials required to skip tuning.
        """
        n_samples = n_samples if n_samples is not None else train_df.shape[0]
        study_name = f"{self.batchNormType}-{self.dataName}-{n_samples}" if study_name is None else study_name
        storage_url = self.storage_url if storage_url is None else storage_url
        
        # Create an Optuna study for hyperparameter optimization
        study = optuna.create_study(direction='maximize',
                                    storage=storage_url, 
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
        
        # @TODO: selectively drop hyperparameters
        # excluded_params_by_size = {
        #     100: ["batch_size", "dropout", "activation"],
        #     500: ["batch_size"]
        # } --> add to class attribute
        # filtered_config = filter_hyperparams_by_size(new_config, n_samples, self.exclusion_dict)
        
        # Optimize the study using the objective function
        print(f"Starting hyperparam tuning: {list(new_config.keys())}")
        
        study.optimize(
            lambda trial: self._objective(trial,
                                        train_df=train_df,
                                        config=new_config,
                                        fixed_params=fixed_params,
                                        n_splits=n_splits),
            n_trials=n_trials, 
            n_jobs=n_jobs,
            timeout=timeout
        )
        
        # Ensure not all trials failed
        successful_trials = [t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE]
        if len(successful_trials) == 0:
            raise ValueError("All trials failed. Consider changing the hyperparameter grid and retrying.")
        
        # Save the best hyperparameters
        self.best_params = {**fixed_params, **study.best_params}
        print(f"Found best hyperparameters: {self.best_params}")
        
        return study
    
    def train_over_subsets(self,
                           subset_sizes=[100, 500, 1000, 2000, 5000, 10000],
                           runs_per_size=[20,20,20,20,20,20],
                           splits_per_size=[3,5,5,10,10,10],
                           trials_per_size=[30,30,30,30,30,30],
                           trial_threshold=30,
                           n_jobs=1,
                           is_tune=False,
                           is_save=False,
                           fileName="model_results.csv"):
        
        n_train, model_time, model_train_cind, model_test_cind, model_train_brier, model_test_brier = [],[],[],[],[],[]
        
        for (n, n_run, n_splits, n_trials) in zip(subset_sizes, runs_per_size, splits_per_size, trials_per_size): 
            
            run_time, run_train_cind, run_test_cind, run_train_brier, run_test_brier = [],[],[],[],[]
            
            print(f"Running for N={n}...")
            
            # Apply early stopping based on training sample size
            if hasattr(self, "early_stop_per_size"):
                stop_params = self.early_stop_per_size.get(str(n), {})
                self.patience = stop_params.get("patience", 30)
                self.min_delta = stop_params.get("min_delta", 1e-3)
            
            for run in range(n_run):
                if n < self.train_df.shape[0]:
                    train_sub, _ = train_test_split(self.train_df,
                                                    train_size=n, # train_size= n/self.train_df.shape[0], 
                                                    shuffle=True,
                                                    random_state=run,
                                                    stratify=self.train_df[[self.status_col, self.batch_col]])
                else:
                    train_sub = self.train_df
                
                ## update 6/3/2025: vary test set across iterations
                if self.test_size < self.test_df.shape[0]:
                    test_sub, _ = train_test_split(self.test_df,
                                                train_size=self.test_size, # train_size= n/self.train_df.shape[0], 
                                                shuffle=True,
                                                random_state=run,
                                                stratify=self.test_df[[self.status_col, self.batch_col]])
                else:
                    test_sub = self.test_df
                
                if run == 0:
                    if is_tune:     
                        self.tune_hyperparameters(
                            train_sub, n_splits=n_splits, n_samples=n, n_trials=n_trials, n_jobs=n_jobs, trial_threshold=trial_threshold
                        )
                        best_params = self.best_params
                    
                    else:
                        try:
                            study_name = f"{self.batchNormType}-{self.dataName}-{n}"
                            study = optuna.load_study(study_name=study_name, storage=self.storage_url)
                            best_params = study.best_params
                            print(f"Retrieved best hyperparameters from {study_name}: {best_params}")
                        except Exception as e:
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
                    
                    
                duration, train_brier, test_brier, train_c, test_c = self.train_mlp(
                    train_df=train_sub, val_df=test_sub, params=best_params, verbose=False, calculate_brier=True
                )
                
                n_train.append(n)
                model_time.append(duration)
                model_train_cind.append(train_c)
                model_test_cind.append(test_c)
                model_train_brier.append(train_brier)
                model_test_brier.append(test_brier)
                
                run_train_cind.append(train_c)
                run_test_cind.append(test_c)
                run_train_brier.append(train_brier)
                run_test_brier.append(test_brier)
                run_time.append(duration)
            
            print(
                f"(Avg. runtime: {np.mean(run_time):.2f}s)   |\
                (C-index)  Train: {round(np.mean(run_train_cind),3)}, Test: {round(np.mean(run_test_cind),3)}   |\
                (Brier)  Train: {round(np.nanmean(run_train_brier),3)}, Test: {round(np.nanmean(run_test_brier),3)} (Mean)\n"
            )                
            
        model_results = pd.DataFrame({
            'n train': n_train, 
            'train time': model_time,
            'train C': model_train_cind, 
            'test C': model_test_cind,
            'train brier': model_train_brier,
            'test brier': model_test_brier
        })
        
        if is_save:
            self.write(model_results=model_results, fileName=fileName)
            
        return model_results