from __future__ import annotations
from typing import Dict, Any
import time
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn_pandas import DataFrameMapper

import torch
import torchtuples as tt
from pss.pycox.models.cox import CoxPH, CoxPHStratified, StratifiedDataset
from pss.pycox.evaluation.eval_surv import EvalSurv
from pss.run_models.common import TunablePipelineBase, ResultsWriterMixin, _model_string, _default_hp_space


class DeepSurvPipeline(TunablePipelineBase, ResultsWriterMixin):
    """
    DeepSurv suvival risk prediction models with built-in hyperparameter tuning (optuna). 

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
        train_df: pd.DataFrame,
        test_df: pd.DataFrame,
        batchNormType: str | None = None,
        dataName: str | None = None,
        is_stratified: bool = False,
        time_col: str = "time",
        status_col: str = "status",
        batch_col: str = "batch.id",
        hyperparameters: Dict[str, Any] | None = None,
        storage_url: str = "sqlite:///deepsurv-torch-hp-log.db"
    ):
        modelString = _model_string("deepsurv-torch", is_stratified)
        hyperparameters = _default_hp_space("deepsurv-torch") if hyperparameters is None else hyperparameters
        super().__init__(
            batchNormType=batchNormType,
            dataName=dataName, 
            modelString=modelString,
            hyperparameters=hyperparameters,
            storage_url=storage_url
            )
        self.train_df = train_df
        self.test_df = test_df
        self.time_col = time_col
        self.status_col = status_col
        self.batch_col = batch_col
        self.is_stratified = is_stratified
        self.patience = 30
        self.min_delta = 1e-3
        self._best_params = {}

    # --------------------------------------------------
    #  Optuna objective (k-fold mean validation C-index) 
    # --------------------------------------------------
    def _objective(self, df, n_splits, params, fixed_params=None, **_):
  
        full_params = {**fixed_params, **params} if fixed_params else params
        
        kf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
        # stratify by "status_batch" for balanced assignment
        stratify_labels = df[[self.status_col, self.batch_col]].astype(str).agg("_".join, axis=1).values
        
        val_c_indexes = []
        for train_index, val_index in kf.split(df, stratify_labels):
            tr_df, val_df = df.iloc[train_index, :], df.iloc[val_index, :]

            *_, val_c = self._trainer(tr_df, val_df, 
                                     params=full_params,
                                     verbose=False)
            val_c_indexes.append(float(val_c))
            
        return float(np.mean(val_c_indexes))

    def _trainer(self, 
        train_df: pd.DataFrame,
        val_df: pd.DataFrame, 
        params: Dict[str, Any] | None = None,
        patience: int = None, 
        min_delta: float = None, 
        verbose: bool = False
    ):
        """MLP trining function that works for both hyperparameter tuning and model training.
        
        Args:
            train_df (pandas.DataFrame): training data.
            val_df (pandas.DataFrame): testing/validatoin data.
            params (dict): neural network training parameters.
                - batch_size (int): default = 128
                    Number of samples used for calculating loss and updating weights.  
                - epochs (int): default = 500 
                    Number of times the entire dataset is being processed.
                - num_nodes (list): default = [32,16]
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
                    Normalize weights in each mini batch.
                - output_bias (bool): default = True
            patience (int): The number of epochs allowed for no improvement in loss before stopping the network training.
            min_delta (float): Minimum improvement for early stopping.
            verbose (bool): Whether to print training logs, default = False.
                        
        Returns:
            duration: float 
                Training time in seconds.
            tr_c_index, val_c_index: float
                Harrell's concordance index scores.
            tr_brier, val_brier: float 
                Integrated Brier scores calculated at 50 time points
                selected evenly throughout min and max survival time in the train and test data.
        """
        # ==================== Prepare data ====================
        def _preprocess_data(self, df, mapper=None, fit_scaler=True):
            """ 
            Applies StandardScaler() to the input data features and re-format outcome to tuples. 
            """
            survival_cols = [self.time_col, self.status_col]
            covariate_cols = [col for col in df.columns if col not in survival_cols]
            
            # Transform features (miRNA expression)
            if fit_scaler or mapper is None:
                standardize = [([col], StandardScaler()) for col in covariate_cols]
                mapper = DataFrameMapper(standardize)
                
                x = mapper.fit_transform(df[covariate_cols]).astype('float32')
            else:
                x = mapper.transform(df[covariate_cols]).astype('float32')
                
            # Prepare labels (survival data)
            y = (df[self.time_col].values, df[self.status_col].values)
            
            return x, y, mapper
       
        train_df = train_df.drop(columns=self.batch_col)
        val_df = val_df.drop(columns=self.batch_col)
        x_train, y_train, mapper = _preprocess_data(train_df, fit_scaler=True)
        x_val, y_val, _ = _preprocess_data(val_df, mapper=mapper, fit_scaler=False)
        
        batch_ids_train = train_df[self.batch_col].to_numpy().reshape(-1)
        batch_ids_val = val_df[[self.batch_col]].to_numpy().reshape(-1)
        
        # ==================== GPU integration ====================
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        
        # Convert numpy arrays to torch tensors and move to GPU
        x_train         = torch.from_numpy(x_train).to(device)
        x_val           = torch.from_numpy(x_val).to(device)
        durations_train = torch.from_numpy(y_train[0]).float().to(device)
        durations_val   = torch.from_numpy(y_val[0]).float().to(device)
        events_train    = torch.from_numpy(y_train[1]).float().to(device)
        events_val      = torch.from_numpy(y_val[1]).float().to(device)
        batch_ids_train = torch.from_numpy(batch_ids_train).long().to(device)
        batch_ids_val   = torch.from_numpy(batch_ids_val).long().to(device)
        
        # ==================== Set hyperparameters ====================
        if params is None and self.best_params is None:
            print("No hyperparameters provided. Training with defaults.")
            params = {}
        elif params is None:
            params = self.best_params
            
        input_size = x_train.shape[1]
        output_size = 1
        num_nodes = params.get("num_nodes", [32,16])            # Default num of layers & nodes
        dropout = params.get("dropout", 0.1)                    # Default dropout rate
        learning_rate = params.get("learning_rate", 1e-3)       # Default learning rate
        batch_size = params.get("batch_size", 128)              # Default batch size
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
        
        # ==================== Build Neural Net ====================
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
        
        # Early stopping settings
        patience = self.patience if patience is None else patience
        min_delta = self.min_delta if min_delta is None else min_delta
        callbacks = [tt.callbacks.EarlyStopping(patience=patience, min_delta=min_delta)]
        
        # ==================== Train Model ====================
        start = time.time() # Record iteration start time
        if self.is_stratified:
            train_dataset = StratifiedDataset(x_train, durations_train, events_train, batch_ids_train)
            val_dataset = StratifiedDataset(x_val, durations_val, events_val, batch_ids_val)
            # train_dataset = torch.utils.data.TensorDataset(x_train, durations_train, events_train, batch_ids_train)
            # val_dataset   = torch.utils.data.TensorDataset(x_val, durations_val, events_val, batch_ids_val)
            train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
            val_loader   = torch.utils.data.DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
            
            model = CoxPHStratified(net, optimizer)
            _ = model.fit_dataloader(
                train_loader,
                epochs=epochs,
                callbacks=callbacks,
                verbose=verbose,
                val_dataloader=val_loader)
        else:
            y_train = (durations_train, events_train)
            y_val = (durations_val, events_val)
            
            model = CoxPH(net, optimizer)
            _ = model.fit(
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
        # Convert torch tensors back to numpy objects for evaluation
        durations_train_np = durations_train.detach().cpu().numpy()
        durations_val_np   = durations_val.detach().cpu().numpy()
        events_train_np    = events_train.detach().cpu().numpy()
        events_val_np      = events_val.detach().cpu().numpy()

        # Compute baseline hazards (single stratum or per-batch)
        if self.is_stratified:
            _ = model.compute_baseline_hazards(input=x_train, target=(durations_train, events_train), batch_ids=batch_ids_val)   
            tr_surv  = model.predict_surv_df(x_train, batch_ids=batch_ids_train)
            val_surv = model.predict_surv_df(x_val, batch_ids=batch_ids_val)
        else:
            _ = model.compute_baseline_hazards(input=x_train, target=(durations_train, events_train)) 
            tr_surv  = model.predict_surv_df(x_train)
            val_surv = model.predict_surv_df(x_val)  
            
        # Initialize EvalSurv objects 
        tr_ev = EvalSurv(tr_surv, durations_train_np, events_train_np, censor_surv='km')
        val_ev = EvalSurv(val_surv, durations_val_np, events_val_np, censor_surv='km')
        
        # Concordance index (C-index) --------------
        if self.is_stratified:
            tr_c_index  = tr_ev.stratified_concordance_td(batch_indices=batch_ids_train) 
            val_c_index = val_ev.stratified_concordance_td(batch_indices=batch_ids_val) 
        else:
            tr_c_index, _  = tr_ev.concordance_td() 
            val_c_index, _ = val_ev.concordance_td() 
        
        # Integrated Brier score (IBS) -------------
        min_surv = np.ceil(max(np.min(durations_train_np), np.min(durations_val_np)))
        max_surv = np.floor(min(np.max(durations_train_np), np.max(durations_val_np)))
        times = np.linspace(min_surv, max_surv, 20)
        
        if self.is_stratified:
            tr_brier  = tr_ev.stratified_integrated_brier_score(time_grid=times, batch_indices=batch_ids_train) 
            val_brier = val_ev.stratified_integrated_brier_score(time_grid=times, batch_indices=batch_ids_val)
        else:
            tr_brier  = tr_ev.integrated_brier_score(time_grid=times) 
            val_brier = val_ev.integrated_brier_score(time_grid=times) 

        return float(duration), float(tr_brier), float(val_brier), float(tr_c_index), float(val_c_index)
    
    # -------------------------------------------------------------
    # Train/Eval once using MLP _trainer() with self._best_params()
    # -------------------------------------------------------------
    def train(
        self,
        params: Dict[str, Any] | None = None,
        *,
        train_df: pd.DataFrame | None = None,
        test_df: pd.DataFrame | None = None) -> Dict[str, float]:
        
        tr_df = self.train_df if train_df is None else train_df
        te_df = self.test_df if test_df is None else test_df

        params = params or self.get_best_params()

        return self._trainer(tr_df, te_df, params=params, verbose=False)
