import os
import time 
import numpy as np

os.chdir("..")
from pycox.models.cox import CoxPH
from pycox.evaluation.eval_surv import EvalSurv

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn_pandas import DataFrameMapper

from abc import ABC, abstractmethod
import optuna
import torchtuples as tt

def parse_hyperparameters(trial, config=None):
    """
    Convert a dictionary of hyperparameter ranges into Optuna-compatible format.
    
    Args:
    - trial (optuna.trial.Trial): The Optuna trial object to sample parameters from.
    - hyperparam_dict (dict): A dictionary with hyperparameter names as keys. Values should be a dictionary specifying the type of parameter ('int', 'float', or 'categorical')
      and the range or choices.
      
    Returns:
    - dict: A dictionary of sampled hyperparameters compatible with the given ranges.
    """
    # config = self.hyperparameters if config is None else config
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


class DeepSurvPipeline(ABC):
    def __init__(self, train_df, test_df, 
                time_col='time', status_col='status',
                dataName=None, hyperparameters=None, random_state=42):
        self.train_df=train_df
        self.test_df=test_df
        self.n = train_df.shape[0]
        self.time_col = time_col
        self.status_col = status_col
        self.dataName = dataName
        self.hyperparameters = hyperparameters
        self.random_state = random_state
        self.model = None
        self.modelString = 'deepsurv-torch'
    
    def _preprocess_data(self, df):
        """Use StandardScaler() to transform input miRNA data.
        Return transformed features and labels (survival time and censor status) separately.
        """
        survival_cols = [self.time_col, self.status_col]
        covariate_cols = [col for col in df.columns if col not in survival_cols]
        standardize = [([col], StandardScaler()) for col in covariate_cols]
        x_mapper = DataFrameMapper(standardize)

        # transform features (miRNA expression)
        x = x_mapper.fit_transform(df[covariate_cols]).astype('float32')
        # prepare labels (survival data)
        get_target = lambda df: (df[self.time_col].values, df[self.status_col].values)
        y = get_target(df)
        
        return x, y
    
    def train_mlp(self, train_df, val_df, params, patience, min_delta, verbose=False, is_tuning=False):
        """Training function that works for both tuning and main training.

        Args:
            x_train (pandas DataFrame): _description_
            y_train (pandas DataFrame): _description_
            params (dictionary): neural network training parameters
                - 'batch_size': batch size
                - 'epochs': number of training epochs
                - 'num_nodes': hidden layer size
                - 'learning_rate'
                - 'dropout'
                - ...
        """
        # Set hyperparameters with default values
        num_nodes = params.get("num_nodes", [16,16])           # Default to [32, 32] if not provided
        dropout = params.get("dropout", 0.5)                    # Default dropout rate
        learning_rate = params.get("learning_rate", 1e-3)       # Default learning rate
        batch_size = params.get("batch_size", 32)               # Default batch size
        epochs = params.get("epochs", 100)                      # Default number of epochs
        batch_norm = params.get("batch_norm", True)             # Default batch normalization
        output_bias = params.get("output_bias", True)           # Default output bias
        
        # Prepare data 
        x_train, y_train = self._preprocess_data(train_df)
        val = self._preprocess_data(val_df)
        x_val, y_val = val
        input_size = x_train.shape[1]
        output_size = 1
        
        # =================== Build Neural Net ===================
        ## define network
        net = tt.practical.MLPVanilla(in_features=input_size,
                                    out_features=output_size,
                                    num_nodes=num_nodes,
                                    dropout=dropout, 
                                    batch_norm=batch_norm,
                                    output_bias=output_bias)
        # define optimizer 
        optimizer = tt.optim.Adam(weight_decay=0.01)
        model = CoxPH(net, optimizer)
        model.optimizer.set_lr(learning_rate)

        # =================== Train Model =====================
        callbacks = [tt.callbacks.EarlyStopping(patience=patience, min_delta=min_delta)]
        
        start = time.time() # Record iteration start time
        log = model.fit(x_train, y_train,
                        batch_size=batch_size,
                        epochs=epochs,
                        callbacks=callbacks, 
                        verbose=verbose,
                        val_data=val,
                        val_batch_size=batch_size
                        )
        stop = time.time() # Record time when training finished
        duration = round(stop - start, 2)
        
        # =================== Evaluation ===================
        _ = model.compute_baseline_hazards()
        val_surv = model.predict_surv_df(x_val)
        val_c_index  = EvalSurv(val_surv, y_val[0], y_val[1], censor_surv='km').concordance_td() 
        
        tr_surv = model.predict_surv_df(x_train)
        tr_c_index  = EvalSurv(tr_surv, y_train[0], y_train[1], censor_surv='km').concordance_td() 
        
        if not is_tuning:
            self.model = model
            self.optimizer = optimizer
            
        return duration, tr_c_index, val_c_index
    
    
    def _objective(self, trial, train_df, config, n_splits=5):
        """Perform K-Fold CV and return the average validation loss across folds.""" 
        
        # Define search space
        self.hyperparameters = config if config is not None else self.hyperparameters
        if self.hyperparameters is None:
            config = {
                "num_nodes": trial.suggest_categorical("num_nodes", [[32, 32], [16, 16]]),
                "dropout": trial.suggest_categorical("dropout", [0.1, 0.5, 0.9]),
                "learning_rate": trial.suggest_loguniform("learning_rate", 1e-5, 1e-2),
                "batch_size": trial.suggest_categorical("batch_size", [16, 32, 64]),
                # "epochs": trial.suggest_int("epochs", 10, 50)  
            }
        else:
            config = parse_hyperparameters(trial, config=self.hyperparameters)

        kf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=self.random_state)
        val_c_indexes = []

        for train_index, val_index in kf.split(train_df, train_df[self.status_col]):
            tr_df, val_df = train_df.iloc[train_index,:], train_df.iloc[val_index,:]

            # Train and evaluate on the current fold -- validation score
            *_, val_c_index = self.train_mlp(tr_df, val_df,
                                            params=config,
                                            patience=10, min_delta=1e-2,
                                            verbose=False,
                                            is_tuning=True)
            val_c_indexes.append(val_c_index)
        
        c_index=np.mean(val_c_indexes)
        return c_index
    
    
    def tune_hyperparameters(self, train_df, 
                            config=None,
                            storage_url=None, study_name=None,
                             n_splits=5, n_trials=20, n_jobs=2):
        """Method to perform hyperparameter tuning with K-Fold CV."""
        storage_url = f"sqlite:///{self.dataName}-hp-log.db" if storage_url is None else storage_url
        study_name = f"deepsurv_hp_{self.dataName}" if study_name is None else study_name
        
        # Create an Optuna study for hyperparameter optimization
        study = optuna.create_study(direction='maximize', storage=storage_url, study_name=study_name,
                                    load_if_exists=True)
        
        # Optimize the study using the objective function
        study.optimize(lambda trial: self._objective(trial,
                                                    train_df=train_df,
                                                    config=config,
                                                    n_splits=n_splits),
                    n_trials=n_trials, n_jobs=n_jobs)
        
        # Save the best hyperparameters
        self.best_params = study.best_params
        print(f"Best trial parameters: {self.best_params}")
        return study
        
        
    def train_with_best_params(self, val_size=0.2, 
                            params=None, 
                            patience=25, min_delta=1e-2,
                            verbose=True, print_scores=False):
        """
        Args:
            subset (_type_): _description_
            batch_sizes (_type_): _description_
            val_size (float, optional): _description_. Defaults to 0.2.
            kwargs (_type_, optional): _description_. Defaults to None.
            time_col (str, optional): _description_. Defaults to 'time'.
            status_col (str, optional): _description_. Defaults to 'status'.
            verbose (bool, optional): _description_. Defaults to True.
            print_scores (bool, optional): _description_. Defaults to False.

        Returns:
            _type_: _description_
        """
        if params is None and self.best_params is None:
            raise ValueError("No hyperparameters found. Run tune_hyperparameters() first or define them manually.")
        elif params is None:
            params = self.best_params
        
        # ===================== Prepare Data =======================
        tr_df, val_df = train_test_split(self.train_df, 
                                        test_size=val_size,
                                        shuffle=True, random_state=self.random_state,
                                        stratify=self.train_df[self.status_col])
        
        # ===================== Train Model =====================
        duration, train_c, _ = self.train_mlp(tr_df, val_df, params=params, 
                                            patience=patience, min_delta=min_delta,
                                            verbose=verbose, is_tuning=False)

        # ===================== Evaluation =====================
        x_test, y_test = self._preprocess_data(self.test_df)
        te_surv = self.model.predict_surv_df(x_test)
        test_c  = EvalSurv(te_surv, y_test[0], y_test[1], censor_surv='km').concordance_td()
        
        if print_scores:
            print(f"N={self.n} Training time ({duration}s): Train C-Index: {round(train_c,3)} | Test C-index: {round(test_c,3)}")
            
        return duration, train_c, test_c
