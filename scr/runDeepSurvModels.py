import os 
from datetime import date
import time
from abc import ABC, abstractmethod

from .utils import *


class DLSurvModels(ABC):
    def __init__(self, subset, epochs, train_df, test_df, dataName, hyperparams=None):
        """
        This class implements and trains DL/ML survival models with
        a variety of methods available as choices.

        Args:
            subset (_type_): _description_
            epochs (_type_): _description_
        """
        self.subset      = subset
        self.epochs      = epochs
        self.train_df    = train_df
        self.test_df     = test_df
        self.dataName    = dataName
        self.hyperparams = hyperparams
        self.model       = None
        self.modelString = None
        self.date        = date.today().strftime("%m%d%y")
    
    @abstractmethod
    def build_model(self):
        pass
        
    def tune_hyperparameters(self):
        pass
    
    def run_deepsurvModel(self):
        pass
    
    def write(self, model_results, out_dir=None, fileName='model.results.txt'):
        out_dir=os.path.join('models', self.dataName, self.modelString) if out_dir is None else out_dir
        os.makedirs(out_dir, exist_ok=True)
        
        # save results as txt or csv file
        if 'txt' in fileName:
            model_results.to_csv(os.path.join(out_dir,fileName), sep='\t')
        elif 'csv' in fileName:
            model_results.to_csv(os.path.join(out_dir,fileName), index=False)
        else:
            print('Please specify a file name with either a txt or csv extension.')  



class DeepSurvModel(DLSurvModels):
    
    def build_model(self):
        self.modelString = 'deepsurv'
        
        if self.hyperparams is None:
            self.hyperparams = {}
            
            
    def run_deepsurvModel(self, kwargs, log_dir=None, print_scores=False):
        os.chdir("DeepSurv/deepsurv")
        import deepsurv
        from deepsurv_logger import DeepSurvLogger, TensorboardLogger
        import lasagne
        import utils 
        import viz
        os.chdir('../..')
        
        
        # prepare data
        train_data = dataframe_to_deepsurv_ds(self.train_df)
        test_data = dataframe_to_deepsurv_ds(self.test_df)
        
        # Initialize results dictionaries for various training sizes
        metrics_dict = {}
        model_time, train_scores, test_scores = [],[],[]
        
        for n, epoch in zip(self.subset, self.epochs):
            #=============== Prepare data ===================
            train_data_sub = {}
            if n < len(val):
                for key, val in train_data.items():
                    val_sub, _ = train_test_split(val,
                                            train_size=n/len(val), 
                                            shuffle=True, random_state=42,
                                            stratify=train_data['e'])
                    train_data_sub[key] = val_sub
            else:
                train_data_sub = train_data
                
            #=============== Defining the model ===============
            if kwargs is None:
                kwargs = {'learning_rate': 1e-5 , 'n_in': n}
            else:
                kwargs['n_in'] = n
            ds = deepsurv.DeepSurv(**kwargs)
            
            # if log_dir is not None:
            #     experiment_name = f'deepsurv-{n}-samples'
            #     logger = TensorboardLogger(experiment_name, logdir=log_dir)
            # else:
            #     logger = None
            update_fn=lasagne.updates.nesterov_momentum # The type of optimizer to use. \
                                                        # Check out http://lasagne.readthedocs.io/en/latest/modules/updates.html \
                                                        # for other optimizers to use
            
            #=============== Training the model ===================
            start = time.time() # Record iteration start time
            metrics = ds.train(train_data_sub, n_epochs=epoch, logger=None, update_fn=update_fn)
            metrics_dict[str(n)] = metrics
            stop = time.time() # Record time when training finished
            duration = round(stop - start, 2)
            
            #=============== Get the final metrics ===============
            train_score = ds.get_concordance_index(train_data_sub['x'], train_data_sub['t'], train_data_sub['e'])
            test_score = ds.get_concordance_index(test_data['x'], test_data['t'], test_data['e'])
            train_scores.append(train_score)
            test_scores.append(test_score)
            model_time.append(duration)
            
            if print_scores:
                print(f'{n} training samples:\tTrain C-Index: {train_score}  | Test C-Index: {test_score}')
            
        model_results = pd.DataFrame({
            'n train': self.subset, 
            'train time': model_time,
            'train score': train_scores, 
            'test score': test_scores}
        )
        return metrics_dict, model_results



class DeepSurvKModel(DLSurvModels):
    def build_model(self):
        self.modelString = 'deepsurvk'
        if self.hyperparams is None:
            self.hyperparams = {}
            
    def run_deepsurvModel(self, 
                        params,
                        batch_size=None,
                        nan_stopping=True,
                        early_stopping=True,
                        save_model=True,
                        tensorboard=True,
                        patience = 20, min_delta = 1, #early stopping params
                        print_scores=False,):
        import deepsurvk
        import tensorflow as tf
        from keras.callbacks import TensorBoard, ModelCheckpoint, EarlyStopping, TerminateOnNaN
        
        X_train, Y_train, E_train, = dataframe_to_deepsurv_ds(self.train_df, is_deepsurvk=True)
        X_test, Y_test, E_test = dataframe_to_deepsurv_ds(self.test_df, is_deepsurvk=True)
        
        history_dict, dsk_dict = {},{}
        model_time, train_scores, test_scores = [],[],[]
        
        for n, epoch in zip(self.subset, self.epochs):
            ################ prepare training subset ###############
            if n < X_train.shape[0]:
                x_train,_, y_train,_, e_train,_ = train_test_split(
                    X_train, Y_train, E_train, 
                    train_size=n/X_train.shape[0], 
                    shuffle=True, random_state=42, stratify=E_train
                )
            else:
                x_train = X_train
                y_train = Y_train
                e_train = E_train
            n_patients_train = x_train.shape[0]
            n_features = x_train.shape[1]
            
            ################# Create model #################
            dsk = deepsurvk.DeepSurvK(n_features=n_features, E=e_train, **params) 
            
            # callbacks = deepsurvk.common_callbacks()
            callbacks = []
            if nan_stopping:
                callbacks.append(TerminateOnNaN())
            if early_stopping:
                callbacks.append(
                    EarlyStopping(monitor='loss', patience=patience, min_delta=min_delta))
            if save_model:
                path_model = os.path.join('models',self.dataName,self.modelString,f'{n}-{self.date}.h5')
                callbacks.append(
                    ModelCheckpoint(path_model, monitor='loss', save_best_only=True, mode='min'))
            if tensorboard:
                tb_log_dir = os.path.join('models', self.modelString, 'tb')
                tb_cb = TensorBoard(log_dir=tb_log_dir, histogram_freq=1, write_grads=True)
                callbacks.append(tb_cb)
            
            
            ################ Train model #################
            start = time.time() # Record iteration start time
            history_dict[str(n)] = dsk.fit(x_train, y_train, 
                            batch_size=batch_size,
                            epochs=epoch, 
                            callbacks=callbacks,
                            shuffle=False)
            dsk_dict[str(n)] = dsk    
            stop = time.time() # Record time when training finished
            duration = round(stop - start, 2)
            
            
            ############ Get evaluation scores ############
            y_pred_train = np.exp(-dsk.predict(x_train))
            c_index_train = deepsurvk.concordance_index(y_train, y_pred_train, e_train)
            
            Y_pred_test = np.exp(-dsk.predict(X_test))
            c_index_test = deepsurvk.concordance_index(Y_test, Y_pred_test, E_test)
            
            train_scores.append(c_index_train)
            test_scores.append(c_index_test)
            model_time.append(duration)
            
            if print_scores:
                print(f'{n} training samples:\tTrain C-Index: {c_index_train}  | Test C-Index: {c_index_test}')
        
        model_results = pd.DataFrame({
            'n train': self.subset, 
            'train time': model_time,
            'train score': train_scores, 
            'test score': test_scores}
            )   
        return history_dict, dsk_dict, model_results
    
