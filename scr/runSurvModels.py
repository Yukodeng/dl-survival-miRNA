################ Baseline and Machine learning survival models #################
# Date         Notes
# 30Dec2024    Add function to allow multiple runs for each training sample size
# 03Feb2025    Add functions to save/load sksurv models and training time;
#              Add integrated Brier score calculation

import os
import time
import pickle
from tqdm import tqdm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import KFold, GridSearchCV, train_test_split
from sksurv.metrics import integrated_brier_score
from abc import ABC, abstractmethod
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.svm import FastKernelSurvivalSVM
from sksurv.ensemble import RandomSurvivalForest

class SurvivalModelPipeline(ABC):
    def __init__(self, X_train, y_train, X_test, y_test, dataName=None, hyperparameters=None):
        self.X_train = X_train
        self.y_train = y_train
        self.X_test = X_test  
        self.y_test = y_test
        self.dataName = dataName
        self.subset = [50, 500, 1000, 2000, 5000]
        self.runs = [20, 20, 20, 10, 5]
        self.hyperparameters = hyperparameters
        self.model = None
        self.modelString = None
        
    @abstractmethod
    def build_model(self):
        pass

    def tune_hyperparameters(self, n_split, plotting_param=None):
        cv = KFold(n_splits=n_split, shuffle=True, random_state=42)
        model_cv = GridSearchCV(self.model,
            param_grid=self.hyperparameters,
            cv=cv,
            error_score=0.5,
            n_jobs=-1,
        ).fit(self.X_train, self.y_train)
        
        self.model = model_cv.best_estimator_
        
        if plotting_param is not None:
            cv_results = pd.DataFrame(model_cv.cv_results_)
            try:
                params = cv_results[f'param_{plotting_param}'].map(lambda x: x[0])
            except TypeError:
                params = cv_results[f'param_{plotting_param}']
            mean = cv_results.mean_test_score
            std = cv_results.std_test_score
            
            _, ax = plt.subplots(figsize=(6,4))
            ax.plot(params, mean)
            ax.fill_between(params, mean - std, mean + std, alpha=0.15)
            ax.set_xscale("log")
            ax.set_xlabel(plotting_param)
            ax.set_ylabel("concordance index")
            ax.axvline(model_cv.best_params_[plotting_param][0] 
                            if isinstance(model_cv.best_params_[plotting_param], list)
                            else model_cv.best_params_[plotting_param], 
                        c="C1")
            ax.grid(True)
        return model_cv
    
    
    def train(self, subset=None, runs=None, save_model=True):
        """
        Args:
            subset (list, optional): _description_. Defaults to None.
            save_model (bool, optional): _description_. Defaults to True.
        Returns:
            data frame of evaluation scores of different training sizes
        """
        subset = self.subset if subset is None else subset
        runs   = self.runs if runs is None else runs        
        n_train, model_time, model_train_cind, model_test_cind, model_train_brier, model_test_brier = [],[],[],[],[],[]
        
        for n, n_run in zip(subset, runs):
            print(f"{n} training samples...")
            
            run_time, run_train_cind, run_test_cind, run_train_brier, run_test_brier = [],[],[],[],[]
            for run in tqdm(range(n_run)):
                #=============== Prepare data =================
                if n < self.X_train.shape[0]:
                    x_train, _, y_train, _ = train_test_split(
                        self.X_train, self.y_train,
                        train_size=n / self.X_train.shape[0], 
                        shuffle=True, random_state=run, 
                        stratify=[val[0] for val in self.y_train]
                    )
                else:
                    x_train = self.X_train
                    y_train = self.y_train
                    
                #================ Load model ==================
                model_dir = os.path.join('models', self.dataName, self.modelString, f'n-{n}')
                model_path = os.path.join(model_dir, f'run-{run}-model.pkl')
                # print(model_path) ### for sanity check

                if os.path.exists(model_path):
                    with open(model_path, 'rb') as f:
                        model_data = pickle.load(f)
                        model = model_data['model']
                        duration = model_data['training_time']
                else:
                    #============== Fit model ================
                    start = time.time() # Record iteration start time
                    model = self.model.fit(x_train, y_train)  # Train model with subset of training set
                    stop = time.time() # Record time when training finished
                    duration = round(stop - start, 2)
                    #============= Save model ================
                    if save_model:
                        os.makedirs(model_dir, exist_ok=True)
                        model_data = {
                            'model' : model,
                            'training_time': duration
                        }
                        with open(model_path, 'wb') as f:
                            pickle.dump(model_data, f)
                            
                #============= Get evaluation metrics =============
                ## C index scores
                train_c = model.score(x_train, y_train)
                test_c = model.score(self.X_test, self.y_test)
                
                ## Integrated Brier scores
                min_surv = np.ceil(max(np.min([val[1] for val in y_train]), np.min([val[1] for val in self.y_test])))
                max_surv = np.floor(min(np.max([val[1] for val in y_train]), np.max([val[1] for val in self.y_test])))
                times = np.linspace(min_surv, max_surv, 50)

                survs_tr = model.predict_survival_function(x_train)
                preds_tr = np.asarray([[fn(t) for t in times] for fn in survs_tr])
                train_brier = integrated_brier_score(y_train, y_train, preds_tr, times)

                survs_te = model.predict_survival_function(self.X_test)
                preds_te = np.asarray([[fn(t) for t in times] for fn in survs_te])
                try:
                    test_brier = integrated_brier_score(y_train, self.y_test, preds_te, times)
                except ValueError:
                    test_brier = np.nan

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

            print(f"N={n} Training time ({np.mean(run_time)}s):",
                f"\tC index:  Train ({round(np.mean(run_train_cind),3)})  Test ({round(np.mean(run_test_cind),3)})  |",
                f"  Brier score:  Train ({round(np.nanmean(run_train_brier),3)})  Test ({round(np.nanmean(run_test_brier),3)})")                
                
        # ============== Save results ===============
        model_results = pd.DataFrame({
            'n train': n_train, 
            'train time': model_time,
            'train C': model_train_cind, 
            'test C': model_test_cind,
            'brier score': model_train_brier,
            'test brier': model_test_brier
        })     
        return model_results
    
    
    def write(self, model_results, output_dir=None, fileName='model.results.txt'):
        output_dir = os.path.join('models', self.dataName, self.modelString) if output_dir is None else output_dir
        os.makedirs(output_dir, exist_ok=True)
        write_file_path = os.path.join(output_dir, fileName)
        
        # save results as text file
        if 'txt' in write_file_path:
            model_results.to_csv(write_file_path, sep='\t')
        elif 'csv' in output_dir:
            model_results.to_csv(write_file_path, index=False)
        else:
            print('Please specify a file name with either a txt or csv extension.')  


class CoxPHElasticNetModel(SurvivalModelPipeline):
    def build_model(self, **kwargs):
        # [2/3/2025] (YD) Ensure fit_baseline_model is set to True
        self.model = CoxnetSurvivalAnalysis(**kwargs, fit_baseline_model=True)
        self.modelString = 'coxnet'
        
        if self.hyperparameters is None:
            self.hyperparameters = {
                'alphas': [ [2.0**i] for i in np.arange(-6, 3, 2)],
                'l1_ratio': np.linspace(0.2, 1, 5),
            }
class SVMSurvivalModel(SurvivalModelPipeline):
    def build_model(self, **kwargs):
        self.model = FastKernelSurvivalSVM(**kwargs)
        self.modelString = 'svm'

        if self.hyperparameters is None:
            self.hyperparameters = {
                'alpha': 2.0**np.arange(-10, 2, 2),
                'kernel': ['linear', 'polynomial', 'cosine','rbf'], #[‘additive_chi2’, ‘chi2’, ‘linear’, ‘poly’, ‘polynomial’, ‘rbf’, ‘laplacian’, ‘sigmoid’, ‘cosine’]
                'rank_ratio': np.linspace(0, 1, 6)
            }
class RandomSurvivalForestModel(SurvivalModelPipeline):
    def build_model(self, **kwargs):
        self.model = RandomSurvivalForest(**kwargs)
        self.modelString = 'rsf'
        
        if self.hyperparameters is None:
            self.hyperparameters = {
                'n_estimators': [100, 200],
                'max_depth': [3, 5, 10],
                'min_samples_split': [5, 10, 15],
            }



# class SurvivalModelPipeline(ABC):
#     def __init__(self, X_train, y_train, X_test, y_test, dataName=None, hyperparameters=None, random_state=42):
#         self.X_train=X_train
#         self.y_train=y_train
#         self.X_test=X_test  
#         self.y_test=y_test
#         self.dataName = dataName
#         self.subset = [50,200,500,1000,2000,5000,8000]
#         self.hyperparameters = hyperparameters
#         self.random_state = random_state
#         self.model = None
#         self.modelString = None
        
    
#     @abstractmethod
#     def build_model(self):
#         pass
    
#     # def preprocess_data(self):
#     #     self.X_train = self.scaler.fit_transform(self.X_train)
#     #     self.X_test = self.scaler.transform(self.X_test)
    
#     def tune_hyperparameters(self, n_split, plotting_param=None):
#         cv = KFold(n_splits=n_split, shuffle=True, random_state=self.random_state)
#         model_cv = GridSearchCV(self.model,
#             param_grid=self.hyperparameters,
#             cv=cv,
#             error_score=0.5,
#             n_jobs=-1,
#         ).fit(self.X_train, self.y_train)
        
#         self.model = model_cv.best_estimator_
        
#         if plotting_param is not None:
#             cv_results = pd.DataFrame(model_cv.cv_results_)
#             try:
#                 params = cv_results[f'param_{plotting_param}'].map(lambda x: x[0])
#             except TypeError:
#                 params = cv_results[f'param_{plotting_param}']
#             mean = cv_results.mean_test_score
#             std = cv_results.std_test_score
            
#             _, ax = plt.subplots(figsize=(6,4))
#             ax.plot(params, mean)
#             ax.fill_between(params, mean - std, mean + std, alpha=0.15)
#             ax.set_xscale("log")
#             ax.set_xlabel(plotting_param)
#             ax.set_ylabel("concordance index")
#             ax.axvline(model_cv.best_params_[plotting_param][0] 
#                             if isinstance(model_cv.best_params_[plotting_param], list)
#                             else model_cv.best_params_[plotting_param], 
#                         c="C1")
#             ax.grid(True)
#         return model_cv
    
    
#     def train(self, subset=None, save_results=True, output_dir=None):
#         """_summary_

#         Args:
#             subset (list, optional): _description_. Defaults to None.
#             save_model (bool, optional): _description_. Defaults to True.

#         Returns:
#             _type_: _description_
#         """
#         subset = self.subset if subset is None else subset        
#         model_dict = {}
#         model_time, model_train_scores, model_test_scores = [],[],[]
#         for n in subset:
#             #=============== prepare data =================
#             if n < self.X_train.shape[0]:
#                 x_train,_, y_train,_ = train_test_split(
#                     self.X_train, self.y_train,
#                     train_size=n/self.X_train.shape[0], 
#                     shuffle=True, random_state=42, stratify=[val[0] for val in self.y_train]
#                 )
#             else:
#                 x_train = self.X_train
#                 y_train = self.y_train
                
#             #================ Fit model ===================
#             start = time.time() # Record iteration start time
#             # Train the model with subset of training set
#             model = self.model.fit(x_train, y_train)
#             stop = time.time() # Record time when training finished
#             duration = round(stop - start, 2)
            
#             #=========== Get evaluation metrics ===========
#             train_sc = model.score(x_train, y_train)
#             test_sc = model.score(self.X_test, self.y_test)
            
#             print(f"N={n} Training time ({duration}s): Train C-Index: {round(train_sc,3)} | Test C-index: {round(test_sc,3)}")
#             model_train_scores.append(train_sc)
#             model_test_scores.append(test_sc)
#             model_time.append(duration)
            
#             model_dict[str(n)] = model
            
#         # ============== Save results ===============
#         model_results = pd.DataFrame({
#             'n train': subset, 
#             'train time': model_time,
#             'train score': model_train_scores, 
#             'test score': model_test_scores}
#         )
#         if save_results:
#             output_dir=os.path.join('models', self.dataName, self.modelString) if output_dir is None else output_dir
#             os.makedirs(output_dir, exist_ok=True)
#             write_file_path = os.path.join(output_dir, 'model.results.txt')
#             # save results as text file
#             if 'txt' in write_file_path:
#                 model_results.to_csv(write_file_path, sep='\t')
#             elif 'csv' in output_dir:
#                 model_results.to_csv(write_file_path, index=False)
#             else:
#                 print('Please specify a file name with either a txt or csv extension.')  
            
#         return model_results, model_dict
    
    
#     def write(self, model_dict, output_dir=None):
        
#         for n, model in model_dict.items():
#             output_dir = os.path.join('models', self.dataName, self.modelString) if output_dir is None else output_dir
#             os.makedirs(output_dir, exist_ok=True)
            
#             write_file_path = os.path.join(output_dir, f'{n}.train.models.pkl')
            
#             with open(write_file_path,'wb') as f:
#                 pickle.dump(model, f) 



# class CoxPHModel(SurvivalModelPipeline):
#     def build_model(self, **kwargs):
#         self.model = CoxPHSurvivalAnalysis(**kwargs)
#         self.modelString = 'coxph'

#         if self.hyperparameters is None:
#             self.hyperparameters = {}
            
            
# class CoxPHElasticNetModel(SurvivalModelPipeline):
#     def build_model(self, **kwargs):
#         self.model = CoxnetSurvivalAnalysis(**kwargs)
#         self.modelString = 'coxnet'
        
#         if self.hyperparameters is None:
#             self.hyperparameters = {
#                 'alphas': [ [2.0**i] for i in np.arange(-6, 3, 2)],
#                 'l1_ratio': np.linspace(0.2, 1, 5),
#             }


# class SVMSurvivalModel(SurvivalModelPipeline):
#     def build_model(self, **kwargs):
#         self.model = FastKernelSurvivalSVM(**kwargs)
#         self.modelString = 'svm'

#         if self.hyperparameters is None:
#             self.hyperparameters = {
#                 'alpha': 2.0**np.arange(-10, 2, 2),
#                 'kernel': ['linear', 'polynomial', 'cosine'], #[‘additive_chi2’, ‘chi2’, ‘linear’, ‘poly’, ‘polynomial’, ‘rbf’, ‘laplacian’, ‘sigmoid’, ‘cosine’]
#                 'rank_ratio': np.linspace(0, 1, 6)
#             }


# class RandomSurvivalForestModel(SurvivalModelPipeline):
#     def build_model(self, **kwargs):
#         self.model = RandomSurvivalForest(**kwargs)
#         self.modelString = 'rsf'
        
#         if self.hyperparameters is None:
#             self.hyperparameters = {
#                 'n_estimators': [100,500,1000],
#                 'max_depth': [3, 5, 10],
#                 'min_samples_split': [5, 10, 15],
#             }


# class GradientBoostingSurvivalModel(SurvivalModelPipeline):
#     def build_model(self, **kwargs):
#         self.model = GradientBoostingSurvivalAnalysis(**kwargs)
#         self.modelString = 'gb'

#         if self.hyperparameters is None:
#             self.hyperparameters = {
#                 'learning_rate': [0.01, 0.1],
#                 'n_estimators': [100, 200],
#                 'max_depth': [3, 5, 10],
#             }