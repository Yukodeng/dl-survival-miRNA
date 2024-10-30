import os
import time
from .utils import *
import pickle

from abc import ABC, abstractmethod
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, Normalizer
from sklearn.model_selection import GridSearchCV, KFold

from sksurv.linear_model import CoxPHSurvivalAnalysis, CoxnetSurvivalAnalysis
from sksurv.ensemble import RandomSurvivalForest, GradientBoostingSurvivalAnalysis


class SurvivalModelPipeline(ABC):
    def __init__(self, X_train, y_train, X_test, y_test, dataName=None, hyperparameters=None, random_state=42):
        # self.train_data = train_data
        # self.test_data = test_data
        self.X_train=X_train
        self.y_train=y_train
        self.X_test=X_test  
        self.y_test=y_test
        self.dataName = dataName
        self.subset = [50,200,500,1000,2000,5000,8000]
        self.hyperparameters = hyperparameters
        self.random_state = random_state
        self.model = None
        self.modelString = None
        
    
    @abstractmethod
    def build_model(self):
        pass
    
    # def preprocess_data(self):
    #     self.X_train = self.scaler.fit_transform(self.X_train)
    #     self.X_test = self.scaler.transform(self.X_test)
    
    def tune_hyperparameters(self, n_split, plotting_param=None):
        cv = KFold(n_splits=n_split, shuffle=True, random_state=self.random_state)
        model_cv = GridSearchCV(self.model,
            param_grid=self.hyperparameters,
            cv=cv,
            error_score=0.5,
            # n_jobs=-1,
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
    
    
    def train(self, subset=None, save_results=True, output_dir=None):
        """_summary_

        Args:
            subset (list, optional): _description_. Defaults to None.
            save_model (bool, optional): _description_. Defaults to True.

        Returns:
            _type_: _description_
        """
        subset = self.subset if subset is None else subset        
        model_dict = {}
        model_time, model_train_scores, model_test_scores = [],[],[]
        for n in subset:
            #=============== prepare data =================
            if n < self.X_train.shape[0]:
                x_train,_, y_train,_ = train_test_split(
                    self.X_train, self.y_train,
                    train_size=n/self.X_train.shape[0], 
                    shuffle=True, random_state=42, stratify=[val[0] for val in self.y_train]
                )
            else:
                x_train = self.X_train
                y_train = self.y_train
                
            #================ Fit model ===================
            start = time.time() # Record iteration start time
            # Train the model with subset of training set
            model = self.model.fit(x_train, y_train)
            stop = time.time() # Record time when training finished
            duration = round(stop - start, 2)
            
            #=========== Get evaluation metrics ===========
            train_sc = model.score(x_train, y_train)
            test_sc = model.score(self.X_test, self.y_test)
            
            print(f"N={n} Training time ({duration}s): Train C-Index: {round(train_sc,3)} | Test C-index: {round(test_sc,3)}")
            model_train_scores.append(train_sc)
            model_test_scores.append(test_sc)
            model_time.append(duration)
            
            model_dict[str(n)] = model
            
        # ============== Save results ===============
        model_results = pd.DataFrame({
            'n train': subset, 
            'train time': model_time,
            'train score': model_train_scores, 
            'test score': model_test_scores}
        )
        if save_results:
            output_dir=os.path.join('models', self.dataName, self.modelString) if output_dir is None else output_dir
            os.makedirs(output_dir, exist_ok=True)
            write_file_path = os.path.join(output_dir, 'model.results.txt')
            # save results as text file
            if 'txt' in write_file_path:
                model_results.to_csv(write_file_path, sep='\t')
            elif 'csv' in output_dir:
                model_results.to_csv(write_file_path, index=False)
            else:
                print('Please specify a file name with either a txt or csv extension.')  
            
        return model_results, model_dict
    
    
    def write(self, model_dict, output_dir=None):
        
        for n, model in model_dict.items():
            output_dir = os.path.join('models', self.dataName, self.modelString) if output_dir is None else output_dir
            os.makedirs(output_dir, exist_ok=True)
            
            write_file_path = os.path.join(output_dir, f'{n}.train.models.pkl')
            
            with open(write_file_path,'wb') as f:
                pickle.dump(model, f) 



class CoxPHModel(SurvivalModelPipeline):
    def build_model(self, **kwargs):
        self.model = CoxPHSurvivalAnalysis(**kwargs)
        self.modelString = 'coxph'

        if self.hyperparameters is None:
            self.hyperparameters = {}
            
            
class CoxPHElasticNetModel(SurvivalModelPipeline):
    def build_model(self, **kwargs):
        self.model = CoxnetSurvivalAnalysis(**kwargs)
        self.modelString = 'coxnet'
        
        if self.hyperparameters is None:
            self.hyperparameters = {
                'alphas': [[0.1, 0.5, 1.0]],
                'l1_ratio': [0.1, 0.5, 0.9],
            }


class RandomSurvivalForestModel(SurvivalModelPipeline):
    def build_model(self, **kwargs):
        self.model = RandomSurvivalForest(**kwargs)
        self.modelString = 'rsf'
        
        if self.hyperparameters is None:
            self.hyperparameters = {
                'n_estimators': [100,500,1000],
                'max_depth': [3, 5, 10],
                'min_samples_split': [5, 10, 15],
            }


class GradientBoostingSurvivalModel(SurvivalModelPipeline):
    def build_model(self, **kwargs):
        self.model = GradientBoostingSurvivalAnalysis(**kwargs)
        self.modelString = 'gb'

        if self.hyperparameters is None:
            self.hyperparameters = {
                'learning_rate': [0.01, 0.1],
                'n_estimators': [100, 200],
                'max_depth': [3, 5, 10],
            }