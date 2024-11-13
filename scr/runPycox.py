from scr.pycox.models.cox import *
from sklearn.base import BaseEstimator

import torchtuples as tt


    
class PyCoxDeepSurv(BaseEstimator):
    def __init__(
        self,
        learning_rate=1e-4,
        batch_norm=True,
        dropout=0.0,
        num_nodes=[32, 16],
        batch_size=128,
        epochs=50,
    ):
        self.learning_rate = learning_rate
        self.batch_norm = batch_norm
        self.dropout = dropout
        self.num_nodes = num_nodes
        self.batch_size = batch_size
        self.epochs = epochs
        
    def build_model(self):
        self.modelString = 'pycox'
        if self.hyperparams is None:
            self.hyperparams = {}
            
            
    def fit(self, X, y):
        self.net_ = tt.practical.MLPVanilla(
            X.shape[1],
            self.num_nodes,
            self.epochs,
            self.batch_norm,
            self.dropout,
            output_bias=True,
        )
        self.deepsurv_ = CoxPH(self.net_, tt.optim.Adam)
        self.deepsurv_.optimizer.set_lr(self.learning_rate)

        # Sklearn needs the y inputs to be arranged as a matrix with each row
        # corresponding to an example but CoxPH needs a tuple with two arrays?
        y_ = (y[:, 0], y[:, 1])

        callbacks = [tt.callbacks.EarlyStopping()]
        log = self.deepsurv_.fit(
            X,
            y_,
            self.batch_size,
            self.epochs,
            verbose=False,
        )

        return self

    def score(self, X, y):
        _ = self.deepsurv_.compute_baseline_hazards()
        surv = self.deepsurv_.predict_surv_df(X)

        ev = EvalSurv(
            surv,
            y[:, 0],  # time to event
            y[:, 1],  # event
            censor_surv="km",
        )

        return ev.concordance_td()