import os
import warnings
import numpy as np
import pandas as pd
import torch
import torchtuples as tt

# Import the entire base module
from . import base
from . import loss as Loss
from . import utils

import torchtuples.callbacks as cb
# from torchtuples.tupletree import tuplefy


def search_sorted_idx(array, values):
    '''For sorted array, get index of values.
    If value not in array, give left index of value.
    '''
    n = len(array)
    idx = np.searchsorted(array, values)
    idx[idx == n] = n-1 # We can't have indexes higher than the length-1
    not_exact = values != array[idx]
    idx -= not_exact
    if any(idx < 0):
        warnings.warn('Given value smaller than first value')
        idx[idx < 0] = 0
    return idx


class _CoxBase(base.SurvBase):
    duration_col = 'duration'
    event_col = 'event'

    def fit(self, input, target, batch_size=256, epochs=1, callbacks=None, verbose=True,
            num_workers=0, shuffle=True, metrics=None, val_data=None, val_batch_size=8224,
            **kwargs):
        """Fit  model with inputs and targets. Where 'input' is the covariates, and
        'target' is a tuple with (durations, events).
        
        Arguments:
            input {np.array, tensor or tuple} -- Input x passed to net.
            target {np.array, tensor or tuple} -- Target [durations, events]. 
        
        Keyword Arguments:
            batch_size {int} -- Elements in each batch (default: {256})
            epochs {int} -- Number of epochs (default: {1})
            callbacks {list} -- list of callbacks (default: {None})
            verbose {bool} -- Print progress (default: {True})
            num_workers {int} -- Number of workers used in the dataloader (default: {0})
            shuffle {bool} -- If we should shuffle the order of the dataset (default: {True})
            **kwargs are passed to 'make_dataloader' method.
    
        Returns:
            TrainingLogger -- Training log
        """
        self.training_data = tt.tuplefy(input, target)
        return super().fit(input, target, batch_size, epochs, callbacks, verbose,
                           num_workers, shuffle, metrics, val_data, val_batch_size,
                           **kwargs)
    
    # [UPDATE 10/27] Add argument batch_id to allow different baseline hazard calculation for stratified NN 
    def _compute_baseline_hazards(self, input, df, batch_ids, max_duration, batch_size, eval_=True, num_workers=0):
        raise NotImplementedError

    def target_to_df(self, target):
        durations, events = tt.tuplefy(target).to_numpy()
        df = pd.DataFrame({self.duration_col: durations, self.event_col: events}) 
        return df

    def compute_baseline_hazards(self, input=None, target=None, batch_ids=None, max_duration=None, sample=None, batch_size=8224,
                                set_hazards=True, eval_=True, num_workers=0):
        """Computes the Breslow estimates form the data defined by `input` and `target`
        (if `None` use training data).

        Typically call
        model.compute_baseline_hazards() after fitting.
        
        Keyword Arguments:
            input  -- Input data (train input) (default: {None})
            target  -- Target data (train target) (default: {None})
            max_duration {float} -- Don't compute estimates for duration higher (default: {None})
            sample {float or int} -- Compute estimates of subsample of data (default: {None})
            batch_size {int} -- Batch size (default: {8224})
            set_hazards {bool} -- Set hazards in model object, or just return hazards. (default: {True})
        
        Returns:
            pd.Series -- Pandas series with baseline hazards. Index is duration_col.
        """
        if (input is None) and (target is None):
            if not hasattr(self, 'training_data'):
                raise ValueError("Need to give a 'input' and 'target' to this function.")
            input, target = self.training_data
            
        df = self.target_to_df(target)#.sort_values(self.duration_col)
        if sample is not None:
            if sample >= 1:
                df = df.sample(n=sample)
            else:
                df = df.sample(frac=sample)
                
        input = tt.tuplefy(input).to_numpy().iloc[df.index.values]
        # update 10/27: add batch ids
        batch_ids = batch_ids[df.index.values] if batch_ids is not None else np.repeat(1, len(df))
        base_haz = self._compute_baseline_hazards(input, df, batch_ids, max_duration, batch_size,
                                                  eval_=eval_, num_workers=num_workers)
        if set_hazards:
            self.compute_baseline_cumulative_hazards(set_hazards=True, baseline_hazards_=base_haz)
        return base_haz

    def compute_baseline_cumulative_hazards(self, input=None, target=None, batch_ids=None, max_duration=None, sample=None,
                                            batch_size=8224, set_hazards=True, baseline_hazards_=None,
                                            eval_=True, num_workers=0):
        """See `compute_baseline_hazards. This is the cumulative version."""
        
        if ((input is not None) or (target is not None)) and (baseline_hazards_ is not None):
            raise ValueError("'input', 'target' and 'baseline_hazards_' can not both be different from 'None'.")
        if baseline_hazards_ is None:
            baseline_hazards_ = self.compute_baseline_hazards(input, target, batch_ids, max_duration, sample, batch_size,
                                                             set_hazards=False, eval_=eval_, num_workers=num_workers)
        
        assert all([bh.index.is_monotonic_increasing for bh in baseline_hazards_.values()]),\
            'Need index of baseline_hazards_ to be monotonic increasing, as it represents time.'
            
        bch = {k: bh.cumsum().rename('baseline_cumulative_hazards') for k, bh in baseline_hazards_.items()}
        # bch = (baseline_hazards_
        #         .cumsum()
        #         .rename('baseline_cumulative_hazards'))
        
        if set_hazards:
            self.baseline_hazards_ = baseline_hazards_
            self.baseline_cumulative_hazards_ = bch
            
        return bch


    def predict_cumulative_hazards(self, input, batch_id=None, max_duration=None, batch_size=8224, verbose=False,
                                   baseline_hazards_=None, eval_=True, num_workers=0):
        """See `predict_survival_function`."""
        
        if type(input) is pd.DataFrame:
            input = self.df_to_input(input)
        if baseline_hazards_ is None:
            if not hasattr(self, 'baseline_hazards_'):
                raise ValueError('Need to compute baseline_hazards_. E.g run `model.compute_baseline_hazards()`')
            baseline_hazards_ = self.baseline_hazards_
            
        if batch_id is None:
            raise ValueError("Need to specify `batch_id` to select the correct baseline hazard")
        
        assert baseline_hazards_[batch_id].index.is_monotonic_increasing,\
            'Need index of baseline_hazards_ to be monotonic increasing, as it represents time.'
            
        return self._predict_cumulative_hazards(input, max_duration, batch_size, verbose, baseline_hazards_[batch_id],
                                                eval_, num_workers=num_workers)

    def _predict_cumulative_hazards(self, input, max_duration, batch_size, verbose, baseline_hazards_,
                                    eval_=True, num_workers=0):
        raise NotImplementedError


    def predict_surv_df(self, input, batch_id=None, max_duration=None, batch_size=8224, verbose=False, baseline_hazards_=None,
                        eval_=True, num_workers=0):
        """Predict survival function for `input`. S(x, t) = exp(-H(x, t))
        Require computed baseline hazards.

        Arguments:
            input {np.array, tensor or tuple} -- Input x passed to net.

        Keyword Arguments:
            max_duration {float} -- Don't compute estimates for duration higher (default: {None})
            batch_size {int} -- Batch size (default: {8224})
            baseline_hazards_ {pd.Series} -- Baseline hazards. If `None` used `model.baseline_hazards_` (default: {None})
            eval_ {bool} -- If 'True', use 'eval' mode on net. (default: {True})
            num_workers {int} -- Number of workers in created dataloader (default: {0})

        Returns:
            pd.DataFrame -- Survival estimates. One columns for each individual.
        """
        cum_haz = self.predict_cumulative_hazards(input, batch_id, max_duration, batch_size, verbose, baseline_hazards_,
                                                  eval_, num_workers)
        return np.exp(-cum_haz)

    def predict_surv(self, input, batch_id=None, max_duration=None, batch_size=8224, numpy=None, verbose=False,
                     baseline_hazards_=None, eval_=True, num_workers=0):
        """Predict survival function for `input`. S(x, t) = exp(-H(x, t))
        Require compueted baseline hazards.

        Arguments:
            input {np.array, tensor or tuple} -- Input x passed to net.

        Keyword Arguments:
            max_duration {float} -- Don't compute estimates for duration higher (default: {None})
            batch_size {int} -- Batch size (default: {8224})
            numpy {bool} -- 'False' gives tensor, 'True' gives numpy, and None give same as input
                (default: {None})
            baseline_hazards_ {pd.Series} -- Baseline hazards. If `None` used `model.baseline_hazards_` (default: {None})
            eval_ {bool} -- If 'True', use 'eval' mode on net. (default: {True})
            num_workers {int} -- Number of workers in created dataloader (default: {0})

        Returns:
            pd.DataFrame -- Survival estimates. One columns for each individual.
        """
        surv = self.predict_surv_df(input, batch_id, max_duration, batch_size, verbose, baseline_hazards_,
                                    eval_, num_workers)
        surv = torch.from_numpy(surv.values.transpose())
        return tt.utils.array_or_tensor(surv, numpy, input)

    def save_net(self, path, **kwargs):
        """Save self.net and baseline hazards to file.

        Arguments:
            path {str} -- Path to file.
            **kwargs are passed to torch.save

        Returns:
            None
        """
        path, extension = os.path.splitext(path)
        if extension == "":
            extension = '.pt'
        super().save_net(path+extension, **kwargs)
        if hasattr(self, 'baseline_hazards_'):
            self.baseline_hazards_.to_pickle(path+'_blh.pickle')

    def load_net(self, path, **kwargs):
        """Load net and hazards from file.

        Arguments:
            path {str} -- Path to file.
            **kwargs are passed to torch.load

        Returns:
            None
        """
        path, extension = os.path.splitext(path)
        if extension == "":
            extension = '.pt'
        super().load_net(path+extension, **kwargs)
        blh_path = path+'_blh.pickle'
        if os.path.isfile(blh_path):
            self.baseline_hazards_ = pd.read_pickle(blh_path)
            self.baseline_cumulative_hazards_ = self.baseline_hazards_.cumsum()

    def df_to_input(self, df):
        input = df[self.input_cols].values
        return input
    

class _CoxPHBase(_CoxBase):
    def _compute_baseline_hazards(self, input, df_target, batch_ids=None, max_duration=None,
                                  batch_size=8224, eval_=True, num_workers=0):
                                 # self, input, df_target, max_duration, batch_size, batch_ids=None, eval_=True, num_workers=0):
        if max_duration is None:
            max_duration = np.inf
            
        # Compute exp(log-risk)
        expg = np.exp(self.predict(input, batch_size, True, eval_, num_workers=num_workers)).flatten()
        df_target = df_target.assign(expg=expg)
        
        # --------------- Non-stratified version ----------------
        if batch_ids is None:
            # Here we are computing when expg when there are no events.
            #   Could be made faster, by only computing when there are events.
            return  (df_target
                        .groupby(self.duration_col)
                        .agg({'expg': 'sum', self.event_col: 'sum'})
                        .sort_index(ascending=False)
                        .assign(expg=lambda x: x['expg'].cumsum())
                        .pipe(lambda x: x[self.event_col] / x['expg'])
                        .fillna(0.)
                        .iloc[::-1]
                        .loc[lambda x: x.index <= max_duration]
                        .rename('baseline_hazards'))
        
        # ----------- Stratified version (per batch) ------------

        batch_ids = np.asarray(batch_ids)
        unique_batches = np.unique(batch_ids)
        baseline_hazards_dict = {}

        for b in unique_batches:
            mask = (batch_ids == b)
            df_b = df_target.loc[mask]
            if df_b[self.event_col].sum() == 0:
                continue  # skip strata with no events

            base_haz_b = (
                df_b.groupby(self.duration_col)
                    .agg({'expg': 'sum', self.event_col: 'sum'})
                    .sort_index(ascending=False)
                    .assign(expg=lambda x: x['expg'].cumsum())
                    .pipe(lambda x: x[self.event_col] / x['expg'])
                    .fillna(0.)
                    .iloc[::-1]
                    .loc[lambda x: x.index <= max_duration]
                    .rename(f'baseline_hazards_batch_{b}')
            )
            baseline_hazards_dict[b] = base_haz_b  

        return baseline_hazards_dict

    def _predict_cumulative_hazards(self, input, max_duration, batch_size, verbose, baseline_hazards_,
                                    eval_=True, num_workers=0):
        '''
        Predict cumulative hazard H(t|x) for either:
         - a single baseline hazard (non-stratified)
         - a batch-specific baseline_hazards_ (if baseline_hazards_ is a dict)
        '''
        max_duration = np.inf if max_duration is None else max_duration
        expg = np.exp(self.predict(input, batch_size, True, eval_, num_workers=num_workers)).reshape(1, -1)

        # --------------- Non-stratified version ---------------
        if isinstance(baseline_hazards_, pd.Series):
            if baseline_hazards_ is self.baseline_hazards_:
                bch = self.baseline_cumulative_hazards_
            else:
                bch = self.compute_baseline_cumulative_hazards(set_hazards=False, 
                                                            baseline_hazards_=baseline_hazards_)
            bch = bch.loc[lambda x: x.index <= max_duration]
            expg = np.exp(self.predict(input, batch_size, True, eval_, num_workers=num_workers)).reshape(1, -1)
            return pd.DataFrame(bch.values.reshape(-1, 1).dot(expg), 
                                index=bch.index)
        
        # ----------- Stratified version (per batch) ------------
        elif isinstance(baseline_hazards_, dict):
            cumulative_hazards_dict = {}
            for b, haz in baseline_hazards_.items():
                bch = haz.cumsum().loc[lambda x: x.index <= max_duration]
                cumulative_hazards_dict[b] = pd.DataFrame(
                    bch.values.reshape(-1, 1).dot(expg), index=bch.index
                )
            return cumulative_hazards_dict

        else:
            raise ValueError("Invalid baseline_hazards_ type: expected pd.Series or dict of Series.")
        

    def partial_log_likelihood(self, input, target, g_preds=None, batch_size=8224, eps=1e-7, eval_=True,
                               num_workers=0):
        '''Calculate the partial log-likelihood for the events in datafram df.
        This likelihood does not sample the controls.
        Note that censored data (non events) does not have a partial log-likelihood.

        Arguments:
            input {tuple, np.ndarray, or torch.tensor} -- Input to net.
            target {tuple, np.ndarray, or torch.tensor} -- Target labels.

        Keyword Arguments:
            g_preds {np.array} -- Predictions from `model.predict` (default: {None})
            batch_size {int} -- Batch size (default: {8224})
            eval_ {bool} -- If 'True', use 'eval' mode on net. (default: {True})
            num_workers {int} -- Number of workers in created dataloader (default: {0})

        Returns:
            Partial log-likelihood.
        '''
        df = self.target_to_df(target)
        if g_preds is None:
            g_preds = self.predict(input, batch_size, True, eval_, num_workers=num_workers)
        return (df
                .assign(_g_preds=g_preds)
                .sort_values(self.duration_col, ascending=False)
                .assign(_cum_exp_g=(lambda x: x['_g_preds']
                                    .pipe(np.exp)
                                    .cumsum()
                                    .groupby(x[self.duration_col])
                                    .transform('max')))
                .loc[lambda x: x[self.event_col] == 1]
                .assign(pll=lambda x: x['_g_preds'] - np.log(x['_cum_exp_g'] + eps))
                ['pll'])


class CoxPH(_CoxPHBase):
    """Cox proportional hazards model parameterized with a neural net.
    This is essentially the DeepSurv method [1].

    The loss function is not quite the partial log-likelihood, but close.    
    The difference is that for tied events, we use a random order instead of 
    including all individuals that had an event at that point in time.

    Arguments:
        net {torch.nn.Module} -- A pytorch net.
    
    Keyword Arguments:
        optimizer {torch or torchtuples optimizer} -- Optimizer (default: {None})
        device {str, int, torch.device} -- Device to compute on. (default: {None})
            Preferably pass a torch.device object.
            If 'None': use default gpu if available, else use cpu.
            If 'int': used that gpu: torch.device('cuda:<device>').
            If 'string': string is passed to torch.device('string').

    [1] Jared L. Katzman, Uri Shaham, Alexander Cloninger, Jonathan Bates, Tingting Jiang, and Yuval Kluger.
        Deepsurv: personalized treatment recommender system using a Cox proportional hazards deep neural network.
        BMC Medical Research Methodology, 18(1), 2018.
        https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-018-0482-1
    """
    def __init__(self, net, optimizer=None, device=None, loss=None):
        if loss is None:
            loss = Loss.CoxPHLoss()
        super().__init__(net, loss, optimizer, device)


### Update 7/1/2025
class StratifiedDataset(torch.utils.data.Dataset):
    def __init__(self, x, durations, events, batch_ids):
        self.x = x
        self.durations = durations
        self.events = events
        self.batch_ids = batch_ids

    def __len__(self):
        return len(self.x)

    def __getitem__(self, idx):
        return self.x[idx], self.durations[idx], self.events[idx], self.batch_ids[idx]

class CoxPHStratified(_CoxPHBase):
    """Cox proportional hazards model parameterized with a neural net.

    The loss function is not quite the partial log-likelihood, but close.    
    The difference is that for we stratify events by batch (strata) when
    calculating partial log likelihood.

    Arguments:
        net {torch.nn.Module} -- A pytorch net.
    
    Keyword Arguments:
        optimizer {torch or torchtuples optimizer} -- Optimizer (default: {None})
        device {str, int, torch.device} -- Device to compute on. (default: {None})
            Preferably pass a torch.device object.
            If 'None': use default gpu if available, else use cpu.
        loss {function} -- Loss function to use, default is stratified_cox_ph_loss.
    """
    def __init__(self, net, optimizer=None, device=None, loss=None):
        self.batch_ids = None
        if loss is None:
            loss = Loss.CoxPHLossStratified()
        super().__init__(net, loss, optimizer, device)
        
    def compute_metrics(self, data, metrics=None):
        if metrics is None:
            metrics = self.metrics
        if self.loss is None and self.loss in metrics.values():
            raise RuntimeError("Need to set `self.loss`.")
        x, durations, events, batch_ids = data
        log_h = self.net(x)
        
        if torch.isnan(log_h).any():
            print("NaNs detected in log_h during compute_metrics()")
            return {name: float('nan') for name in metrics}
        
        return {name: metric(log_h, durations, events, batch_ids) for name, metric in metrics.items()}
    
    def fit_dataloader(self, dataloader, epochs=1, callbacks=None, verbose=True, metrics=None, val_dataloader=None):
        """
        Custom training loop for CoxPHStratified that reads (x, duration, event, batch_id) from DataLoader.
        
        Args:
            dataloader (DataLoader): training data loader returning 4-tuples.
            epochs (int): number of training epochs.
            callbacks (list): optional callbacks.
            verbose (bool): print training progress.
            metrics (dict): optional metrics.
            val_dataloader (DataLoader): optional validation dataloader (can be normal 2-tuple).
        
        Returns:
            TrainingLogger
        """
        self._setup_train_info(dataloader)
        self.metrics = self._setup_metrics(metrics)
        self.log.verbose = verbose
        self.val_metrics.dataloader = val_dataloader
        
        if callbacks is None:
            callbacks = []
        self.callbacks = cb.TrainingCallbackHandler(
            self.optimizer, self.train_metrics, self.log, self.val_metrics, callbacks
        )
        self.callbacks.give_model(self)

        stop = self.callbacks.on_fit_start()
        for _ in range(epochs):
            if stop:
                break
            stop = self.callbacks.on_epoch_start()
            if stop:
                break
            
            for x, durations, events, batch_ids in dataloader:
                stop = self.callbacks.on_batch_start()
                if stop:
                    break
                
                self.optimizer.zero_grad()
                log_h = self.net(x)
                loss = self.loss(log_h, durations, events, batch_ids)
                self.batch_loss = loss
                self.batch_metrics = {"loss": loss}
                
                loss.backward()
                stop = self.callbacks.before_step()
                if stop:
                    break
                self.optimizer.step()
                stop = self.callbacks.on_batch_end()
                if stop:
                    break
            else:
                stop = self.callbacks.on_epoch_end()
        self.callbacks.on_fit_end()
        
        return self.log

    def predict_surv_df(self, input, batch_ids=None, max_duration=None,
                        batch_size=8224, verbose=False, baseline_hazards_=None,
                        eval_=True, num_workers=0):
        """
        Predict survival curves. If stratified baselines have been computed,
        pass `batch_ids` (1D array-like, length = n_samples) to use the
        correct baseline hazard per stratum. Returns a pd.DataFrame of shape
        (time_grid x n_samples) with a common time index.

        Backward compatible:
        - If `batch_ids` is None and a single baseline (pd.Series) is stored, we
          fall back to _CoxBase.predict_surv_df behavior.
        - If baseline hazards are a dict (stratified) but batch_ids is None,
          we raise a clear error.
        """
        
        if baseline_hazards_ is None:
            if not hasattr(self, 'baseline_hazards_'):
                raise ValueError("Need to compute baseline hazards. Run `model.compute_baseline_hazards(...)` first.")
            baseline_hazards_ = self.baseline_hazards_
            
        # -------- Non-stratified (use base implementation) ---------
        if isinstance(baseline_hazards_, pd.Series):
            
            return super().predict_surv_df(input, max_duration=max_duration,
                                           batch_size=batch_size, verbose=verbose,
                                           baseline_hazards_=baseline_hazards_,
                                           eval_=eval_, num_workers=num_workers)
        
        # Stratified: need batch_ids
        if not isinstance(baseline_hazards_, dict):
            raise ValueError("`baseline_hazards_` must be either a pd.Series or a dict of Series.")

        if batch_ids is None:
            raise ValueError("Stratified prediction requires `batch_ids` for each sample.")
        
        
        expg = np.exp(self.predict(input, batch_size, True, eval_, num_workers=num_workers)).reshape(-1)
        
        # Collect all survival times across batches
        all_times = None
        per_batch_surv = {}
        for b, haz in baseline_hazards_.items():
            bch = haz.cumsum()   # H_0,b(t)
            per_batch_surv[b] = bch   
            all_times = bch.index if all_times is None else all_times.union(bch.index)
        all_times = all_times.sort_values()

        # Fill in predicted survial along time grid
        out = np.empty((len(all_times), input.shape[0]), dtype=float)
        for b in np.unique(batch_ids):
            idx = np.where(batch_ids == b)[0]
            bch = per_batch_surv[b]
            # align H0b onto union grid with forward-fill
            H0b = bch.reindex(all_times, method='ffill').fillna(0.0).values
            # survival = exp(- H0_b(t) * exp(g(x)))
            out[:, idx] = np.exp(- np.outer(H0b, expg[idx]))
            
        return pd.DataFrame(out, index=all_times)
