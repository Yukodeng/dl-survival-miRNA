import numpy as np
from sklearn.preprocessing import StandardScaler, Normalizer

# ----------------  Utility Functions --------------------

def dataframe_to_deepsurv_ds(df, 
                            is_deepsurvk=False,
                            transformer='Normalizer',
                            event_col='status', time_col='time'):
    """Process input data (pandas DataFrames) to formats that are 
    recognizable by DeepSurv (dictionary) or DeepSurvK (arrays) for model training.

    Args:
        df (pandas DataFrame): input data
        is_deepsurvk (bool, optional): whether to return format for DeepSurvK. Defaults to False.
        time_col (str, optional): name of the time column. Defaults to 'time'.
        status_col (str, optional): name of the event status column. Defaults to 'status'.
        
    Returns:
    """
    # Extract the event and time columns as numpy arrays
    e = df[event_col].values.astype(np.int32)
    t = df[time_col].values.astype(np.float32)

    # Extract the patient's predictors as a numpy array
    x_df = df.drop([event_col, time_col], axis = 1)
    x = x_df.values.astype(np.float32)
    
    # Return the deepsurv data format
    if is_deepsurvk:
        ####NOTE: standardization might lead to invalid loss calculations
        ##because of very small time after transformation.
        #####Q: justification for also sorting the Y time? 
        ##(P.S. this results in all 0's due to very rightly-skewed distribution of survival times)
        if transformer == 'Normalizer': 
            ##### Normalization (0-1)
            x = Normalizer().fit_transform(x)
            t = Normalizer().fit_transform(t.reshape(1,-1)).reshape(-1)
        elif transformer == 'StandardScaler': # Standardization
            x = StandardScaler().fit_transform(x)
            t = StandardScaler().fit_transform(t.reshape(1,-1))
        # Sorting
        sort_idx = np.argsort(t)[::-1]
        x = x[sort_idx]
        t = t[sort_idx]
        e = e[sort_idx]
        
        return x, t, e 
    else:
        return {
            'x' : x,
            'e' : e,
            't' : t
        }
        
    
def dataframe_to_scikitsurv_ds(df, time_col='time', status_col='status'):
    """Convert input data (pandas DataFrame) to scikit-survival 
    compatible format for model training.

    Args:
        df (pandas DataFrame): input data
        time_col (str, optional): name of the column. Note that column values are 0 (censor) and 1 (event). Defaults to 'time'.
        status_col (str, optional): _description_. Defaults to 'status'.

    Returns:
        X: predictors (numpy array)
        y: survival outcomes (numpy structured array)
    """
    df.loc[:,'status_bool'] = [True if e==1 else False for e in df[status_col]]
    X_cols = [col for col in df.columns if col not in [time_col, status_col, 'status_bool']]

    X = np.array(df.loc[:, X_cols])
    y = np.array(df[["status_bool", time_col]].apply(tuple, 1), 
                dtype=[(status_col, '?'), (time_col, '<f8')])

    return X, y

