import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler, Normalizer
from sklearn.model_selection import train_test_split

def load_simulate_survival_data(file_path=None,
                                folder='', 
                                keywords=['10000'], 
                                initial_split=False, test_size=0.2, random_state=42, 
                                time_col='time',status_col='status',
                                save_data=False):
    # prepare simulated survival data
    keywords = [keywords] if not isinstance(keywords, list) else keywords
    keywords = [str(kw) for kw in keywords]
    
    if file_path is None:
        data_folder = os.path.join('data', folder)
        found = False
        for file in os.listdir(data_folder):
            if all([kw in file for kw in keywords]):
                surv_df = pd.read_csv(
                    os.path.join(data_folder, file)).reset_index(drop=True)[[time_col,status_col]]
                found = True       
        if not found:
            raise FileExistsError("No data file found with given query keywords!")
    else:
        try:
            surv_df = pd.read_csv(file_path, index_col=0).reset_index(drop=True)
        except:
            raise FileNotFoundError("Data path not found!")
    
    # load gene expression data and concat into 1 data frame
    x_df = pd.read_csv(os.path.join("data", "simulate_survival_10000.csv")).reset_index(drop=True)
    data_df = pd.concat([x_df, surv_df], axis=1)
    
    
    if initial_split:
        # FIRST TIME: split data into train (0.8) and test (0.2) set
        data = np.asarray(data_df)
        train_df, test_df = train_test_split(data, 
                                            test_size=test_size, 
                                            shuffle=True, random_state=random_state,
                                            stratify=data_df[status_col])
        pd.Series(train_df.index).to_csv(
            os.path.join("data", "train_index.csv"), index=False)
        pd.Series(test_df.index).to_csv(
            os.path.join("data", "test_index.csv"), index=False)
        
    train_ind = pd.read_csv(os.path.join('data', "train_index.csv")).iloc[:,0]
    test_ind = pd.read_csv(os.path.join('data', "test_index.csv")).iloc[:,0]
    
    train_df = data_df.iloc[train_ind,:]
    test_df = data_df.iloc[test_ind,:]
    
    if save_data:
        # Directly output train and test data
        train_df.to_csv(
            os.path.join("data", folder, "simulate_survival_train.csv"), index=False)
        test_df.to_csv(
            os.path.join("data", folder, "simulate_survival_test.csv"), index=False)
        
    return train_df, test_df


def load_simulate_results(dataFolderName, 
                        subset=[50, 200, 500, 1000, 2000, 5000, 8000],
                        modelnames=['coxnet','rsf','gb','deepsurvk']):
    
    results = pd.DataFrame({'n train': subset})
    for mdl in modelnames:
        file_dir = os.path.join('models', dataFolderName, mdl, 'model.results.txt')
        try:
            result_df = pd.read_table(file_dir, index_col=0)
        except:
            continue
        result_df.columns = ["_".join((col, mdl)) if col!='n train' else col for col in result_df.columns.tolist() ]
        results = results.merge(result_df, on='n train', how='outer')# right_on=result_df.columns[0])

    return results


############# DeepSurv Utility Functions ############

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
                    dtype=[('Status', '?'), ('Survival_in_months', '<f8')])

    return X, y



########## Scikit-surv Util Functions #########

def plot_coefficients(coefs, n_highlight=10): 
    """plot coefficent shrinkage at different alpha (regularization parameter) levels.

    Args:
        coefs (pandas dataframe): _description_
        n_highlight (int, optional): _description_. Defaults to 10.
    """
    _, ax = plt.subplots(figsize=(9, 6))
    n_features = coefs.shape[0]
    alphas = coefs.columns
    for row in coefs.itertuples():
        ax.semilogx(alphas, row[1:], ".-", label=row.Index)

    alpha_min = alphas.min()
    top_coefs = coefs.loc[:, alpha_min].map(abs).sort_values().tail(n_highlight)
    for name in top_coefs.index:
        coef = coefs.loc[name, alpha_min]
        plt.text(alpha_min, coef, name + "   ", horizontalalignment="right", verticalalignment="center")

    ax.yaxis.set_label_position("right")
    ax.yaxis.tick_right()
    ax.grid(True)
    ax.set_xlabel("alpha")
    ax.set_ylabel("coefficient")