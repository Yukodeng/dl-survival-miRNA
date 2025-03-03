import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler, Normalizer
from sklearn.model_selection import train_test_split

def load_simulate_survival_data(file_path=None,
                                folder='', 
                                keywords=[''],
                                N = 20000, 
                                # initial_split=True,
                                test_size=0.1, 
                                val_split=False, val_size=0.25, 
                                random_state=42, 
                                time_col='time', status_col='status',
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
            raise FileExistsError("No data file found with the given query keywords")
    else:
        try:
            surv_df = pd.read_csv(file_path, index_col=0).reset_index(drop=True)
        except:
            raise FileNotFoundError("Data path not found!")
    
    # load gene expression data and concat into 1 data frame
    x_df = pd.read_csv(os.path.join("data", f"simulate_survival_{N}.csv")).reset_index(drop=True)
    data_df = pd.concat([x_df, surv_df], axis=1)
    
    # if initial_split:
    # FIRST TIME: split data into train and test set
    data = np.asarray(data_df)
    train_df, test_df = train_test_split(data_df, 
                                        test_size=test_size,
                                        shuffle=True, random_state=random_state,
                                        stratify=data_df[status_col])
    train_df = train_df.reset_index(drop=True)
    test_df  = test_df.reset_index(drop=True)

    if val_split:
        train_df, val_df = train_test_split(train_df, 
                                        test_size=val_size,
                                        shuffle=True, random_state=random_state,
                                        stratify=train_df[status_col])
        train_df = train_df.reset_index(drop=True)
        val_df = val_df.reset_index(drop=True)
    
    if save_data:
        # Directly output train and test data
        train_df.to_csv(
            os.path.join("data", folder, "simulate_survival_train.csv"), index=False)
        test_df.to_csv(
            os.path.join("data", folder, "simulate_survival_test.csv"), index=False)
        if val_split:
            val_df.to_csv(
                os.path.join("data", folder, "simulate_survival_val.csv"), index=False)
    if val_split:
        return train_df, val_df, test_df
    else:
        return train_df, test_df


############## summary functions ###############
modelname_dict = {'coxnet':'CoxPH',
                'svm': 'SSVM',
                'rsf': 'RSF',
                'deepsurv-torch': "DeepSurv"}

def load_simulate_results(dataFolderName, 
                        modelnames=['coxnet','svm','rsf','deepsurv-torch'],
                        fileName = 'model.results.10runs.txt'):
    
    results = pd.DataFrame({'n train': []})
    for mdl in modelnames:
        file_dir = os.path.join('models', dataFolderName, mdl, fileName)
        try:
            result_df = pd.read_table(file_dir, index_col=0)
        except:
            result_df = pd.read_table(os.path.join('models', dataFolderName, mdl, 'model.results.txt'), index_col=0)
        result_df['model'] = modelname_dict[mdl]
        results = pd.concat([results, result_df], axis=0)

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
        

def plot_simulation_data(train_df, test_df):
    '''Function to plot simulated data.'''
    # observe data
    print("Event rate in train set: %f" % (sum(train_df['status']==1) / train_df.shape[0]))
    print("Event rate in test set: %f" %  (sum(test_df['status']==1) / test_df.shape[0]))
    print('Survival time distribution:')
    _, ax = plt.subplots(figsize=(3,3))
    ax.hist(train_df['time'], label='train')
    # ax.hist(val_df['time'],   label='val', alpha=.8)
    ax.hist(test_df['time'], label='test', alpha=0.6)
    ax.legend(fontsize=12)
    plt.show()
    
    
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