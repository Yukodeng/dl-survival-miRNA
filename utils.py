import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns
import glob
import logging
from sklearn.preprocessing import StandardScaler, Normalizer


# --------------------- Load Data ---------------------------

class StreamToLogger:
    """Redirects sys.stdout and sys.stderr to a logger."""
    def __init__(self, logger, log_level=logging.INFO):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ''
        
    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.log_level, line.rstrip())
            
    def flush(self):
        pass  # Needed for compatibility


def create_logger(log_path, verbose_modules=['sklearn_pandas']):
    """
    Set up logger and redirect stdout/stderr to it.
    
    Args:
        config (dict): Config dict containing 'batchNormType' and 'dataName'.
        log_dir (str): Directory where logs are saved.
        verbose_modules (list): Optional list of module names to quiet.
    
    Returns:
        logging.Logger: Configured logger object
    """ 
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_path),
            logging.StreamHandler()
        ]
    )
    
    logger = logging.getLogger(__name__)
    
    # Quiet specified submodules
    if verbose_modules is None:
        verbose_modules = []
    for mod in verbose_modules:
        logging.getLogger(mod).setLevel(logging.WARNING)
        
    # Redirect print() output
    sys.stdout = StreamToLogger(logger, logging.INFO)
    sys.stderr = StreamToLogger(logger, logging.ERROR)
    
    return logger



def load_simulate_survival_data(batchNormType,
                                dataName,
                                keywords='',
                                keep_batch=False,
                                batch_col='batch.id'):    
    # prepare simulated survival data
    keywords = [keywords] if not isinstance(keywords, list) else keywords
    
    surv_folder = os.path.join('data', batchNormType, dataName)
    mirna_folder = os.path.join('data', batchNormType)
    
    found_surv = False
    for file in os.listdir(surv_folder):
        if all([kw in file for kw in keywords+['simSurvival', 'train']]):
            surv_train = pd.read_csv(os.path.join(surv_folder, file))
            found_surv = True      
        if all([kw in file for kw in keywords+['simSurvival', 'test']]):
            surv_test = pd.read_csv(os.path.join(surv_folder, file))
            found_surv = True    
    if not found_surv:
        raise FileExistsError("No survival data file found with the given query keywords")
    
    found_mirna = False
    for file in os.listdir(mirna_folder):
        if all([kw in file for kw in keywords+['simGeneExp', 'train']]):
            x_train = pd.read_csv(os.path.join(mirna_folder, file))
            found_mirna = True      
        if all([kw in file for kw in keywords+['simGeneExp', 'test']]):
            x_test = pd.read_csv(os.path.join(mirna_folder, file))
            found_mirna = True    
    if not found_mirna:
        raise FileExistsError("No gene expression data file found with the given query keywords")
    
    train_df = pd.concat([surv_train, x_train], axis=1)
    test_df = pd.concat([surv_test, x_test], axis=1)
    
    if not keep_batch:
        train_df = train_df.drop(columns=[batch_col])
        test_df = test_df.drop(columns=[batch_col])
    
    return train_df, test_df

# def load_simulate_survival_data(surv_file_path=None,
#                                 miRNA_file_path=None,
#                                 N=20000,
#                                 folder='', 
#                                 keywords='',
#                                 # initial_split=True,
#                                 test_size=0.1, 
#                                 val_split=False, val_size=0.25, 
#                                 random_state=42, 
#                                 time_col='time', status_col='status',
#                                 save_data=False):
#     if surv_file_path is None:
#         data_folder = os.path.join('data', folder)
#         # prepare simulated survival data
#         keywords = [keywords] if not isinstance(keywords, list) else keywords
#         keywords = [str(kw) for kw in keywords]
    
#         found = False
#         for file in os.listdir(data_folder):
#             if all([kw in file for kw in keywords]):
#                 surv_df = pd.read_csv(
#                     os.path.join(data_folder, file)).reset_index(drop=True)[[time_col,status_col]]
#                 found = True       
#         if not found:
#             raise FileExistsError("No data file found with the given query keywords")
#     else:
#         try:
#             surv_df = pd.read_csv(surv_file_path, index_col=0).reset_index(drop=True)
#         except:
#             raise FileNotFoundError("Data path not found!")

#     # load gene expression data and concat into 1 data frame
#     file_path = f"simulate_survival_{N}" if miRNA_file_path is None else miRNA_file_path
#     x_df = pd.read_csv(os.path.join("data", file_path)).reset_index(drop=True)
#     data_df = pd.concat([x_df, surv_df], axis=1)
#     # if initial_split:
#     # FIRST TIME: split data into train and test set
#     data = np.asarray(data_df)
#     train_df, test_df = train_test_split(data_df, 
#                                         test_size=test_size,
#                                         shuffle=True, random_state=random_state,
#                                         stratify=data_df[status_col])
#     train_df = train_df.reset_index(drop=True)
#     test_df  = test_df.reset_index(drop=True)

#     if val_split:
#         train_df, val_df = train_test_split(train_df, 
#                                         test_size=val_size,
#                                         shuffle=True, random_state=random_state,
#                                         stratify=train_df[status_col])
#         train_df = train_df.reset_index(drop=True)
#         val_df = val_df.reset_index(drop=True)
    
#     if save_data:
#         # Directly output train and test data
#         train_df.to_csv(
#             os.path.join("data", folder, "simulate_survival_train.csv"), index=False)
#         test_df.to_csv(
#             os.path.join("data", folder, "simulate_survival_test.csv"), index=False)
#         if val_split:
#             val_df.to_csv(
#                 os.path.join("data", folder, "simulate_survival_val.csv"), index=False)
#     if val_split:
#         return train_df, val_df, test_df
#     else:
#         return train_df, test_df


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


# ----------------- Plotting Functions ---------------------------

def plot_precisionSeq_results_onepanel(normType = "None", ylim = [0.45,0.95], save_plot=False):
    '''Function to load simulation results by model for a specified normalization method (e.g., TC).  
    '''

    # --- 1) Initialize results data frame -------
    results = pd.DataFrame(columns=['n train', 
                                    'train time',
                                    'train C',
                                    'test C',
                                    'train brier',
                                    'test brier',
                                    'model type',
                                    'data type'])
    
    batchNormTypes = [f"{s}_norm{normType}" for s in [
        "BE00Asso00",
        "BE10Asso00",
        "BE11Asso00",
        "BE10Asso10",
        "BE11Asso10",
        "BE11Asso11"
    ]]
    modelTypes = ['oracle-nl', 'oracle-linear', 'lasso', 'deepsurv-torch','rsf','svm','gb']
    dataNames = ['linear-moderate', 'linear-weak','nl-quadratic', 'nl-shiftquad', 'nl-sine', 'nl-interaction']
    
    # Read in performance scores
    for batchNormType in batchNormTypes:
        for dataName in dataNames:
            for model in modelTypes:
                
                try:
                    file_dir = glob.glob(
                        os.path.join('models', batchNormType, dataName, model, f"model_results_*25.csv"))
                    filename = [f for f in file_dir if not os.path.basename(f).startswith("model_results_2025")][-1]
                    result_df = pd.read_csv(filename)
                except: 
                    try:
                        filename = glob.glob(
                            os.path.join('models', batchNormType, dataName, model, f"model_results_2025*.csv")
                        )[-1]
                        result_df = pd.read_csv(filename)
                    except: 
                        result_df = pd.DataFrame()
                        print(f"no results found for {batchNormType}-{dataName} {model} model.")
                
                if model == 'rsf':
                    file_dir_5000 = os.path.join('models', batchNormType, dataName, model, f"model_results_5000.csv")
                    file_dir_8000 = os.path.join('models', batchNormType, dataName, model, f"model_results_8000.csv")
                    # file_dir_10000 = os.path.join('models', batchNormType, dataName, model, f"model_results_10000.csv")
                    
                    if os.path.exists(file_dir_5000):
                        result_df = pd.concat([result_df, pd.read_csv(file_dir_5000)])
                    if os.path.exists(file_dir_8000):
                        result_df = pd.concat([result_df, pd.read_csv(file_dir_8000)])
                    # if os.path.exists(file_dir_10000):
                    #     result_df = pd.concat([result_df, pd.read_csv(file_dir_10000)])
                    
                if model in ['lasso', 'oracle-linear', 'oracle-nl']:
                    result_df.columns = ["n train","train C","test C"]
                    
                result_df['model type'] = model
                result_df['data type'] = dataName
                result_df['batchnorm type'] = batchNormType
                results = pd.concat([results, result_df], axis=0)
        
    # Return None if no results found for the condition
    if (results.shape[0] == 0):
        print(f"No results found for {batchNormType}.")
    
    
    # --- 2) Set up facet grid --------------------
    sns.set_theme(style="whitegrid", font_scale=1)
    colors = {
        'oracle-nl':  "#3E3E3E", #"#664848" # Brown
        'oracle-linear': "#5B5B5B",     # Medium gray
        'lasso': "#186fcc",             # Muted blue
        'deepsurv-torch': "#d01c43",    # Soft red
        'rsf': '#733381',
        'svm': "#589C48",
        'gb': "#F58024"
    }
    
    g = sns.FacetGrid(
        results, 
        col="data type", col_order=dataNames,
        row='batchnorm type', row_order=batchNormTypes,
        height=3.25, aspect=1.5, 
        sharey=True, sharex=True
    )
    
    # --- 3) Plot individual lines and corresponding error bars ---
    for i, batchNormType in enumerate(batchNormTypes):
        for j, dataName in enumerate(dataNames):
            
            ax = g.axes[i, j]
            subset = (results[
                    (results["batchnorm type"] == batchNormType) &
                    (results["data type"] == dataName)]
                .groupby(['data type', 'n train', 'model type'])
                .agg(mean=('test C', 'mean'),
                    sd=('test C', 'std'))
                .reset_index()
            ) 
            
            for model, sub in subset.groupby('model type'):
                
                if sub['mean'].isna().all(): # skip oracle-nl for linear simulations
                    continue 
                
                model = 'oracle-nl' if (model=='oracle-linear') & ('linear' in dataName) else model
                ax.errorbar(
                    sub['n train'], sub['mean'], yerr=sub['sd'], label=model, color=colors[model], 
                    linestyle=":" if model=="oracle-linear" else "-",
                    linewidth=2.25 if model=="oracle-linear" else 1.75,
                    marker='o', markersize=5, alpha=0.8
                )
                
                subtitleString = re.sub("_norm.*","",batchNormType)
                ax.set_title(f"{subtitleString} {dataName}")
                ax.set_ylim(ylim)
                ax.set_xticks([100,1000,2000,5000,8000,10000])
                ax.set_xticklabels(['100','1k','2k','5k','8k','10k'])
                ax.tick_params(axis='x', rotation=0)
                # if dataName == 'linear-moderate':
                #     handles, labels = ax.get_legend_handles_labels()
                
    color_handles = [
        mlines.Line2D([], [], color=c, marker="o", linewidth=2,
                    linestyle=":" if name=="oracle-linear" else "-", 
                    label=name.replace("deepsurv-torch","DeepSurv")
                            .replace("lasso","LASSO")
                            .replace("oracle-linear","Oracle (linear)")
                            .replace("oracle-nl","Oracle")
                            .replace("rsf","Random Survival Forest")
                            .replace("svm","Support Vector Machine")
                            .replace("gb","Gradient Boosting"))
        for name, c in colors.items()
    ]
    g.figure.suptitle(f"Normalization: {normType}")
    g.set_axis_labels("Training size", "Test C-index")
    g.figure.legend(
        color_handles, [h.get_label() for h in color_handles],
        loc="lower center", bbox_to_anchor=(0.5, -0.05), frameon=False,
        ncol=4, title_fontsize=15, fontsize=15, title="Model Type"
    )
    
    plt.tight_layout()
    if save_plot:
        today = date.today().strftime("%m%d%y")
        plot_path = os.path.join("results", "plots", today)
        os.makedirs(plot_path, exist_ok=True)
        plt.savefig(os.path.join(plot_path, f"norm-{normType}.png"), bbox_inches="tight", dpi=300)
    plt.show()
    
    return(results)


# ---------------- DeepSurv Utility Functions --------------------

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
                dtype=[(status_col, '?'), (time_col, '<f8')])

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