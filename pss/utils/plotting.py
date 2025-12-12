import os
import re
import glob
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns


# ----------------- Plotting Functions ---------------------------

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
        today = datetime.today().strftime("%m%d%y")
        plot_path = os.path.join("results", "plots", today)
        os.makedirs(plot_path, exist_ok=True)
        plt.savefig(os.path.join(plot_path, f"norm-{normType}.png"), bbox_inches="tight", dpi=300)
    plt.show()
    
    return(results)


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