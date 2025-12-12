# PRECISION.seq.survival: A simulation study

This repository provides documentation and source code for PRECISION.seq.survival, an extension of the PRECISION.seq (**P**ai**RE**d mi**C**rorna analys**I**s of differential expres**SION** for sequencing) analytic tool to assessment and comparison of different depth normalization methods and modeling techniques in the context of **survival prediction**–a problem crucial to clinical oncology research.

## Overview

------------------------------------------------------------------------

The axes of experimental comparisons of data generation and survival modeling scenarios are summarized below (Note: for each scenario, we ran modeling for **20 iterations**).

### **Simulation Scenarios: Batch Type & Harmonization**

-   **Batch effect** and **batch-survival association** (`BE{he_train}{he_test}Asso{sort_train}{sort_test}`)

    -   BE00Asso00

    -   BE10Asso00

    -   BE10Asso10

    -   BE11Asso00

    -   BE11Asso10

    -   BE11Asso11

-   **Marker-survival association** (linear/nonlinear): 10 true miRNA markers were pre-selected based on expression abundance observed in the real-world sarcoma tumor samples and used for survival time simulation.

    -   Moderate linear association (**linear-moderate**): 10 markers with strong signals

    -   Weak linear association (**linear-weak**): 30 markers with weak, spread-out signals (20 additional markers randomly drawn from the remaining markers)

    -   Nonlinear quadratic association (**nl-quadratic**): 10 markers standardized

    -   *Qudratic with shifted center by median (**nl-shiftquad**) (**optional**)*

    -   Sine transformation (**nl-sine**): 10 true markers with sine transformation

    -   Marker-marker **interaction** (**nl-interaction**): 10 markers

### **Depth Normalization & Survival Modeling**

-   **Normalization methods**: Input should be the count data and output also count data

    -   None (Raw)

    -   Total Count (TC)

    -   Upper Quartile (UQ)

    -   Median (median)

    -   Trimmed Mean of M-values (TMM)

    -   DESeq

    -   Quantile Normalization (QN)\*

-   **Survival risk prediction models**:

    -   Baseline CoxPH models (R)

        -   **Oracle model** (i.e. Cox regression using the correct true predictors and the correct functional form)

        -   **Oracle-linear model** (Coxregression using the true predictors and the linear form, which is not necessarily the correct functional form)

        -   CoxPH **Lasso** regression (lambda is lambda.min selected by k-fold cross-validation with cv.glmnet)

    -   Established ML methods (Python {[scikit-survival](https://scikit-survival.readthedocs.io/en/stable/index.html)} package)

        -   Random Survival Forest (**RSF**): sksurv.RandomSurvivalForest()

        -   Gradient Boosting Survival Analysis (**SGB**): sksurv.GradientBoostingSurvivalAnalysis()

        -   Survival Support Vector Machine (**SSVM**): sksurv.FastKernelSurvivalSVM()

    -   Deep learning method:

        -   **DeepSurv**: PyTorch implementation of the CoxPH-based multi-layer feed-forward neural network ([DeepSurv](https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-018-0482-1); <https://github.com/jaredleekatzman/DeepSurv>)

-   We additionally test out the effect of training sample size on prediction accuracy. For each simulation condition and survival modeling algorithm, a training sample size of 100, 500, 1000, 2000, 5000, and 10000† was utilized.

*\*Input for QN is log2-transformed count values*

*†For random survival forest (RSF) models, due to long runtime for large samples, sample size of 8000 was used in place of 10000.*

## Requirements

------------------------------------------------------------------------

The project partly deploy To reproduce the analyses in this project, first install Python virtual environment with the `requirements.txt` file using the following:

``` bash
# Step 1: Create virtual environment
python3.10 -m venv <env-name>  

# Step 2: Activate virtual environment

    ## OS/Linux:
    source <env-name>/bin/activate

    ## Windows:
    <env-name>\Scripts\activate.bat

# Step 3: Install precision.seq.survival (pss) package
pip install -e .

# (Optional) Any extra Python dependencies
pip install -r requirements.txt   
```

## Examples

------------------------------------------------------------------------

### Machine Learning / Deep Learning Pipelines

The Python scripts support survival risk prediction modeling using the three machine learning models–**RSF, SGB, and SSVM**--and the deep learning **DeepSurv** model. Users need to input the training (`train_df`) and testing (`test_df`) datasets.

|  **time** | **status** | **batch.id** | **hsa.let.7a.3..1** | **hsa.let.7a..2..1** | **hsa.let.7b.1** | **hsa.let.7b..1** | **hsa.let.7c.1** |
|--------:|--------:|--------:|--------:|--------:|--------:|--------:|--------:|
| 18.433152 |          1 |            1 |           13.661359 |             7.973968 |         9.670266 |          5.803146 |        10.571855 |
|  4.052384 |          0 |            1 |            7.936586 |             4.287322 |         7.573842 |          1.784844 |         8.501820 |
|  3.377271 |          1 |            1 |           13.458518 |             9.188944 |        13.210251 |          8.940502 |        10.639829 |
|  7.692236 |          1 |            1 |            9.236919 |             3.635254 |         7.946339 |          2.516852 |         6.670609 |
|  0.333537 |          1 |            1 |           11.573337 |             6.200822 |         8.600134 |          5.083641 |         8.708288 |

The datasets should be formatted as pandas data frames (`pandas.DataFrame`), with rows representing tumor samples. Columns should include the survival time (continuous), survival status (0=censor, 1=event), and batch ID (categorical), with the rest storing the log-transformed, normalized miRNA marker expression values.

Once the dataset is loaded, we can initialize the modeling class object to implement any of the survival risk model of interest. For example, the SSVM model:

``` python
ml = MLSurvivalPipeline(
    model_type="ssvm",  # ML model type (`rsf`, `ssvm`, `sgb`)
    train_df=train_df,
    test_df=test_df, 
    batchNormType="BE00Asso00_normTC", # batch type and normalization type
    dataName='linear-moderate',        # marker-survival association type
    is_stratified=True                 # whether to apply batch-stratification in modeling 
)

results = train_over_subsets( 
    ml,
    subset_sizes=[100, 500, 1000, 2000, 5000, 10000], # number of training samples
    runs_per_size=[20, 20, 20, 20, 20, 20], # number of iterations for each training size
    splits_per_size=[3, 5, 5, 10, 10, 10],  # number of folds in cross-validation for each training size
    is_tune=False,
    is_save=True,
    file_name=...
)
```

### Hyperparameter Tuning

The ML/DL survival modeling pipelines also support hyperparameter tuning through integration with the `optuna` Python library. To enable tuning, specify the hyperparameter search grid (`hyperparameters`) and set `is_tune=True`.

By default, each training size condition will be tuned separately via a 30-trial semi-random search, and all tuning trials will be automatically saved to a SQLite database. (`*.db`) You can customize the file name through the `storage_url` argument and the number of trials for each training size condition through the `trials_per_size` argument. You can also change the `trial_threshold` parameter which let the pipeline bypass the tuning procedure if there is enough historically completed trials exist in the SQLite database.

``` python
hyperparameters = {
    "learning_rate": {"type": "float", "low": 1e-4, "high": 1e-2, "log": True},
    "num_nodes": {"type": "categorical", "choices": [[64], [32,16], [64,32], [128,64], [64,64,32],[32,32,16]]},
    "dropout": {"type": "float", "low": 0.1, "high": 0.5}
}

dl = DeepSurvPipeline(
    train_df, test_df, 
    batchNormType=..., 
    dataName=...,
    is_stratified=True,                 
    hyperparameters=hyperparameters,     # dictionary of optuna-compatible hyperparameter search grid 
    storage_url="sqlite:///example.db",  # SQLite databse with tuning trial information for individual studies
)

results = train_over_subsets(
    dl,
    subset_sizes=...,    
    runs_per_size=...,   
    splits_per_size=..., 
    trials_per_size=...  # num of trials (length should match with training size list) 
    trial_threshold=..,  # minimum number of existing trials to allow skipping tuning
    is_tune=True         # remember to set to True to allow optuna hyperparameter tuning
    is_save=True,        # save resutls
    file_name=...        # default to "model_results_mmddyy.csv"
)
```

## Survival Models Summary Table *(for internal use)*

| Model                                         | Software       | Package                                              | Pipeline Implementation         | C-index                        | Brier Score                                | Stratification                                                            | Notes                                                         |
|---------|---------|---------|---------|---------|---------|---------|---------|
| **CoxPH (Oracle)**                            | R v4.4.1       | `coxph` {survival}                                   | `aug-batch-survival-pipeline.R` | `summary` {base}               | *(pending implementation)*                 | BatMan                                                                    |                                                               |
| **Lasso-Penalized CoxPH (Lasso)**             | R v4.4.1       | `cv.glmnet` {glmnet}                                 | `aug-batch-survival-pipeline.R` | `concordance.index` {survcomp} | *(pending implementation)*                 | BatMan                                                                    |                                                               |
| **Random Survival Forests (RSF)**             | Python v3.10.0 | `RandomSurvivalForest` {scikit-survival}             | `runSurvModels.py`              | `score` {scikit-survival}      | `integrated_brier_score` {scikit-survival} | batch-informed modeling                                                   |                                                               |
| **Gradient Boosting Survival Analysis (SGB)** | Python v3.10.0 | `GradientBoostingSurvivalAnalysis` {scikit-survival} | `runSurvModels.py`              | `score` {scikit-survival}      | `integrated_brier_score` {scikit-survival} | batch-informed modeling                                                   |                                                               |
| **Survival Support Vector Machine (SSVM)**    | Python v3.10.0 | `FastKernelSurvivalSVM` {scikit-survival}            | `runSurvModels.py`              | `score` {scikit-survival}      | NA                                         | batch-informed modeling                                                   |                                                               |
| **DeepSurv-PyTorch (DeepSurv)**               | Python v3.10.0 | `pycox` {PyTorch}; {optuna}                          | `runDeepSurvModels.py`          | `concordance_td` {pycox}       | `integrated_brier_score` {pycox}           | StratifiedCoxPH (manual implementation of stratified partial loglik loss) | performance is lower than gradient boosting in most scenarios |

## References

TBD