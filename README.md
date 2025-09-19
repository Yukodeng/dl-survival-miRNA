# Simulation Study

### **Data Harmonization & Survival Modeling Comparison**

Axes of comparison (note: run **20 iterations** for each scenario):

-   Batch effect and batch-survival sorting scenarios (BE{he_train}{he_test}Asso{sort_train}{sort_test})

    -   BE00Asso00

    -   BE10Asso00

    -   BE10Asso10

    -   BE11Asso00

    -   BE11Asso10

    -   BE11Asso01

    -   BE11Asso11

    -   BE11Asso1−1

-   Normalization methods: Input should be the count data and output also count data

    -   None

    -   Total Count (TC)

    -   Upper Quartile (UQ)

    -   Median

    -   TMM

    -   Quantile Normalization (input is log-count)

    -   DESeq

-   Survival time-gene association (linear/nonlinear)

    -   linear-moderate (10 strong genes)

    -   linear-weak (30 weak genes)

    -   nonlinear-quadratic

    -   nonlinear-shiftquad

    -   nonlinear-sine

    -   nonlinear-interaction

-   Training sample sizes:

    -   N train = 100, 500, 1000, 2000, 5000, 10000

### **Survival Models for Comparison**

-   Baseline CoxPH models (R)

    -   **Oracle model** (i.e. Cox regression using the correct true predictors and the correct functional form)

    -   **Oracle-linear model** (Coxregression using the true predictors and the linear form, which is not necessarily the correct functional form)

    -   CoxPH **Lasso** regression (lambda is lambda.min selected by k-fold cross-validation with cv.glmnet)

-   Established ML methods (Python {[scikit-survival](https://scikit-survival.readthedocs.io/en/stable/index.html)} package)

    -   Random Survival Forest (RSF): sksurv.RandomSurvivalForest

    -   Gradient Boosting Survival Analysis (SGB): sksurv.GradientBoostingSurvivalAnalysis

    -   Survival Support Vector Machine (SSVM): sksurv.FastKernelSurvivalSVM

-   Deep learning method:

    -   **DeepSurv-torch**: PyTorch implementation of the CoxPH-based multi-layer feedforward neural network ([DeepSurv](https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-018-0482-1); <https://github.com/jaredleekatzman/DeepSurv>)

## Survival Models Summary Table (for internal use)

| Model | Software | Package | Pipeline Implementation | C-index | Brier Score | Stratification | Notes |
|---------|---------|---------|---------|---------|---------|---------|---------|
| **CoxPH (Oracle)** | R v4.4.1 | `coxph` {survival} | `aug-batch-survival-pipeline.R` | `summary` {base} | *(pending implementation)* | BatMan *(pending implementation)* |  |
| **Lasso-Penalized CoxPH (Lasso)** | R v4.4.1 | `cv.glmnet` {glmnet} | `aug-batch-survival-pipeline.R` | `concordance.index` {survcomp} | *(pending implementation)* | BatMan *(pending implementation)* |  |
| **Random Survival Forests (RSF)** | Python v3.10.0 | `RandomSurvivalForest` {scikit-survival} | `runSurvModels.py` | `score` {scikit-survival} | `integrated_brier_score` {scikit-survival} | *(pending implementation)* |  |
| **Gradient Boosting Survival Analysis (GBSA)** | Python v3.10.0 | `GradientBoostingSurvivalAnalysis` {scikit-survival} | `runSurvModels.py` | `score` {scikit-survival} | `integrated_brier_score` {scikit-survival} | *(pending implementation)* |  |
| **Survival Support Vector Machine (SSVM)** | Python v3.10.0 | `FastKernelSurvivalSVM` {scikit-survival} | `runSurvModels.py` | `score` {scikit-survival} | NA | *(pending implementation)* |  |
| **DeepSurv-PyTorch** | Python v3.10.0 | `pycox` {PyTorch}; {optuna} | `runDeepSurvModels.py` | `concordance_td` {pycox} | `integrated_brier_score` {pycox} | StratifiedCoxPH (manual stratified partial loglik loss) | performance is lower than gradient boosting in most scenarios |
