# ------------------------------------------------------------------
# Attempt Log
# Date      Description
# 11APR2025 MAF 10 batches separately (10000 per batch)
# 12DEC2025 CVAE generated 10 batches (500 per batch)
# 23DEC2025 CVAE generated 10 batches (10000 / 20000 per batch)
# 09JAN2026 CVAE generated 15 batches (1000 per batch) & 
#           VAE re-generated batch 6 and 8 (10000 per batch)
# 26JAN2026 CVAE generated 15 batches (1000 per batch)
# 18FEB2026 CVAE generated 15 batches with partial universal constant BE (1000 per batch) 
# ------------------------------------------------------------------

import os
import random
import numpy as np
import pandas as pd
import umap
from sklearn.preprocessing import StandardScaler
from scipy.stats import iqr
import matplotlib.pyplot as plt
import seaborn as sns


# Load Data --------------------------------------------------------

attempt = "02182026"
output_dir = os.path.join("results", "plots", "evaluation", attempt)
os.makedirs(output_dir, exist_ok = True)

N = 20000 #500   # number of samples per batch
P = 538          # number of markers
n_batch = 15
rwd_max = 22.516867 # rwd max expression (log2)
rwd_sd = 4.039547   # rwd sd of all (log2) counts
gnames = pd.read_csv(os.path.join("data", "augmentation_miRNA_names_538.csv"))['gene'].tolist()      # gene names
gindex = pd.read_csv(os.path.join("data", "augmentation_miRNA_names_538.csv"))['index_py'].tolist()  # gene index
gnames_no1 = [g[:-2] for g in gnames] # cleaned marker names


# pd.DataFrame(gindex).to_csv("gindex.csv", row)


# * True markers ----------------------------------
nonzero_genes = ["hsa.miR.1277.3p.1", "hsa.miR.1277.5p.1", "hsa.miR.133a.2..1", "hsa.miR.144..1" , 
                 "hsa.miR.148b..1" ,  "hsa.miR.181c..1" ,  "hsa.miR.204.1", "hsa.miR.33a.1",
                 "hsa.miR.598.1", "hsa.miR.887.1" ]


# * Benchmark data ---------------------------------

rwd_clean = pd.read_csv("raw-data/MSKCC_sarcoma_benchmark_data.csv")
rwd_batch = pd.read_csv("raw-data/MSKCC_sarcoma_test_data.csv")
rwd_clean = np.log2(rwd_clean.loc[:, gnames_no1] + 1)
rwd_batch = np.log2(rwd_batch.loc[:, gnames_no1] + 1)
rwd_clean.columns = rwd_batch.columns = gnames


# * Parametric simulation --------------------------

# ## [added 12/8] extra steps for parametrically simulated data with *ALL* 1033 columns
# par_clean = pd.read_csv(os.path.join("raw-data", "with_without_batch_parametric_15_batch_02182026.csv")).iloc[:, :1033]
# par_batch = pd.read_csv(os.path.join("raw-data", "with_without_batch_parametric_15_batch_02182026.csv")).iloc[:, 1033:2066]
# batch_id  = pd.read_csv(os.path.join("raw-data", "with_without_batch_parametric_15_batch_02182026.csv")).iloc[:, -1]
# pd.concat([
#     par_clean.iloc[:, gindex],
#     par_batch.iloc[:, gindex],
#     batch_id
# ], axis=1).to_csv(os.path.join("raw-data", "with_without_batch_parametric_15_batch_538_markers_02182026.csv"),index=False)

## [NOTE] use below code if pre-filtered to have 538 markers
## [NOTE] !!THE LATEST ATTEMPT STORES CLEAN DATA FIRST, DIRTY DATA LATER!!
par = pd.read_csv(os.path.join("raw-data", 
                               "with_without_batch_parametric_15_batch_538_markers_02182026.csv"
                  ))
par_clean = np.log2(par.iloc[:54, :P] + 1)
par_batch = np.log2(par.iloc[:, P:(2*P)] + 1)

par_clean.columns = par_batch.columns = gnames


# * Augmented data ---------------------------------

# ## Attempt (11APR2025): MAF
# aug_data = pd.DataFrame()
# for i in range(n_batch):
#     # read in each batch
#     dat_sub = pd.read_csv(
#         os.path.join("raw-data", "with_without_batch_maf", f"batch-{i+1}_epochES_maf_generated.csv"),
#         header=None
#     ).iloc[:N, :(P*2)]
#     # merge batches
#     dat_sub.columns = gnames + gnames
#     aug_data = pd.concat([aug_data, dat_sub], axis=0)

aug_count = pd.read_csv(
        # os.path.join("raw-data", "with_without_batch_CVAE_500perbatch_20251212.csv")     # <- Attempt 12DEC2025: 500 per batch
        # os.path.join("raw-data", "with_without_batch_CVAE_10000perbatch_20251223.csv")   # <- Attempt 23DEC2025
        # os.path.join("raw-data", "with_without_15_batch_CVAE_1000perbatch_01092026.csv") # <- all 1033 markers; need subsetting
        # os.path.join("raw-data", "with_without_15_batch_CVAE_1000perbatch_01262026.csv") # <- all 1033 markers; need subsetting
        # os.path.join("raw-data", "with_without_15_batch_CVAE_1000perbatch_02182026.csv") # <- Attempt 18FEB2026; 538 markers
        os.path.join("raw-data", "with_without_15_batch_CVAE_20000perbatch_02182026.csv") # <- Attempt 18FEB2026; 538 markers
)

# # [NOTE] 09JAN2026 regenerated batch 6 and 8
# aug_batch6 = pd.read_csv( os.path.join("raw-data", "batch6_VAE_01092026.csv"))
# aug_batch8 = pd.read_csv( os.path.join("raw-data", "batch8_VAE_01092026.csv"))
# aug_count = (
#     pd.concat([
#         aug_count.loc[(aug_count["batch.id"] != 6) & (aug_count["batch.id"] != 8), :],
#         aug_batch6,
#         aug_batch8], axis=0)
#     .set_index('batch.id', drop=False)
#     .loc[sorted(aug_count['batch.id'].unique())]
#     .reset_index(drop=True)
#     .drop(columns = 'batch.id')
# )

# # [NOTE] If augmented data has all 1033 columns, run below to subset to selected markers
# aug_count = pd.concat([
#     aug_count.iloc[:, gindex],                  # clean
#     aug_count.iloc[:, [i+1033 for i in gindex]] # dirty
# ], axis=1)

# [NOTE] if generated data are based on per-filtered 538 markers
aug_count = aug_count.iloc[:, :-1]

## Cap extremely large values
def cap_and_add_noise(x, cap, sd, rng=None):
    x = np.asarray(x, dtype=float)
    new_cap = cap + 0.1 * sd
    rng = np.random.default_rng() if rng is None else rng

    mask = np.isinf(x) | (x > new_cap)
    xs = x[mask]
    if xs.size > 0:
        noise = rng.normal(loc=0.0, scale=0.01 * sd, size=xs.size)
        xs_hat = np.minimum(xs, new_cap + noise)
        x[mask] = xs_hat
    return x

aug_data_log = np.log2(aug_count + 1)
aug_data_cap = pd.DataFrame(np.stack(
    aug_data_log.apply(cap_and_add_noise, cap=rwd_max, sd=rwd_sd, axis=1)
))
aug_data_cap.columns = gnames + gnames



## Set color palette
palette = { # sns.color_palette("husl", 2)
    "RWD":       "#fcc200", # orange #"#33aaff", # Blue 
    "Augmented": "#1e90ff", # blue   
    "Parametric": "#9457eb" # purple #"#2a8000"  # green  #ff75dd"  # pink
}



# Per-Batch Assessment ----------------------------------------------------------------


# * Median vs. IQR plot (MSKCC vs. Augmented) ------------------

fig, _ = plt.subplots(n_batch, 3, figsize=(11, 3.6*n_batch), constrained_layout=True)
for i in range(n_batch):
    print(f"Loading batch {i+1}: row {i*N} to {(i+1)*N}..")

    # Read in each batch
    # [NOTE] LATEST ATTEMPT HAS CLEAN DATA FIRST, DIRTY DATA SECOND
    aug_clean = aug_data_cap.iloc[i*N:(i+1)*N, :P]
    aug_batch = aug_data_cap.iloc[i*N:(i+1)*N:, P:]

    par_batch_i = par_batch.iloc[i*54:(i+1)*54, ]

    per_marker_clean = pd.DataFrame({"Type": [], 
                                    "median": [],
                                    "IQR": []})
    per_marker_clean['Marker'] = list(gnames + gnames)
    per_marker_clean['Type']   = list(np.repeat('RWD', P)) + list(np.repeat('Augmented', P))
    per_marker_clean['median'] = list(rwd_clean.apply(np.median)) + list(aug_clean.apply(np.median))
    per_marker_clean['IQR']    = list(rwd_clean.apply(iqr)) + list(aug_clean.apply(iqr))
    per_marker_clean['true'] = ['True' if g in nonzero_genes else 'None' for g in per_marker_clean.Marker]

    per_marker_batch = pd.DataFrame({"Type": [], 
                                    "median": [],
                                    "IQR": []})
    per_marker_batch['Marker'] = list(gnames + gnames + gnames)
    per_marker_batch['Type']   = list(np.repeat('RWD', P)) + list(np.repeat('Parametric', P)) + list(np.repeat('Augmented', P))
    per_marker_batch['median'] = list(rwd_batch.apply(np.median)) + list(par_batch_i.apply(np.median)) + list(aug_batch.apply(np.median))
    per_marker_batch['IQR']    = list(rwd_batch.apply(iqr)) + list(par_batch_i.apply(iqr)) + list(aug_batch.apply(iqr))
    per_marker_batch['true'] = ['True' if g in nonzero_genes else 'None' for g in per_marker_batch.Marker]

    ## Scatter Plot
    ax1 = plt.subplot(n_batch, 3, i*3+1)
    ax2 = plt.subplot(n_batch, 3, i*3+2)
    ax3 = plt.subplot(n_batch, 3,(i+1)*3)

    ## Clean
    sns.scatterplot(data=per_marker_clean,
                    x="median", y="IQR", hue="Type",
                    palette=palette, alpha=0.12, marker='o', s=15, ax=ax1)
    sns.scatterplot(data=per_marker_clean.loc[per_marker_clean['true']=='True',:],
                    x="median", y="IQR", hue="Type", 
                    palette=palette, alpha=1, marker='X', s=28, ax=ax1)
    ax1.legend([], [], frameon=False)
    ax1.set_title(f"Batch {i+1} (Clean)")
    ax1.set_ylim(0, 7)

    ## Batch (MSK vs. Aug)
    sns.scatterplot(data=per_marker_batch.loc[(per_marker_batch['Type']=='RWD') | 
                                              (per_marker_batch['Type']=='Augmented'), :],
                    x="median", y="IQR", hue="Type",
                    palette=palette, alpha=0.12, marker='o', s=15, ax=ax2)
    sns.scatterplot(data=per_marker_batch.loc[((per_marker_batch['Type']=='RWD') | 
                                               (per_marker_batch['Type']=='Augmented')) & 
                                             (per_marker_batch['true']=='True'), :],
                    x="median", y="IQR", hue="Type", 
                    palette=palette, alpha=1, marker='X', s=28, ax=ax2)
    ax2.legend([], [], frameon=False)
    ax2.set_title(f"Batch {i+1} (Batch): MSKCC vs. Augmented")
    ax2.set_ylim(0, 7)


    ## Batch (Parametric vs. Aug)
    sns.scatterplot(data=per_marker_batch.loc[(per_marker_batch['Type']=='Parametric') | 
                                              (per_marker_batch['Type']=='Augmented'), :],
                    x="median", y="IQR", hue="Type",
                    palette=palette, alpha=0.12, marker='o', s=15, ax=ax3)
    sns.scatterplot(data=per_marker_batch.loc[((per_marker_batch['Type']=='Parametric') | 
                                               (per_marker_batch['Type']=='Augmented')) & 
                                             (per_marker_batch['true']=='True'), :],
                    x="median", y="IQR", hue="Type", 
                    palette=palette, alpha=1, marker='X', s=28, ax=ax3)
    ax3.legend([], [], frameon=False)
    ax3.set_title(f"Batch {i+1} (Batch): Input vs. Augmented")
    ax3.set_ylim(0, 7)

    handles = ax2.get_legend_handles_labels()[0] + [ax3.get_legend_handles_labels()[0][0]]
    labels  = ax2.get_legend_handles_labels()[1] + [ax3.get_legend_handles_labels()[1][0]]
        
fig.legend(handles, labels, loc='lower center', ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.01))
# fig.subplots_adjust(top=0.8)
fig.suptitle("Median vs. IQR", fontsize=16)
# fig.tight_layout()
plt.show()
# save plot 
plt.savefig(os.path.join(output_dir, f"scatter-median-IQR-per-batch-{attempt}.png"),
            dpi=300, bbox_inches='tight')



# * Mean vs. Stdev plot (MSKCC vs. Augmented) ----------------------

fig, _ = plt.subplots(n_batch, 3, figsize=(11, 3.6*n_batch), constrained_layout = True)

for i in range(n_batch):
    
    print(f"Loading batch {i+1}: row {i*N} to {(i+1)*N}..")

    # Read in each batch
    # [NOTE] LATEST ATTEMPT HAS CLEAN DATA FIRST, DIRTY DATA SECOND
    aug_clean = aug_data_cap.iloc[i*N:(i+1)*N, :P]
    aug_batch = aug_data_cap.iloc[i*N:(i+1)*N:, P:]
    par_batch_i = par_batch.iloc[i*54:(i+1)*54, :]

    per_marker_clean = pd.DataFrame({"Type": [], 
                                    "Mean": [],
                                    "Stdev": []})
    per_marker_clean['Marker'] = list(gnames + gnames)
    per_marker_clean['Type']   = list(np.repeat('RWD', P)) + list(np.repeat('Augmented', P))
    per_marker_clean['Mean']   = list(rwd_clean.apply(np.mean)) + list(aug_clean.apply(np.mean))
    per_marker_clean['Stdev']  = list(rwd_clean.apply(np.std)) + list(aug_clean.apply(np.std))
    per_marker_clean['true']   = ['True' if g in nonzero_genes else 'None' for g in per_marker_clean.Marker]

    per_marker_batch = pd.DataFrame({"Type": [], 
                                    "Mean": [],
                                    "Stdev": []})
    per_marker_batch['Marker'] = list(gnames + gnames + gnames)
    per_marker_batch['Type']   = list(np.repeat('RWD', P)) + list(np.repeat('Parametric', P)) + list(np.repeat('Augmented', P))
    per_marker_batch['Mean']   = list(rwd_batch.apply(np.mean)) + list(par_batch_i.apply(np.mean)) + list(aug_batch.apply(np.mean))
    per_marker_batch['Stdev']  = list(rwd_batch.apply(np.std)) + list(par_batch_i.apply(np.std)) + list(aug_batch.apply(np.std))
    per_marker_batch['true']   = ['True' if g in nonzero_genes else 'None' for g in per_marker_batch.Marker]

    ## Scatter Plot
    ax1 = plt.subplot(n_batch, 3, i*3+1)
    ax2 = plt.subplot(n_batch, 3, i*3+2)
    ax3 = plt.subplot(n_batch, 3, (i+1)*3)

    ## Clean
    sns.scatterplot(data=per_marker_clean,
                    x="Mean", y="Stdev", hue="Type",
                    palette=palette, alpha=0.15, marker='o', s=15, ax=ax1)
    sns.scatterplot(data=per_marker_clean.loc[per_marker_clean['true']=='True',:],
                    x="Mean", y="Stdev", hue="Type", 
                    palette=palette, alpha=1, marker='X', s=28, ax=ax1)
    ax1.legend([], [], frameon=False)
    ax1.set_title(f"Batch {i+1} (Clean)")
    ax1.set_ylim(0, 5)

    ## Batch (MSK vs. Aug)
    sns.scatterplot(data=per_marker_batch.loc[(per_marker_batch['Type']=='RWD') | 
                                            (per_marker_batch['Type']=='Augmented'), :],
                    x="Mean", y="Stdev", hue="Type",
                    palette=palette, alpha=0.15, marker='o', s=15, ax=ax2)
    sns.scatterplot(data=per_marker_batch.loc[((per_marker_batch['Type']=='RWD') | 
                                            (per_marker_batch['Type']=='Augmented')) &
                                            (per_marker_batch['true']=='True'), :],
                    x="Mean", y="Stdev", hue="Type", 
                    palette=palette, alpha=1, marker='X', s=28, ax=ax2)
    ax2.legend([], [], frameon=False)
    ax2.set_title(f"Batch {i+1} (Batch): MSKCC vs. Augmented")
    ax2.set_ylim(0, 5)

    ## Batch (Parametric vs. Aug)
    sns.scatterplot(data=per_marker_batch.loc[(per_marker_batch['Type']=='Parametric') | 
                                            (per_marker_batch['Type']=='Augmented'), :],
                    x="Mean", y="Stdev", hue="Type",
                    palette=palette, alpha=0.12, marker='o', s=15, ax=ax3)
    sns.scatterplot(data=per_marker_batch.loc[((per_marker_batch['Type']=='Parametric') | 
                                                (per_marker_batch['Type']=='Augmented')) & 
                                            (per_marker_batch['true']=='True'), :],
                    x="Mean", y="Stdev", hue="Type", 
                    palette=palette, alpha=1, marker='X', s=28, ax=ax3)
    ax3.legend([], [], frameon=False)
    ax3.set_title(f"Batch {i+1} (Batch): Input vs. Augmented")
    ax3.set_ylim(0, 5)

    handles = ax2.get_legend_handles_labels()[0] + [ax3.get_legend_handles_labels()[0][0]]
    labels  = ax2.get_legend_handles_labels()[1] + [ax3.get_legend_handles_labels()[1][0]]
      
fig.legend(handles, labels, loc='lower center', ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.01))
fig.suptitle("Mean vs. Stdev", fontsize=16)
plt.show()
# save plot 
plt.savefig(os.path.join(output_dir, f"scatter-mean-stdev-per-batch-{attempt}.png"),
            dpi=300, bbox_inches='tight')




# UMAP -------------------------------------------------------------------

random.seed(123)
n_downsample = 54

# * MSKCC vs. Parametric vs. Augmented --------------------------------

fig, _ = plt.subplots(n_batch, 3, figsize=(11, 3.6*n_batch), constrained_layout = True)
for i in range(n_batch):
  
    print(f"Loading batch {i+1}..downsampling...")
    
    ## Downsampling augmented data (N=540)
    idx = random.sample([i for i in range(N)], n_downsample)
    aug_sub = aug_data_cap.iloc[[i*N + id for id in idx], :]
    ## [NOTE] LATEST ATTEMPT HAS CLEAN DATA FIRST, DIRTY DATA SECOND
    aug_clean = aug_sub.iloc[:, :P]
    aug_batch = aug_sub.iloc[:, P:]

    par_batch_sub = par_batch.iloc[i*54:(i+1)*54, :]
    
    rwd_clean["Type"] = "RWD"
    aug_clean["Type"] = "Augmented"
    clean = pd.concat([rwd_clean, aug_clean]).reset_index(drop=True)

    rwd_batch["Type"] = "RWD"
    par_batch_sub["Type"] = "Parametric"
    aug_batch["Type"] = "Augmented"

    batch = pd.concat([rwd_batch, par_batch_sub, aug_batch]).reset_index(drop=True)
        
    ## Create umap projection
    umap_obj = umap.UMAP(n_components=2, random_state=42)

    clean_scaled = StandardScaler().fit_transform(clean.drop(columns='Type')) 
    embedding_clean = umap_obj.fit_transform(clean_scaled)
    embedding_clean = pd.concat([
        pd.DataFrame(embedding_clean, columns = ['UMAP1', "UMAP2"]),
        clean[['Type']] ], axis=1)

    batch_scaled = StandardScaler().fit_transform(batch.drop(columns='Type')) 
    embedding_batch = umap_obj.fit_transform(batch_scaled)
    embedding_batch = pd.concat([
        pd.DataFrame(embedding_batch, columns = ['UMAP1', "UMAP2"]),
        batch[['Type']]] , axis=1)

    ## UMAP plots
    ax1 = plt.subplot(n_batch, 3, i*3+1)
    ax2 = plt.subplot(n_batch, 3, i*3+2)
    ax3 = plt.subplot(n_batch, 3, (i+1)*3)
    
    # Clean
    sns.scatterplot(data=embedding_clean,#.loc[embedding_clean['Type']=='Augmented',:], 
                    x="UMAP1", y="UMAP2", hue="Type", 
                    palette=palette, alpha=0.75, marker='o', s=20, ax=ax1)
    ax1.legend([], [], frameon=False)
    ax1.set_title(f"Batch {i+1} (Clean)")

    # Batch (MSKCC vs. Aug)
    sns.scatterplot(data=embedding_batch.loc[(embedding_batch['Type']=='RWD')|
                                            (embedding_batch['Type']=='Augmented'),:], 
                    x="UMAP1", y="UMAP2", hue="Type", 
                    palette=palette, alpha=0.75, marker='o', s=20, ax=ax2)
    ax2.legend([], [], frameon=False)
    ax2.set_title(f"Batch {i+1} (Batch): MSKCC vs. Aug")

    # Batch (Par vs. Aug)
    sns.scatterplot(data=embedding_batch.loc[(embedding_batch['Type']=='Parametric') |
                                            (embedding_batch['Type']=='Augmented'), :], 
                    x="UMAP1", y="UMAP2", hue="Type", 
                    palette=palette, alpha=0.75, marker='o', s=20, ax=ax3)
    ax3.legend([], [], frameon=False)
    ax3.set_title(f"Batch {i+1} (Batch): Input vs. Aug")

    handles = ax2.get_legend_handles_labels()[0] + [ax3.get_legend_handles_labels()[0][0]]
    labels  = ax2.get_legend_handles_labels()[1] + [ax3.get_legend_handles_labels()[1][0]]

fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, -0.01), ncol=3, frameon=False)
plt.suptitle("UMAP Projection", fontsize=16)
plt.show()
# save plot 
plt.savefig(os.path.join(output_dir, f"umap-projection-per-batch-downsample-{n_downsample}-{attempt}.png"),
            dpi=300, bbox_inches='tight')




# Batch Effect Assessment -------------------------------------------------------

## extract pure batch effects
rwd_batch = rwd_batch.drop(columns='Type')
rwd_clean = rwd_clean.drop(columns='Type')
rwd_batch_obs = (rwd_batch - rwd_clean)

par_clean = np.log2(par.iloc[:, :P] + 1)
par_batch = np.log2(par.iloc[:, P:(2*P)] + 1)
par_clean.columns = par_batch.columns = gnames
par_batch_obs = (par_batch - par_clean)

aug_clean = aug_data_cap.iloc[:, :P]
aug_batch = aug_data_cap.iloc[:, P:]
aug_clean.columns = aug_batch.columns = gnames
aug_batch_obs = (aug_batch - aug_clean) 


# par_clean.to_csv("par_clean.csv", index=False)
# par_batch.to_csv("par_batch.csv", index=False)

# aug_clean.to_csv("aug_clean.csv", index=False)
# aug_batch.to_csv("aug_batch.csv", index=False)

# par_batch_obs.to_csv("par_batch_obs.csv", index=False)
# aug_batch_obs.to_csv("aug_batch_obs.csv", index=False)


## 0) UMAP of pure batch effects [added 12/17]

par_batch_obs["Type"] = "Parametric"
aug_batch_obs["Type"] = "Augmented"
batch = pd.concat([par_batch_obs, aug_batch_obs]).reset_index(drop=True)

batch_scaled = StandardScaler().fit_transform(batch.drop(columns='Type')) 
embedding_batch = pd.DataFrame(umap_obj.fit_transform(batch_scaled))
embedding_batch.columns = ['UMAP1',"UMAP2"]
embedding_batch = pd.concat([
    embedding_batch.loc[:, ['UMAP1', "UMAP2"]],
    batch[['Type']] ] , axis=1)

## plot
fig, _ = plt.subplots(1, figsize=(5, 4.5))
sns.scatterplot(data=embedding_batch.loc[embedding_batch['Type']=='Augmented',:], 
                x="UMAP1", y="UMAP2", hue="Type", 
                palette=palette, alpha=0.75, marker='o', s=16)
sns.scatterplot(data=embedding_batch.loc[embedding_batch['Type']=='Parametric',:], 
                x="UMAP1", y="UMAP2", hue="Type", 
                palette=palette, alpha=0.75, marker='o', s=16)
plt.legend(frameon=False)
plt.title("Batch Effects: Before vs. After Augmentation")
plt.savefig(os.path.join(output_dir,"umap-batch-effects-par-vs-aug.png"), dpi=300)



## 1) Histogram distribution

## MSKCC only

fig, ax = plt.subplots(1)
ax.scatter(x=rwd_clean.apply(np.std), y=rwd_batch_obs.apply(np.median),
           color='blue', edgecolors='white', alpha=0.8)
plt.xlabel('stdev')
plt.ylabel('median batch')
plt.title("MSKCC Sarcoma Samples")
plt.show()
plt.savefig(os.path.join("results", "plots", "mskcc-stdev-vs-batch-median.png"), 
            dpi=300, bbox_inches='tight')


## Parametric vs. Aug

nrow = int(np.ceil(n_batch / 4))

fig, axes = plt.subplots(nrow, 4, figsize=(15, 3*nrow))
c = 0
for i in range(nrow):
    for j in range(4):
        if c < n_batch:
            axes[i, j].hist((par_batch_obs
                            .iloc[c*54:(c*54+54),:-1]
                            .apply(np.median)), 
                    bins=30, alpha=0.7, color='#e28743', edgecolor='black', label='Parametric')
            axes[i,j].hist((aug_batch_obs
                            .iloc[c*N:(c+1)*N, :-1]
                            .apply(np.median)),
                    bins=30, alpha=0.7, color='#2596be', edgecolor='black', label='Augmented')
            axes[i,j].set_title(f"Batch {c+1}")
            axes[i,j].set_xlim([-5,5])
            axes[i,j].set_ylim([0, 80])
        else:
            axes[i,j].set_axis_off()
        c += 1
handles, lables = axes[0,0].get_legend_handles_labels()
fig.legend(handles, lables, loc="lower right", bbox_to_anchor=(0.5, -0.02), ncol=2, frameon=False)
plt.suptitle("Marker-specific batch effects (log2): Parametric vs. Augmented")
plt.tight_layout()
plt.show()
# save fig
plt.savefig(os.path.join(output_dir, f"histplot-marker-median-batch-per-batch-par-vs-aug-{N}.png"), 
            dpi=300, bbox_inches='tight')



# ## 2) Batch-Marker Stdev correlation
# par_batch_10_obs = par_batch_obs.iloc[ [c*54 for c in range(10)],: ]
# par_batch_std = par_batch_10_obs.apply(np.std).tolist()
# rwd_batch_std = rwd_batch_obs.apply(np.std).tolist()

# fig = plt.figure(figsize=(6,5))
# data = pd.DataFrame({
#     'marker std': rwd_clean.apply(np.std).tolist(),
#     'par batch std': par_batch_std,
#     'rwd batch std': rwd_batch_std
#     })
# sns.scatterplot(data=data, x='rwd batch std', y='par batch std')
# plt.title("Marker-specific batch stdev: RWD vs. Parametric")
# plt.show()



## 2) Inter-marker batch correlaiton

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15,4.5))

## RWD
rwd_clean_true = rwd_clean.loc[:, nonzero_genes]
rwd_batch_true = rwd_batch.loc[:, nonzero_genes]
rwd_clean_true.columns = rwd_batch_true.columns = nonzero_genes
### NOTE: extract batch on the natural log scale
rwd_batch_obs_sel = (rwd_batch_true - rwd_clean_true)
corr_batch_rwd = rwd_batch_obs_sel.corr(method='pearson')
sns.heatmap(corr_batch_rwd, cmap="coolwarm", center=0, vmin=-0.2, vmax=1, ax=ax1)
ax1.set_title("MSKCC inter-marker batch correlatioon")

## parametric
par_clean_true = par_clean.loc[:, nonzero_genes]
par_batch_true = par_batch.loc[:, nonzero_genes]
par_batch_obs_sel = (par_batch_true - par_clean_true)

corr_batch_par = par_batch_obs_sel.corr(method='pearson')
sns.heatmap(corr_batch_par, cmap="coolwarm", center=0, vmin=-0.2, vmax=1, ax=ax2)
ax2.set_title("Parametric inter-marker batch correlatioon")

## Augmented
aug_batch_true = aug_batch.loc[:, nonzero_genes]
aug_clean_true = aug_clean.loc[:, nonzero_genes]
aug_batch_obs_sel = (aug_batch_true - aug_clean_true)

corr_batch_aug = aug_batch_obs_sel.corr(method='pearson')
sns.heatmap(corr_batch_aug, cmap="coolwarm", center=0, vmin=-0.2, vmax=1, ax=ax3)
ax3.set_title("Augmented per-marker batch correlatioon")

plt.tight_layout()
plt.show()
plt.savefig(os.path.join(output_dir, "heatmap-marker-marker-batch-correlation.png"), dpi=300)


#  3) ICC ------------
# [NOTE] Dr. Robert M. Hamer of Virginia Commonwealth University has supplied the code for computing
# the Shrout and Fleiss intraclass coefficients:
# https://web.archive.org/web/20090303152247/http://www.nyu.edu/its/statistics/Docs/intracls.html
#           bms= between target mean square
#           wms= within target mean square
#           jms= mean square for judges(raters)
#           ems= error mean square
#           k  = number of judges

#           bms = ss/df for targets
#           msw = ((ems*edf)+(jms*jdf))/(edf+jdf)
#           wms = msw
#           jms = ss/df for judges (raters)                   
#   sfsingle=(bms-wms)/(bms+(k-1)*wms)            * ICC(1,1)
#   sfrandom=(bms-ems)/
#     ((bms)+((k-1)*ems)+((k*(jms-ems))/n))       * ICC(2,1)
#   sffixed=(bms-ems)/(bms+((k-1)*ems))           * ICC(3,1)
#   sfk=(bms-wms)/bms                             * ICC(1,k)
#   sfrandk=(bms-ems)/(bms+((jms-ems)/n))         * ICC(2,k)
#   sffixedk=(bms-ems)/bms                        * ICC(3,k) 

def calculate_icc(data, value_col='marker', group_col='batch'):
    batches = data[group_col].astype("category")
    y = aug_batch_obs.iloc[:, g]
    
    n = len(y)    # total number
    overall_mean = y.mean()
    
    # Batch-level counts, means, variances
    means_i = data.groupby(batches, observed=False)[value_col].mean().sort_index().values
    vars_i  = data.groupby(batches, observed=False)[value_col].var(ddof=1).sort_index().values
    n_i = batches.value_counts().sort_index().values  # number of samples (judges)
    k = len(n_i)  # number of batches (subjects)
    # sum of squares
    ss_between = np.sum(n_i * (means_i - overall_mean) ** 2) # Between-batch SS
    ss_within = np.sum((n_i - 1) * vars_i)                   # Within-batch SS
    
    # degrees of freedom
    df_between = k - 1
    df_within = n - k
    
    # mean squared of error
    ms_between = ss_between / df_between
    ms_within = ss_within / df_within
    
    n_bar = (n - np.sum(n_i ** 2) / n) / (k - 1)
    
    # ICC(1) - per-marker one-way random effect of batch 
    # [NOTE] checked by Andy on 12/28/2025
    icc = (ms_between - ms_within) / (ms_between + (n_bar - 1) * ms_within)
    
    return icc

# [NOTE] xxx_batch_obs are the pure log2-transformed batch effects (N samples x P markers)
aug_icc_list, par_icc_list = [],[]
for g in range(P):
    aug_icc = calculate_icc(data = pd.DataFrame({
        'marker': aug_batch_obs.iloc[:, g], 
        'batch': np.repeat([i+1 for i in range(n_batch)], N)
        })
    )
    par_icc = calculate_icc(data = pd.DataFrame({
        'marker': par_batch_obs.iloc[:, g],
        'batch': np.repeat([i+1 for i in range(n_batch)], 54)
        })
    )
    aug_icc_list.append(aug_icc)
    par_icc_list.append(par_icc)

aug_l, aug_med, aug_u = np.round(np.quantile(aug_icc_list, q = [0.25, 0.5, 0.75]), 2)
par_l, par_med, par_u = np.round(np.quantile(par_icc_list, q = [0.25, 0.5, 0.75]), 2)

icc_df = pd.DataFrame({
    'marker': gnames,
    'ICC_aug': aug_icc_list,
    'ICC_par': par_icc_list
})

## histogram
fig = plt.figure(figsize=(5,4))
plt.hist(icc_df['ICC_aug'], bins=30, alpha=0.7, color='#1e90ff', edgecolor='black')
plt.show()

## scatter plot
fig = plt.figure(figsize=(4.25,4))
sns.scatterplot(data=icc_df, x='ICC_aug', y='ICC_par', alpha=0.8, s=10)
plt.hlines(y = par_med, xmin=0, xmax=1, linestyles='--', color='#9457eb', label='Parametric (Median)')
plt.vlines(x = aug_med, ymin=0, ymax=1, linestyles='--', color='#1e90ff', label='Augmented (Median)')
plt.xlim(0,1)
plt.ylim(0,1)
plt.legend(fontsize=7)
plt.title('Per-marker ICC (one-way batch random effect)', fontsize=12)
plt.show()
plt.savefig(os.path.join(output_dir,"scatter-per-marker-icc-par-vs-aug.png"), dpi=300)



# Inter-marker correlation ------------------------------------------------------

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))

corr_rwd = rwd_clean.corr(method='pearson')
sorted_corr_rwd = corr_rwd.unstack().sort_values(ascending=False).drop_duplicates()
corr_rwd_sel = corr_rwd[abs(corr_rwd)>=0.8]
sns.heatmap(corr_rwd_sel, cmap="coolwarm", vmin=0.8, vmax=1, center=0.8, ax=ax1)
ax1.set_title("MSKCC Sarcoma Samples")
ax1.set(xlabel=None)

corr_par = par_clean.corr(method='pearson')
sorted_corr_par = corr_par.unstack().sort_values(ascending=False).drop_duplicates()
corr_par_sel = corr_par[abs(corr_par)>=0.8]
sns.heatmap(corr_par_sel, cmap="coolwarm", vmin=0.8, vmax=1, center=0.8, ax=ax2)
ax2.set_title("Parametric Simulation")
ax2.set(xlabel=None)

corr_aug = aug_clean.corr(method='pearson')
sorted_corr_aug = corr_aug.unstack().sort_values(ascending=False).drop_duplicates()
corr_aug_sel = corr_aug[abs(corr_aug)>=0.8]
sns.heatmap(corr_aug_sel, cmap="coolwarm", vmin=0.8, vmax=1, center=0.8, ax=ax3)
ax3.set_title("Augmented Data")
ax3.set(xlabel=None)

plt.suptitle("Inter-marker correlation")
plt.tight_layout()
plt.show()
plt.savefig(os.path.join(output_dir,'heatmap-marker-marker-correlation.png'), dpi=300)



# Boxplot ---------------------------------------------------------

n_downsample = 54

## Assemble data
aug_batch_sel = pd.DataFrame()
for i in range(n_batch):
    ## downsampling augmented data
    idx = random.sample([i for i in range(N)], n_downsample) 
    aug_data_sample = aug_data_cap.iloc[[i*N + id for id in idx], P:(2*P)]
    
    ## select 10 true markers (raw count)
    aug_batch_i = 2**aug_data_sample.loc[:, :]#nonzero_genes]
    
    ## add batch label
    aug_batch_i['Batch'] = f'{i+1}' 
    
    aug_batch_sel = pd.concat([aug_batch_sel, aug_batch_i], axis=0)

aug_batch_sel['Type'] = 'Augmented'

par_batch_sel = 2**par_batch
par_batch_sel['Batch'] = np.repeat([f'{i+1}' for i in range(n_batch)], 54)
par_batch_sel['Type'] = 'Parametric'

rwd_batch_sel = 2**rwd_batch.loc[:,:]# nonzero_genes]
rwd_batch_sel['Batch'] = '0'
rwd_batch_sel['Type'] = 'RWD'

# ## RWD vs. Augmented
# box_data = pd.concat([rwd_batch_sel, aug_batch_sel], axis=0)
# box_data['Read Depth (sub)'] = np.log2(box_data.loc[:, nonzero_genes].sum(axis=1)+1)
# box_data['Read Depth'] = np.log2(box_data.drop(columns=['Batch']).sum(axis=1)+1)

## RWD vs. Parametric
box_data = pd.concat([rwd_batch_sel, par_batch_sel, aug_batch_sel], axis=0)
box_data['Read Depth (sub)'] = np.log2(box_data.loc[:, nonzero_genes].sum(axis=1)+1)
box_data['Read Depth'] = np.log2(box_data.drop(columns=['Batch', 'Type']).sum(axis=1)+1)


# ## Density plot
# fig = plt.figure(figsize=(4,5))
# sns.displot( data = box_data, x='Read Depth', hue='Batch', kind='kde')
# plt.show()

## Boxplot 
fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (13.5,6))
sns.boxplot(data=box_data, x='Batch', y="Read Depth", hue='Type', palette=palette,#'Set2', 
            ax=ax1)
sns.boxplot(data=box_data, x='Batch', y="Read Depth (sub)", hue='Type', palette=palette, #'Set2',
            ax=ax2)
ax1.set_title("All 538 markers")
ax2.set_title("10 true markers")
ax2.set_ylabel("")
ax1.legend([], [], frameon=False)
ax2.legend([], [], frameon=False)
fig.legend(ax2.get_legend_handles_labels()[0], ax2.get_legend_handles_labels()[1],
           loc='lower center', ncol=3, frameon=False)
plt.tight_layout()
plt.show()

plt.savefig(
    os.path.join(output_dir, f'boxplot-sample-{n_downsample}-total-count-per-batch-{N}.png'), 
    dpi=300
)
