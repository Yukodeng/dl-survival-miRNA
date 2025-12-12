import os
import re
import random
import numpy as np
import pandas as pd
import umap
from sklearn.preprocessing import StandardScaler
from scipy.stats import iqr
import matplotlib.pyplot as plt
import seaborn as sns



# Load Data --------------------------------------------------------

N = 10000
P = 538 
n_batch = 10
rwd_max = 22.516867
rwd_sd = 4.039547
gnames = pd.read_csv(os.path.join("data", "augmentation_miRNA_names_538.csv"))['gene'].tolist()      # gene names
gindex = pd.read_csv(os.path.join("data", "augmentation_miRNA_names_538.csv"))['index_py'].tolist()  # gene index
gnames_no1 = [g[:-2] for g in gnames]


# * True markers ----------------------------------
nonzero_genes = ["hsa.miR.1277.3p.1", "hsa.miR.1277.5p.1", "hsa.miR.133a.2..1", "hsa.miR.144..1" , 
                 "hsa.miR.148b..1" ,  "hsa.miR.181c..1" ,  "hsa.miR.204.1", "hsa.miR.33a.1",
                 "hsa.miR.598.1", "hsa.miR.887.1" ]


# * Benchmark data ---------------------------------

rwd_clean = pd.read_csv("raw-data/MSKCC_sarcoma_benchmark_data.csv")
rwd_batch = pd.read_csv("raw-data/MSKCC_sarcoma_test_data.csv")
rwd_clean = np.log2(rwd_clean.loc[:, gnames_no1] + 1)
rwd_batch = np.log2(rwd_batch.loc[:, gnames_no1] + 1)


# * Parametric simulation --------------------------

## [added 12/8] for parametrically simulated data with all 1033 columns
par_clean = pd.read_csv(os.path.join("raw-data", "with_without_batch_parametric.csv")).iloc[:, :1033]
par_batch = pd.read_csv(os.path.join("raw-data", "with_without_batch_parametric.csv")).iloc[:, 1033:2066]
batch_id  = pd.read_csv(os.path.join("raw-data", "with_without_batch_parametric.csv")).iloc[:, -1]
par_clean = np.log2(par_clean.iloc[:, gindex] + 1)
par_batch = np.log2(par_batch.iloc[:, gindex] + 1)

pd.concat([
    par_clean.iloc[:, gindex],
    par_batch.iloc[:, gindex],
    batch_id
], axis=1).to_csv(os.path.join("raw-data", "with_without_batch_parametric_538_markers.csv"),index=False)

## [NOTE] use below code if pre-filtered to have 538 markers
# par = pd.read_csv(os.path.join("raw-data", "with_without_batch_parametric.csv"))
# par_clean = np.log2(par.iloc[:54, P:(2*P)] + 1)
# par_batch = np.log2(par.iloc[:, :P]        + 1)

par_clean.columns = par_batch.columns = gnames


# * Augmented data ---------------------------------

## Augmented data w/ and w/out batch (APR112025)
aug_data = pd.DataFrame()

for i in range(n_batch):
    # read in each batch
    dat_sub = pd.read_csv(
        os.path.join("raw-data", "with_without_batch_maf", f"batch-{i+1}_epochES_maf_generated.csv"),
        header=None
    ).iloc[:N, :(P*2)]
    # merge batches
    dat_sub.columns = gnames + gnames
    aug_data = pd.concat([aug_data, dat_sub], axis=0)

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

aug_data_log = np.log2(aug_data + 1)
aug_data_cap = pd.DataFrame(np.stack(
    aug_data_log.apply(cap_and_add_noise, cap=rwd_max, sd=rwd_sd, axis=1)
))

aug_data_cap.columns = gnames + gnames


## Set color palette
palette = { # sns.color_palette("husl", 2)
    "RWD":       "#e28743", # orange #"#33aaff", # Blue 
    "Augmented": "#2596be", # blue   
    "Parametric":"#2a8000"  # green  #ff75dd"  # pink
}


# Median vs. IQR plot ----------------------------------------------------------------

# * MSKCC vs. Augmented ------------------

aug_batch = aug_data_cap.iloc[:, :P]
aug_clean = aug_data_cap.iloc[:, P:]
aug_batch.columns = aug_clean.columns = gnames # rename gene ids in batch-present data

per_marker_clean = pd.DataFrame({"Type": [], 
                                 "median": [],
                                 "IQR": []})
per_marker_clean['Marker'] = list(gnames + gnames)
per_marker_clean['Type']   = list(np.repeat('RWD', P)) + list(np.repeat('Augmented', P))
per_marker_clean['median'] = list(rwd_clean.apply(np.median)) + \
    list(aug_clean.apply(np.median))
per_marker_clean['IQR']    = list(rwd_clean.apply(iqr)) + \
    list(aug_clean.apply(iqr))
per_marker_clean['true'] = ['True' if g in nonzero_genes else 'None' for g in per_marker_clean.Marker]


per_marker_batch = pd.DataFrame({"Type": [], 
                                "median": [],
                                "IQR": []})
per_marker_batch['Marker'] = list(gnames + gnames)
per_marker_batch['Type']   = list(np.repeat('RWD', P)) + list(np.repeat('Augmented', P))
per_marker_batch['median'] = list(rwd_batch.apply(np.median)) + \
    list(aug_batch.apply(np.median))
per_marker_batch['IQR']    = list(rwd_batch.apply(iqr)) + \
    list(aug_batch.apply(iqr))
per_marker_batch['true'] = ['True' if g in nonzero_genes else 'None' for g in per_marker_batch.Marker]


## Scatter Plot

fig, (ax1, ax2) = plt.subplots(1,2, figsize=(11,5))

## Clean
sns.scatterplot(data=per_marker_clean,
                x="median", y="IQR", hue="Type",
                palette=palette, alpha=0.15, marker='o', s=16, ax=ax1)
sns.scatterplot(data=per_marker_clean.loc[per_marker_clean['true']=='True',:],
                x="median", y="IQR", hue="Type", 
                palette=palette, alpha=0.95, marker='o', s=18, ax=ax1)
ax1.legend([], [], frameon=False)
ax1.set_title("Median vs. IQR (Clean)")

## Batch
sns.scatterplot(data=per_marker_batch,
                x="median", y="IQR", hue="Type",
                palette=palette, alpha=0.15, marker='o', s=16, ax=ax2)
sns.scatterplot(data=per_marker_batch.loc[per_marker_batch['true']=='True',:],
                x="median", y="IQR", hue="Type", 
                palette=palette, alpha=0.95, marker='o', s=18, ax=ax2)
sns.move_legend(ax2, 'right', bbox_to_anchor=(1.4,.5), ncol=1, title=None, frameon=False)
ax2.set_title("Median vs. IQR (Batch)")

plt.xlim(0,20)
plt.ylim(0,20)
plt.tight_layout()
plt.show()

# save plot 
plt.savefig(os.path.join("results", "plots", "median-IQR-downsample.png"),
            dpi=300, bbox_inches='tight')

# data_true = pd.concat([
#  pd.DataFrame({
#     'Median': rwd_batch.loc[:,nonzero_genes].apply(np.median).tolist(),
#     'IQR': rwd_batch.loc[:, nonzero_genes].apply(iqr).tolist(),
#     'Type': 'RWD'}),
#  pd.DataFrame({
#     'Median': aug_batch.loc[:,nonzero_genes].apply(np.median).tolist(),
#     'IQR': aug_batch.loc[:, nonzero_genes].apply(iqr).tolist(),
#     'Type': 'Augmented' })
#  ])
# sns.scatterplot(x='Median', y='IQR', hue='Type', data=data_true, 
#                 palette={"RWD (true)":"#ff8d36", "Augmented (true)": "#01bbff"},
#                 alpha=0.95, marker='o', s=16, ax=ax2)


# * MSKCC vs. Parametric ----------------------

per_marker_batch = pd.DataFrame({"Type": [], 
                                 "median": [],
                                 "IQR": []})
per_marker_batch['Marker'] = gnames + gnames
per_marker_batch['Type']   = np.repeat('RWD', P).tolist() + np.repeat('Parametric', P).tolist()
per_marker_batch['median'] = rwd_batch.apply(np.median).tolist() + \
    par_batch.apply(np.median).tolist()
per_marker_batch['IQR']    = rwd_batch.apply(iqr).tolist() + \
    par_batch.apply(iqr).tolist()
per_marker_batch['true'] = ['True' if g in nonzero_genes else 'None' for g in per_marker_batch.Marker]

## Plot (batch)

fig, ax1 = plt.subplots(1,1, figsize=(7,5))
sns.scatterplot(data=per_marker_batch,
                x="median", y="IQR", hue="Type",
                palette=palette, alpha=0.15, marker='o', s=15, ax=ax1)
sns.scatterplot(data=per_marker_batch.loc[per_marker_batch['true']=='True',:],
                x="median", y="IQR", hue="Type", 
                palette=palette, alpha=1, marker='X', s=25, ax=ax1)

ax1.legend(ax1.get_legend_handles_labels()[0][:2], 
           ax1.get_legend_handles_labels()[1][:2], frameon=False)
ax1.set_title("Median vs. IQR (batch)")
plt.show()



# Batch Effect Plots -------------------------------------------------------

rwd_batch_obs = (rwd_batch - rwd_clean) #/ 2
par_batch_obs = (par_batch - par_clean)
aug_batch_obs = (aug_batch - aug_clean) #/ 2


## 1) Histogram distribution

## 
fig, ax = plt.subplots(1,1)
ax.scatter(x=rwd_clean.apply(np.std), y=rwd_batch_obs.apply(np.median),
           color='blue', edgecolors='white', alpha=0.8)
plt.xlabel('stdev')
plt.ylabel('median batch')
plt.title("MSKCC Sarcoma Samples")
plt.show()
plt.savefig(os.path.join("results", "plots", "mskcc-stdev-vs-batch-median.png"), 
            dpi=300, bbox_inches='tight')


## Parametric vs. MSKCC

fig, axes = plt.subplots(3,4, figsize=(15,9))
c = 0
for i in range(3):
    for j in range(4):
        if c < 10:
            axes[i,j].hist((par_batch_obs
                            .iloc[c*54:(c*54+54),:]
                            .apply(np.median)), 
                    bins=30, alpha=0.7, color='#e28743', edgecolor='black', label='Parametric')
            axes[i,j].hist(rwd_batch_obs.apply(np.median),
                    bins=30, alpha=0.7, color='#2596be', edgecolor='black', label='RWD')
            axes[i,j].set_title(f"Batch {c+1}")
            axes[i,j].set_xlim([-5,4])
            axes[i,j].set_ylim([0,90])
        else:
            axes[i,j].set_axis_off()
        c += 1
handles, lables = axes[0,0].get_legend_handles_labels()
fig.legend(handles, lables, loc="lower right", bbox_to_anchor=(0.65, 0.25))
plt.suptitle("Marker-specific median (log2) batch effects: RWD vs. Parametric")
plt.tight_layout()
plt.show()
 


## Parametric vs. Aug

fig, axes = plt.subplots(3,4, figsize=(15,9))
c = 0
for i in range(3):
    for j in range(4):
        if c < 10:
            axes[i,j].hist((par_batch_obs
                            .iloc[c*54:(c*54+54),:]
                            .apply(np.median)), 
                    bins=20, alpha=0.7, color='#e28743', edgecolor='black', label='RWD')
            # axes[i,j].hist(aug_batch_obs.iloc[c*54:(c*54+54),:].apply(np.median),
            #         bins=30, alpha=0.7, color='#2596be', edgecolor='black', label='Augmented')
            axes[i,j].set_title(f"Batch {c+1}")
            axes[i,j].set_xlim([-4,3])
        else:
            axes[i,j].set_axis_off()
        c += 1
handles, lables = axes[0,0].get_legend_handles_labels()
fig.legend(handles, lables, loc="lower right", bbox_to_anchor=(0.65, 0.25))
plt.suptitle("Distribution of marker-specific batch effects (log2): Input vs. Augmented")
plt.tight_layout()
plt.show()
 
plt.savefig(os.path.join("results", "plots", "per-batch-marker-median-effects-rwd-vs-aug.png"), 
            dpi=300, bbox_inches='tight')


## MSKCC vs. Aug
fig, ax = plt.subplots(1,1)
ax.hist(rwd_batch_obs.apply(np.median),
         bins=30, alpha=0.7, color='#e28743', edgecolor='black', label='RWD')
# ax.hist(aug_batch_obs.apply(np.median),
#          bins=30, alpha=0.7, color='#2596be', edgecolor='black', label='Augmented')
ax.set_title("Distribution of marker-specific batch effects (log2):\nMSKCC vs. Augmented")
ax.legend()
plt.show()

plt.savefig(os.path.join("results", "plots", "per-marker-median-effects-mskcc-vs-aug.png"), 
            dpi=300, bbox_inches='tight')



## 2) Batch-Marker Stdev correlation
   
par_batch_10_obs = par_batch_obs.iloc[ [c*54 for c in range(10)],: ]
par_batch_std = par_batch_10_obs.apply(np.std).tolist()
rwd_batch_std = rwd_batch_obs.apply(np.std).tolist()

fig = plt.figure(figsize=(6,5))
data = pd.DataFrame({
    'marker std': rwd_clean.apply(np.std).tolist(),
    'par batch std': par_batch_std,
    'rwd batch std': rwd_batch_std
    })
sns.scatterplot(data=data, x='rwd batch std', y='par batch std')
plt.title("Marker-specific batch stdev: RWD vs. Parametric")
plt.show()

## 2) Inter-marker batch correlaiton ----------


fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15,5))

## RWD
rwd_clean_true = rwd_clean.loc[:, [g[:-2] for g in nonzero_genes]]
rwd_batch_true = rwd_batch.loc[:, [g[:-2] for g in nonzero_genes]]
rwd_clean_true.columns = rwd_batch_true.columns = nonzero_genes # rename gene ids
rwd_batch_obs_sel = (rwd_batch_true - rwd_clean_true) #/ 2

corr_batch_rwd = rwd_batch_obs_sel.corr(method='pearson')
sns.heatmap(corr_batch_rwd, cmap="coolwarm", center=0, vmin=-0.2, vmax=1, ax=ax1)
ax1.set_title("MSKCC per-marker batch correlatioon")


## parametric
par_clean = np.log2(par.iloc[:, P:(2*P)] + 1)
par_batch = np.log2(par.iloc[:, :P]      + 1)
par_clean.columns = par_batch.columns = gnames
par_clean_true = par_clean.loc[:,nonzero_genes]
par_batch_true = par_batch.loc[:,  nonzero_genes]
par_batch_obs_sel = (par_batch_true - par_clean_true) #/ 2

corr_batch_par = par_batch_obs_sel.corr(method='pearson')
sns.heatmap(corr_batch_par, cmap="coolwarm", center=0, vmin=-0.2, vmax=1, ax=ax2)
ax2.set_title("Input dataset per-marker batch correlatioon")


## Augmented
aug_batch_true = aug_batch.loc[:, nonzero_genes]
aug_clean_true = aug_clean.loc[:, nonzero_genes]
### NOTE: extract batch on the natural log scale
### (optional) scale batch effects by 1/2
aug_batch_obs_sel = (aug_batch_true - aug_clean_true) #/ 2

corr_batch_aug = aug_batch_obs_sel.corr(method='pearson')
sns.heatmap(corr_batch_aug, cmap="coolwarm", center=0, vmin=-0.2, vmax=1, ax=ax3)
ax3.set_title("Augmented per-marker batch correlatioon")

plt.tight_layout()
plt.show()



## 3) Marker-marker correlation

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
plt.tight_layout()
plt.show()

corr_aug = aug_clean.corr(method='pearson')
sorted_corr_aug = corr_aug.unstack().sort_values(ascending=False).drop_duplicates()
corr_aug_sel = corr_aug[abs(corr_aug)>=0.8]
sns.heatmap(corr_aug_sel, cmap="coolwarm", vmin=0.8, vmax=1, center=0.8, ax=ax2)
ax2.set_title("Augmented Data")
ax2.set(xlabel=None)
plt.tight_layout()
plt.show()


# UMAP -------------------------------------------------------------------


## Downsampling augmented data (N=540)
random.seed(123)

aug_data_sel = pd.DataFrame()
for i in range(10):
    idx = random.sample([i for i in range(N)], 54) 
    aug_sub = aug_data_cap.iloc[[i*N + id for id in idx], :]
    aug_data_sel = pd.concat([aug_data_sel, aug_sub], axis=0)

## Separate clean and batch-contaminated data
aug_batch = aug_data_sel.iloc[:, :P]
aug_clean = aug_data_sel.iloc[:, P:]
aug_batch.columns = aug_clean.columns = gnames # rename gene ids in batch-present data

rwd_clean["Type"] = "RWD"
aug_clean["Type"] = "Augmented"
clean = pd.concat([rwd_clean, aug_clean]).reset_index(drop=True)

rwd_batch["Type"] = "RWD"
aug_batch["Type"] = "Augmented"
batch = pd.concat([rwd_batch, aug_batch]).reset_index(drop=True)


## Create umap projection
umap_obj = umap.UMAP(n_components=2, random_state=42)

clean_scaled = StandardScaler().fit_transform(clean.drop(columns='Type')) 
embedding_clean = umap_obj.fit_transform(clean_scaled)
embedding_clean = pd.concat([
    pd.DataFrame(embedding_clean, columns = ['UMAP1', "UMAP2"]),
    clean[['Type']]
    ], axis=1
 )


batch_scaled = StandardScaler().fit_transform(batch.drop(columns='Type')) 
embedding_batch = umap_obj.fit_transform(batch_scaled)
embedding_batch = pd.concat([
    pd.DataFrame(embedding_batch, columns = ['UMAP1', "UMAP2"]),
    batch[['Type']]
    ], axis=1
 )


## UMAP plots

fig, (ax1, ax2) = plt.subplots(1,2, figsize=(11,5))

## UMAP (Clean)
sns.scatterplot(data=embedding_clean[embedding_clean['Type']=='Augmented'], 
                x="UMAP1", y="UMAP2", hue="Type", 
                palette=palette, alpha=0.75, marker='o', s=16, ax=ax1)
sns.scatterplot(data=embedding_clean[embedding_clean['Type']=='RWD'],
                x="UMAP1", y="UMAP2", hue="Type", 
                palette=palette, alpha=0.85, marker='o', s=20, ax=ax1)
ax1.legend([], [], frameon=False)
ax1.set_title("UMAP Projection (Clean)")

## UMAP (Batch)
sns.scatterplot(data=embedding_batch, #[embedding_batch['Type']=='Augmented'], 
                x="UMAP1", y="UMAP2", hue="Type", 
                palette=palette, alpha=0.75, marker='o', s=18, ax=ax2)
# sns.scatterplot(data=embedding_batch[embedding_batch['Type']=='RWD'], 
#                 x="UMAP1", y="UMAP2", hue="Type", 
#                 palette=palette, alpha=0.8, marker='o', s=15, ax=ax2)

sns.move_legend(ax2, 'right', bbox_to_anchor=(1.35,.5), ncol=1, title=None, frameon=False)
ax2.set_title("UMAP Projection (Batch)")

plt.tight_layout()
plt.show()

# save plot 
plt.savefig(os.path.join("results", "plots", "umap-downsample-default.png"),
            dpi=300, bbox_inches='tight')



# Mean vs. Stdev plot ----------------------------------------------------------------

aug_batch = aug_data_cap.iloc[:, :P]
aug_clean = aug_data_cap.iloc[:, P:]
aug_batch.columns = aug_clean.columns = gnames # rename gene ids in batch-present data


per_marker_clean = pd.DataFrame({"Type": [], 
                                "mean": [],
                                "sd": []})
per_marker_clean['Type'] = list(np.repeat('RWD', P)) + list(np.repeat('Augmented', P))
per_marker_clean['mean'] = list(rwd_clean.iloc[:,:-1].apply(np.mean)) + \
    list(aug_clean.apply(np.mean))
per_marker_clean['sd']   = list(rwd_clean.iloc[:,:-1].apply(np.std)) + \
    list(aug_clean.apply(np.std))



per_marker_batch = pd.DataFrame({"Type": [], 
                              "mean": [],
                              "sd": []})
per_marker_batch['Type'] = list(np.repeat('RWD', P)) + list(np.repeat('Augmented', P))
per_marker_batch['mean'] = list(rwd_batch.iloc[:,:-1].apply(np.mean)) + \
    list(aug_batch.apply(np.mean))
per_marker_batch['sd']   = list(rwd_batch.iloc[:,:-1].apply(np.std)) + \
    list(aug_batch.apply(np.std))


## Scatter Plot ------------------

fig, (ax1, ax2) = plt.subplots(1,2, figsize=(11,5))

## Clean
sns.scatterplot(data=per_marker_clean[per_marker_clean['Type']=='Augmented'], 
                x="mean", y="sd", hue="Type", 
                palette=palette, alpha=0.75, marker='o', s=16, ax=ax1)
sns.scatterplot(data=per_marker_clean[per_marker_clean['Type']=='RWD'],
                x="mean", y="sd", hue="Type", 
                palette=palette, alpha=0.85, marker='o', s=20, ax=ax1)
ax1.legend([], [], frameon=False)
ax1.set_title("Mean vs. Stdev (Clean)")

## Batch
sns.scatterplot(data=per_marker_batch,
                x="mean", y="sd", hue="Type", 
                palette=palette, alpha=0.75, marker='o', s=18, ax=ax2)

sns.move_legend(ax2, 'right', bbox_to_anchor=(1.35,.5), ncol=1, title=None, frameon=False)
ax2.set_title("Mean vs. Stdev (Batch)")

plt.tight_layout()
plt.show()

# save plot 
plt.savefig(os.path.join("results", "plots", "mean-stdev-downsample.png"),
            dpi=300, bbox_inches='tight')


