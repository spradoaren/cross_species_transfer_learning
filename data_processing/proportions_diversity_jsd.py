import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
from IPython.display import display, Markdown
import scanpy.external as sce
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import seaborn as sns
import torch
from torch import nn
from torch.utils.data import TensorDataset, DataLoader

human_samples = {"Human Sample 1": folder_path + "human_test_sample_1.h5ad",
    "Human Sample 2": folder_path + "human_test_sample_2.h5ad",
    "Human Sample 3": folder_path + "human_test_sample_3.h5ad",
    "Human Sample 4": folder_path + "human_test_sample_4.h5ad",
    "Human Sample 5": folder_path + "human_test_sample_5.h5ad"}
mouse_samples = {"Mouse 1": folder_path + "mouse_1.h5ad", 
    "Mouse 2": folder_path + "mouse_2.h5ad",
    "Mouse 3": folder_path + "mouse_3.h5ad"}
all_samples = {**human_samples, **mouse_samples}
celltype_col = "cell_type"

# cell type proportions in each dataset
all_proportions = {}
all_cell_types = set()

for name, path in all_samples.items():
    adata = sc.read_h5ad(path) 
    counts = adata.obs[celltype_col].value_counts()
    total_cells = counts.sum()
    proportions = (counts / total_cells) * 100
    all_proportions[name] = proportions
    all_cell_types.update(proportions.index)

final_df = pd.DataFrame(index=sorted(list(all_cell_types)))
for name, prop_series in all_proportions.items():
    final_df[name] = prop_series.reindex(final_df.index, fill_value=0.0)

final_df = final_df.round(2)

# shannon entropy and jensen-shannon divergence

def get_celltype_counts(adata, all_types=None):
    s = adata.obs[celltype_col].value_counts()
    s = s.sort_index()
    if all_types is not None:
        s = s.reindex(all_types).fillna(0)
    return s
mouse_counts = mouse.obs[celltype_col].value_counts().sort_index()
all_types = mouse_counts.index

def _to_probs(counts, eps=1e-12):
    counts = np.asarray(counts, dtype=float)
    counts = counts + eps
    probs = counts / counts.sum()
    return probs

def shannon_entropy(counts, base=2, normalize=True):
    p = _to_probs(counts)
    log_p = np.log(p) / np.log(base)
    H = -np.sum(p * log_p)

    if not normalize:
        return H

    K = len(p)
    H_max = np.log(K) / np.log(base)
    return H / H_max if H_max > 0 else np.nan

def jensen_shannon_divergence(p_counts, q_counts, base=2):
    p = _to_probs(p_counts)
    q = _to_probs(q_counts)
    m = 0.5 * (p + q)
    def kl(a, b):
        return np.sum(a * (np.log(a) - np.log(b)) / np.log(base))
    return 0.5 * kl(p, m) + 0.5 * kl(q, m)

human_counts = {1: get_celltype_counts(humansample_1, all_types=all_types),
    2: get_celltype_counts(humansample_2, all_types=all_types),
    3: get_celltype_counts(humansample_3, all_types=all_types),
    4: get_celltype_counts(humansample_4, all_types=all_types),
    5: get_celltype_counts(humansample_5, all_types=all_types)}

diversities = {}
jsd_values = {}

for s, counts in human_counts.items():
    diversities[s] = shannon_entropy(counts.values, base=2, normalize=True)
    jsd_values[s] = jensen_shannon_divergence(mouse_counts.values, counts.values, base=2)

summary = pd.DataFrame({
    "sample": list(diversities.keys()),
    "diversity": list(diversities.values()),
    "jsd_mouse_vs_human": [jsd_values[s] for s in diversities.keys()]}).set_index("sample")
