import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, ConfusionMatrixDisplay
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import homogeneity_score, adjusted_rand_score, normalized_mutual_info_score, f1_score
import umap.umap_ as umap
from sklearn.preprocessing import LabelEncoder

# replace with actual file paths and subsample pairing
mouse = folder_path + "mouse_1.h5ad"
mouse_data = sc.read_h5ad(mouse)
human = folder_path + "/human_test_sample_1.h5ad"
human_data = sc.read_h5ad(human)

def strip_for_knn_inplace(adata, label_col):
    adata.raw = None
    if hasattr(adata, "layers") and adata.layers:
        for k in list(adata.layers.keys()):
            del adata.layers[k]
    for attr in ["obsm", "varm", "obsp", "varp", "uns"]:
        obj = getattr(adata, attr, None)
        if obj is not None:
            obj.clear()
    keep_cols = {label_col, *extra_obs_to_keep}
    cols_to_keep = [c for c in adata.obs.columns if c in keep_cols]
    adata.obs = adata.obs[cols_to_keep].copy()
    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    adata.X = X.astype(np.float32, copy=False)
    gc.collect()
    return adata

label_col = "cell_type"
mouse = strip_for_knn_inplace(mouse_data, label_col)
human = strip_for_knn_inplace(human_data, label_col)

del mouse_data, human_data
gc.collect()

mouse_labels = set(mouse.obs[label_col])
human_labels = set(human.obs[label_col])

X_train = mouse.X 
y_train = mouse.obs[label_col].to_numpy()

X_test  = human.X
y_test  = human.obs[label_col].to_numpy()

pca = PCA(n_components=100, random_state=0)
X_train_pca = pca.fit_transform(X_train)
X_test_pca  = pca.transform(X_test)

knn = KNeighborsClassifier(
    n_neighbors=5,
    weights="distance",
    n_jobs=-1,
    algorithm="auto")

knn.fit(X_train_pca, y_train)
y_pred = knn.predict(X_test_pca)

cm = confusion_matrix(y_test, y_pred)

homogeneity   = homogeneity_score(y_test, y_pred)
ari = adjusted_rand_score(y_test, y_pred)
nmi = normalized_mutual_info_score(y_test, y_pred)
weighted_f1 = f1_score(y_test, y_pred, average="weighted")

results = pd.DataFrame({
    "homogeneity": homogeneity,
    "ari": ari,
    "nmi": nmi,
    "weighted_f1": weighted_f1}).set_index("sample")

pd.write_csv("results/knn_results.csv")

reducer = umap.UMAP(
    n_components=2,
    n_neighbors=15,
    min_dist=0.3,
    metric="euclidean",
    random_state=0,
)

X_test_umap = reducer.fit_transform(X_test_pca)
le = LabelEncoder()
pred_encoded = le.fit_transform(y_pred)

plt.figure(figsize=(8, 6))
sc = plt.scatter(
    X_test_umap[:, 0],
    X_test_umap[:, 1],
    c=pred_encoded,
    s=1,
    alpha=0.7)
plt.title("Human cells UMAP colored by KNN-predicted cell type")
plt.xticks([])
plt.yticks([])
plt.tight_layout()
plt.show()