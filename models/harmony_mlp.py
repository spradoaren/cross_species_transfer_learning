import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import scanpy.external as sce
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix, accuracy_score, adjusted_rand_score, normalized_mutual_info_score, classification_report
import matplotlib.pyplot as plt
import seaborn as sns
import torch
from torch import nn
from torch.utils.data import TensorDataset, DataLoader
from sklearn.utils.class_weight import compute_class_weight

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# replace with actual file paths and subsample pairing
mouse = folder_path + "mouse_1.h5ad"
mouse_data = sc.read_h5ad(mouse)
human = folder_path + "/human_test_sample_1.h5ad"
human_data = sc.read_h5ad(human)

sc.pp.highly_variable_genes(mouse, n_top_genes=1600, flavor="seurat")
sc.pp.highly_variable_genes(human, n_top_genes=1600, flavor="seurat")

hvg_mask = mouse.var["highly_variable"] | human.var["highly_variable"]
mouse_hvg = mouse[:, hvg_mask].copy()
human_hvg = human[:, hvg_mask].copy()

mouse_hvg.obs["species"] = "mouse"
human_hvg.obs["species"] = "human"

combined = ad.concat([mouse_hvg, human_hvg], join="inner")

sc.tl.pca(combined, n_comps=50, svd_solver="arpack")
sce.pp.harmony_integrate(combined, key="species", basis="X_pca")

is_mouse = combined.obs["species"] == "mouse"
is_human = combined.obs["species"] == "human"

X_mouse = np.asarray(combined.obsm["X_pca_harmony"][is_mouse.values, :])
X_human = np.asarray(combined.obsm["X_pca_harmony"][is_human.values, :])

y_mouse = mouse.obs["cell_type"].astype("category")
y_human = human.obs["cell_type"].astype("category")

mouse_label_names = y_mouse.cat.categories
y_mouse_int = y_mouse.cat.codes.to_numpy()

y_mouse_int = y_mouse.cat.codes.to_numpy()

X_train, X_val, y_train, y_val = train_test_split(X_mouse,
    y_mouse_int,
    test_size=0.2,
    random_state=0,
    stratify=y_mouse_int)

X_train_t = torch.tensor(X_train, dtype=torch.float32)
y_train_t = torch.tensor(y_train, dtype=torch.long)

X_val_t   = torch.tensor(X_val,   dtype=torch.float32)
y_val_t   = torch.tensor(y_val,   dtype=torch.long)

X_human_t = torch.tensor(X_human, dtype=torch.float32)

train_ds = TensorDataset(X_train_t, y_train_t)
val_ds   = TensorDataset(X_val_t,   y_val_t)

train_loader = DataLoader(train_ds, batch_size=256, shuffle=True)
val_loader   = DataLoader(val_ds,   batch_size=512, shuffle=False)

input_dim = X_mouse.shape[1]
num_classes = len(mouse_label_names)

class_weights = compute_class_weight(
    class_weight="balanced",
    classes=np.arange(len(mouse_label_names)),
    y=y_mouse_int)

class_weights_t = torch.tensor(class_weights, dtype=torch.float32).to(device)

class MLPClassifier(nn.Module):
    def __init__(self, input_dim, num_classes):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, 512),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(512, 256),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(256, num_classes))

    def forward(self, x):
        return self.net(x)

model = MLPClassifier(input_dim, num_classes).to(device)
criterion = nn.CrossEntropyLoss(weight=class_weights_t)
optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-4)

def evaluate(loader):
    model.eval()
    correct, total = 0, 0
    with torch.no_grad():
        for xb, yb in loader:
            xb, yb = xb.to(device), yb.to(device)
            logits = model(xb)
            preds = logits.argmax(dim=1)
            correct += (preds == yb).sum().item()
            total += yb.size(0)
    return correct / total

n_epochs = 50
for epoch in range(1, n_epochs + 1):
    model.train()
    epoch_loss = 0.0
    for xb, yb in train_loader:
        xb, yb = xb.to(device), yb.to(device)
        optimizer.zero_grad()
        logits = model(xb)
        loss = criterion(logits, yb)
        loss.backward()
        optimizer.step()
        epoch_loss += loss.item() * yb.size(0)
    epoch_loss /= len(train_ds)
    train_acc = evaluate(train_loader)
    val_acc = evaluate(val_loader)

model.eval()
with torch.no_grad():
    logits_human = model(X_human_t.to(device))
    preds_int = logits_human.argmax(dim=1).cpu().numpy()

pred_labels = np.array(mouse_label_names)[preds_int]

human.obs["predicted_celltype_dl"] = pred_labels
human.obs["true_celltype"] = y_human.astype(str).values

np.unique(pred_labels), np.unique(human.obs["true_celltype"])

true_labels = human.obs["true_celltype"].to_numpy()
pred_dl = human.obs["predicted_celltype_dl"].to_numpy()

labels_sorted = np.unique(true_labels)
cm = confusion_matrix(true_labels, pred_dl, labels=labels_sorted)
cm_df = pd.DataFrame(
    cm,
    index=pd.Index(labels_sorted, name="True"),
    columns=pd.Index(labels_sorted, name="Predicted"))
class_rows = [lbl for lbl in report_df.index
              if lbl not in ["accuracy", "macro avg", "weighted avg"]]

cell_counts = pd.Series(true_labels).value_counts()

valid_classes = [c for c in class_rows if c in cell_counts.index]

f1_per_class = report_df.loc[valid_classes, "f1-score"].to_numpy()
weights = cell_counts[valid_classes] / cell_counts.sum()

weighted_f1 = np.sum(f1_per_class * weights.values)
acc = accuracy_score(true_labels, pred_dl)
ari = adjusted_rand_score(true_labels, pred_dl)
nmi = normalized_mutual_info_score(true_labels, pred_dl)

df_results = pd.DataFrame({
    "accuracy": acc,
    "ari": ari,
    "nmi": nmi,
    "weighted_f1": weighted_f1}).set_index("sample")

pd.write_csv("results/harmony_mlp_results.csv")

df_eval = pd.DataFrame({"true": true_labels, "pred": pred_dl})
acc_by_type = (
    df_eval.groupby("true")
           .apply(lambda x: (x["true"] == x["pred"]).mean() * 100)
           .reset_index(name="percent_correct"))

plt.figure(figsize=(8, 4))
sns.barplot(data=acc_by_type, x="true", y="percent_correct")
plt.xticks(rotation=45, ha="right")
plt.ylabel("Percent Correctly Identified (%)")
plt.xlabel("Cell Type")
plt.title("Deep Learning: Accuracy of Cell Type Prediction (mouse → human)")
plt.tight_layout()
plt.show()
