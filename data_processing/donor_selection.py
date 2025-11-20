#!/usr/bin/env python3

import scanpy as sc
import pandas as pd
import numpy as np

mouse = sc.read_h5ad("data/mouse_dataset.h5ad")

orig = pd.Series({
    "oligodendrocyte": 0.26,
    "neuron": 0.25,
    "astrocyte": 0.09,
    "microglial cell": 0.035,
    "mural cell": 0.004,
    "endothelial cell": 0.002,
    "ependymal cell": 0.002,
    "fibroblast": 0.002
})

lam = 0.7
uniform = pd.Series(1 / len(orig), index=orig.index)
target = (1 - lam) * orig + lam * uniform
target = target / target.sum()

counts = mouse.obs.groupby(["Sample_ID", "cell_type"]).size().unstack(fill_value=0)

selected = []
current = pd.Series(0, index=target.index, dtype=float)

while True:
    best = None
    best_d = None
    for donor in counts.index:
        if donor in selected:
            continue
        cand = current + counts.loc[donor]
        prop = cand / cand.sum()
        prop = prop.reindex(target.index, fill_value=0)
        score = (prop - target).abs().sum()
        if best is None or score < best:
            best = score
            best_d = donor
    if best_d is None or len(selected) >= 50:
        break
    selected.append(best_d)
    current += counts.loc[best_d]

mask = mouse.obs["Sample_ID"].isin(selected)
mouse_sel = mouse[mask].copy()
mouse_sel.write("data/mouse_selected.h5ad", compression="gzip")

print("Selected donors:", len(selected))
print("Final proportions:\n", (current / current.sum()).round(3))


