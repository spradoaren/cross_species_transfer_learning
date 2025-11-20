#!/usr/bin/env python3

import scanpy as sc
import pandas as pd

mouse = sc.read_h5ad("data/mouse_dataset.h5ad")
human = sc.read_h5ad("data/human_dataset.h5ad")

rename = {
    "Neurons": "neuron",
    "Astrocytes": "astrocyte",
    "Oligodendrocytes": "oligodendrocyte",
    "NG/OPC": "oligodendrocyte",
    "Immune": "microglial cell",
    "Endothelial": "endothelial cell",
    "Tanycytes": "ependymal cell",
    "Ependymal": "ependymal cell",
    "Mural": "mural cell",
    "Fibroblast": "fibroblast"
}

mouse.obs["cell_type"] = mouse.obs["Author_Class_Curated"].replace(rename)

# quick checks
print(mouse.obs["cell_type"].value_counts().sort_index())
print(human.obs["cell_type"].value_counts().sort_index())



