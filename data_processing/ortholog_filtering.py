#!/usr/bin/env python3

import pandas as pd
import scanpy as sc
import numpy as np
import os

DATA = "data"

mouse = sc.read_h5ad(os.path.join(DATA, "mouse_selected.h5ad"))
human = sc.read_h5ad(os.path.join(DATA, "human_dataset.h5ad"))
ortho = pd.read_csv(os.path.join(DATA, "orthologs.txt"), sep="\t")

mouse_genes = pd.Index(mouse.var_names.astype(str))
human_raw = pd.Index(human.var_names.astype(str))
human_ens = human_raw.str.split(".", 1).str[0]
human.var["ensembl_novers"] = human_ens

ortho = ortho.rename(columns={
    "Gene stable ID": "mouse_ensembl",
    "Gene name": "mouse_symbol",
    "Human gene stable ID": "human_ensembl",
    "Human gene name": "human_symbol"
}).dropna(subset=["mouse_symbol", "human_ensembl"]).drop_duplicates()

ortho = ortho[
    ortho["mouse_symbol"].isin(mouse_genes) &
    ortho["human_ensembl"].isin(human_ens)
].drop_duplicates(subset=["human_ensembl"])

shared = ortho["human_ensembl"].values
human_order = [human_ens.get_loc(e) for e in shared]

mouse_order = (
    ortho.set_index("human_ensembl")
    .loc[shared]["mouse_symbol"].values
)

mouse_ortho = mouse[:, mouse_order].copy()
human_ortho = human[:, human_order].copy()

mouse_ortho.var["human_ensembl"] = shared
human_ortho.var["ensembl_novers"] = shared

mouse_ortho.write(os.path.join(DATA, "mouse_ortho.h5ad"), compression="gzip")
human_ortho.write(os.path.join(DATA, "human_ortho.h5ad"), compression="gzip")
