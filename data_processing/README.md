# Data Processing

This folder contains the Python scripts necessary for preprocessing, gene-filtering, and strategically sampling single-cell RNA sequencing (scRNA-seq) hypothalamic datasets for cross-species transfer learning.. The output files from this pipeline serve as the standardized inputs for the downstream model evaluation using the KNN (Python), Seurat (R) and Harmony/MLP (Python) methods described in our analysis.

## Data Requirements

The data files are not included due to size. Place the following in `data/`:

- `mouse_dataset.h5ad` - Murine HypoMap data in .h5ad format from [Cambridge Repository](https://www.repository.cam.ac.uk/items/8f9c3683-29fd-44f3-aad5-7acf5e963a75)
- `human_dataset.h5ad` - Human HypoMap data in .h5ad format from [Cambridge Repository](https://www.repository.cam.ac.uk/items/cad1c61a-e4e5-4443-ad11-92e4f48b3861)
- `orthologs.txt` - Mouse-human 1:1 orthologs from Ensembl BioMart (in repo)

## Scripts

### 1. data_retrieval.py
Loads raw mouse and human datasets and standardizes mouse cell-type labels to match human nomenclature.
- Renames mouse cell type labels to match with human
- Produces standardized datasets in memory (no file output)

### 2. donor_selection.py
Selects a subset of mouse donors whose combined cell-type distribution matches a target mixture between the human distribution and a uniform distribution. Use this to test different potential compositions to have more variety in in-sample cell-type diversity and JSD scores.

- Defines target cell type proportions (lambda = 0.7 balances target distribution with uniformity)
- Iteratively selects mouse donors to minimize deviation from target distribution
- Limits selection to maximum 50 donors (can be expanded if no computational limit)
- `lam = 0.7` - Mixing parameter (0 = exact human distribution, 1 = uniform)
- produces `data/mouse_selected.h5ad`: Selected mouse donors with target cell type distribution

### 3. ortholog_filtering.py
Filters mouse and human data to a shared set of 1:1 ortholog genes and aligns their gene order for downstream analysis.

- Loads ortholog mapping from `orthologs.txt`
- Filters to genes present in both datasets
- Ensures 1:1 ortholog mapping
- Aligns gene order between mouse and human datasets
- Adds ortholog metadata to AnnData objects

Produces:
- `data/mouse_ortho.h5ad` - Mouse data filtered to shared orthologs
- `data/human_ortho.h5ad` - Human data filtered to shared orthologs

### 4. proportions_diversity_jsd.py
Computes cell type proportions, Shannon entropy (diversity), and Jensen-Shannon divergence (JSD) between datasets.

- `shannon_entropy()`: Computes normalized Shannon entropy
- `jensen_shannon_divergence()`: Computes JSD between two distributions
- `get_celltype_counts()`: Extracts cell type counts from AnnData objects

**Note:** The scripts expect data files to be in a `data/` directory relative to where the scripts are run
