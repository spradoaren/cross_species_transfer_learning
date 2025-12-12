# Cross-Species Transfer Learning for Single-Cell RNA Sequencing

Please see our Final Report file.

This repository contains code for cross-species cell type annotation using transfer learning methods. The project compares different approaches for transferring cell type labels from mouse to human single-cell RNA sequencing (scRNA-seq) data.

## Project Overview

The goal of this project is to develop and evaluate methods for annotating human cell types using models trained on mouse data. This is particularly useful when labeled human data is scarce but abundant mouse data is available. The project includes:

1. **Data Processing Pipeline**: Preprocessing and alignment of mouse and human scRNA-seq datasets
2. **Model Implementations**: Three different transfer learning approaches
3. **Evaluation Metrics**: Evaluation of model performance
4. **Our Final Report**

## Repository Structure

```
cross_species_transfer_learning/
├── data_processing/           Data preprocessing scripts
│   ├── data_retrieval.py      Load and standardize datasets
│   ├── donor_selection.py     Select mouse donors with target distribution
│   ├── ortholog_filtering.py  Filter to shared orthologs
│   ├── proportions_diversity_jsd.py Compute diversity metrics
│   ├── orthologs.txt          Mouse-human ortholog mapping
│   └── README.md              Detailed data processing documentation
├── models/                    Model implementations
│   ├── harmony_mlp.py        Harmony integration + MLP classifier
│   ├── knn.py                KNN classifier with PCA
│   └── seurat_v4.rmd         Seurat v4 label transfer (R)

## Data Requirements

The data files are not included due to size. Place the following in `data/`:

- `mouse_dataset.h5ad` - Murine HypoMap data in .h5ad format from [Cambridge Repository](https://www.repository.cam.ac.uk/items/8f9c3683-29fd-44f3-aad5-7acf5e963a75)
- `human_dataset.h5ad` - Human HypoMap data in .h5ad format from [Cambridge Repository](https://www.repository.cam.ac.uk/items/cad1c61a-e4e5-4443-ad11-92e4f48b3861)
- `orthologs.txt` - Mouse-human 1:1 orthologs from Ensembl BioMart

## Model Descriptions

### 1. Harmony + MLP
- **Integration**: Uses Harmony to integrate mouse and human data in PCA space
- **Classifier**: Multi-layer perceptron (MLP) with dropout regularization
- **Features**: Uses highly variable genes, then PCA + Harmony integration
- **Training**: Trained on mouse data with class-weighted loss function

### 2. KNN Classifier
- **Dimensionality Reduction**: PCA with 100 components
- **Classifier**: K-nearest neighbors (K=5) with distance weighting
- **Features**: Uses raw gene expression data after PCA
- **Visualization**: UMAP plots comparing predicted vs true labels

### 3. Seurat v4 Label Transfer
- **Method**: Seurat's canonical correlation analysis (CCA) based label transfer
- **Integration**: Uses FindTransferAnchors and TransferData functions
- **Features**: Selects integration features from both datasets
- **Visualization**: UMAP plots comparing predicted vs true labels


