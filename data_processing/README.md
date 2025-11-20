

# Data Processing

The data files are not included due to size. Place the following in `data/`:

- `mouse_dataset.h5ad` (Murine HypoMap data in .h5ad format from https://www.repository.cam.ac.uk/items/8f9c3683-29fd-44f3-aad5-7acf5e963a75)
- `human_dataset.h5ad` (Human HypoMap data in .h5ad format from https://www.repository.cam.ac.uk/items/cad1c61a-e4e5-4443-ad11-92e4f48b3861)
- `orthologs.txt` (mouse–human 1:1 orthologs from Ensembl BioMart)

## Scripts

### 1. data_retrieval.py
Loads raw datasets and standardizes mouse cell-type labels.

### 2. donor_selection.py
Selects a subset of mouse donors whose combined cell-type distribution matches a
target mixture between the human distribution and a uniform distribution.
Outputs `mouse_selected.h5ad`.

### 3. ortholog_filtering.py
Filters mouse and human data to a shared set of 1:1 ortholog genes and aligns
their gene order. Produces:

- `mouse_ortho.h5ad`
- `human_ortho.h5ad`

## Running

python scripts/data_retrieval.py
python scripts/donor_selection.py
python scripts/ortholog_filtering.py

The resulting aligned files can be used for downstream models 
