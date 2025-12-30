# Quick Start Guide: Fixing the cNMF transfer_to_adata Error

## Problem
You're getting this error:
```python
AttributeError: 'cNMF' object has no attribute 'transfer_to_adata'
```

## Solution
Use the extended cNMF class provided in `cnmf_extensions.py` instead of the standard cNMF class.

## Quick Fix (2 steps)

### Step 1: Update your import
**Before:**
```python
from cnmf import cNMF
```

**After:**
```python
from cnmf_extensions import cNMF
```

### Step 2: Use transfer_to_adata
```python
# Load your AnnData object
import scanpy as sc
adata = sc.read_h5ad('./data/integrated/47.integrated.h5ad')

# Transfer cNMF results to AnnData
adata = cnmf_obj.transfer_to_adata(
    adata, 
    K=8,                    # Replace with your K value
    density_threshold=0.1   # Replace with your density_threshold
)

# That's it! Your cNMF results are now in adata
```

## What Gets Added to Your AnnData

After calling `transfer_to_adata`, your AnnData object will have:

```python
# Cell-level usage scores for each program
adata.obsm['X_cnmf_usage']      # Shape: (n_cells, K)

# Gene-level spectra scores for each program
adata.varm['cnmf_spectra']      # Shape: (n_genes, K)

# Top marker genes for each program
adata.uns['cnmf_top_genes']     # DataFrame

# Analysis metadata
adata.uns['cnmf_top_genes_metadata']  # Dict with K, density_threshold, etc.
```

## Accessing Results

```python
# Get usage scores for all cells
usage = adata.obsm['X_cnmf_usage']
print(f"Usage matrix shape: {usage.shape}")

# Get spectra scores for all genes
spectra = adata.varm['cnmf_spectra']
print(f"Spectra shape: {spectra.shape}")

# Get top genes for Program 1
top_genes = adata.uns['cnmf_top_genes']
print("Top 10 genes for Program 1:")
print(top_genes.iloc[:10, 0])
```

## Visualization Example

```python
import scanpy as sc

# Compute UMAP (if not already done)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Plot the first cNMF program on UMAP
sc.pl.umap(adata, color=adata.obsm['X_cnmf_usage'][:, 0], 
           title='cNMF Program 1')
```

## Common Issues

### Issue 1: K or density_threshold not specified
**Error:**
```
ValueError: K (number of programs) must be specified.
```

**Solution:** Provide the K and density_threshold values you used in your cNMF analysis:
```python
adata = cnmf_obj.transfer_to_adata(adata, K=8, density_threshold=0.1)
```

### Issue 2: Number of cells doesn't match
**Error:**
```
ValueError: Number of cells in cNMF results (1000) does not match number of observations in AnnData (1200).
```

**Solution:** Make sure you're using the same AnnData object that was used for the cNMF analysis. If you filtered cells after cNMF, subset your AnnData to match:
```python
# If cNMF was run on a subset of cells
cell_names_in_cnmf = [...]  # Get this from your cNMF analysis
adata_subset = adata[cell_names_in_cnmf, :]
adata_subset = cnmf_obj.transfer_to_adata(adata_subset, K=8, density_threshold=0.1)
```

### Issue 3: Gene names don't match
**Warning:**
```
⚠ WARNING: Only 500/2000 (25.0%) genes match between cNMF results and AnnData.
```

**Solution:** This usually happens when gene IDs are different (e.g., gene symbols vs Ensembl IDs). Convert your gene names to match:
```python
# If cNMF uses gene symbols but AnnData uses Ensembl IDs
# You'll need to convert one to the other before running transfer_to_adata
```

## Need More Help?

- **Full Documentation:** See `CNMF_EXTENSIONS_README.md`
- **Example Script:** See `example_cnmf_usage.py`
- **Run Tests:** `python3 test_cnmf_extensions.py`

## Requirements

Make sure you have these packages installed:
```bash
pip install cnmf scanpy anndata numpy pandas
```

## Complete Example Workflow

Here's a complete example from start to finish:

```python
import scanpy as sc
from cnmf_extensions import cNMF

# 1. Load your data
adata = sc.read_h5ad('./data/your_data.h5ad')

# 2. Initialize cNMF (if not already done)
cnmf_obj = cNMF(output_dir='./cnmf_output', name='my_analysis')

# 3. Run cNMF analysis (if not already done)
# cnmf_obj.prepare(...)
# cnmf_obj.factorize(...)
# cnmf_obj.combine()
# cnmf_obj.consensus(k=8, density_threshold=0.1)

# 4. Transfer results to AnnData - THIS IS THE NEW PART!
adata = cnmf_obj.transfer_to_adata(
    adata,
    K=8,
    density_threshold=0.1,
    n_top_genes=100
)

# 5. Use the results in your downstream analysis
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color=adata.obsm['X_cnmf_usage'][:, 0])

# 6. Save your updated AnnData
adata.write('./data/your_data_with_cnmf.h5ad')
```

## That's It!

You should now be able to use `transfer_to_adata` successfully. If you encounter any other issues, please check the full documentation in `CNMF_EXTENSIONS_README.md`.
