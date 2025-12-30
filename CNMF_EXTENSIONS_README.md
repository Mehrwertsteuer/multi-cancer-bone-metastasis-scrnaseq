# cNMF Extensions for AnnData Integration

This module extends the [cNMF (consensus Non-negative Matrix Factorization)](https://github.com/dylkot/cNMF) library with additional functionality for seamless integration with [AnnData](https://anndata.readthedocs.io/) objects used in single-cell RNA-seq analysis workflows.

## Problem

The standard cNMF library doesn't provide a direct way to transfer cNMF analysis results into AnnData objects. Users typically need to manually load results and add them to their AnnData objects, which can be error-prone.

## Solution

The `cnmf_extensions` module provides an extended `cNMF` class that includes a `transfer_to_adata()` method, making it easy to integrate cNMF results directly into your single-cell analysis pipeline.

## Installation

### Prerequisites

```bash
pip install cnmf scanpy anndata
```

### Using the Extension

Simply place the `cnmf_extensions.py` file in your project directory or Python path.

## Usage

### Basic Usage

```python
from cnmf_extensions import cNMF
import scanpy as sc

# Initialize cNMF (same as before)
cnmf_obj = cNMF(output_dir='./cnmf_output', name='my_analysis')

# ... Run your cNMF analysis pipeline ...
# cnmf_obj.prepare()
# cnmf_obj.factorize()
# cnmf_obj.combine()
# cnmf_obj.consensus()

# Load your AnnData object
adata = sc.read_h5ad('data.h5ad')

# Transfer cNMF results to AnnData - this is the new functionality!
adata = cnmf_obj.transfer_to_adata(
    adata,
    K=8,                    # Number of gene expression programs
    density_threshold=0.1   # Density threshold used in consensus
)

# Access the results
print(adata.obsm['X_cnmf_usage'])      # Cell × Programs usage matrix
print(adata.varm['cnmf_spectra'])      # Gene × Programs spectra scores
print(adata.uns['cnmf_top_genes'])     # Top marker genes per program
```

### Advanced Usage with Custom Keys

```python
# Use custom key names for storing results
adata = cnmf_obj.transfer_to_adata(
    adata,
    K=10,
    density_threshold=0.15,
    n_top_genes=200,                        # Get top 200 genes per program
    norm_usage=True,                        # Normalize usage to sum to 1
    usage_key='X_cnmf_programs',           # Custom key for usage matrix
    spectra_key='cnmf_gene_scores',        # Custom key for gene spectra
    top_genes_key='cnmf_markers'           # Custom key for top genes
)
```

### Integration with Scanpy Workflows

```python
import scanpy as sc
from cnmf_extensions import cNMF

# Load and preprocess data
adata = sc.read_h5ad('raw_data.h5ad')
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Run cNMF analysis
cnmf_obj = cNMF(output_dir='./cnmf_out', name='sample1')
# ... complete cNMF steps ...

# Transfer results to AnnData
adata = cnmf_obj.transfer_to_adata(adata, K=8, density_threshold=0.1)

# Visualize cNMF programs with UMAP
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Plot usage of each program
sc.pl.umap(adata, color=['cNMF_1', 'cNMF_2'], use_raw=False)

# Or directly plot from obsm
import matplotlib.pyplot as plt
fig, axes = plt.subplots(2, 4, figsize=(16, 8))
for i, ax in enumerate(axes.flat):
    if i < adata.obsm['X_cnmf_usage'].shape[1]:
        sc.pl.umap(adata, 
                   color=adata.obsm['X_cnmf_usage'][:, i],
                   ax=ax, 
                   show=False,
                   title=f'Program {i+1}')
plt.tight_layout()
plt.show()
```

## What Gets Added to AnnData

The `transfer_to_adata()` method adds the following to your AnnData object:

1. **`adata.obsm['X_cnmf_usage']`** (default key)
   - Shape: (n_cells, K)
   - Cell-level usage scores for each gene expression program
   - Each row represents a cell, each column represents a program

2. **`adata.varm['cnmf_spectra']`** (default key)
   - Shape: (n_genes, K)
   - Gene-level spectra scores for each program
   - Each row represents a gene, each column represents a program
   - Higher scores indicate stronger association with the program

3. **`adata.uns['cnmf_top_genes']`** (default key)
   - DataFrame with top marker genes for each program
   - Useful for biological interpretation of programs

4. **`adata.uns['cnmf_top_genes_metadata']`** (automatically added)
   - Dictionary containing analysis parameters:
     - `K`: Number of programs
     - `density_threshold`: Threshold used
     - `n_top_genes`: Number of top genes retrieved
     - `norm_usage`: Whether usage was normalized
     - `n_common_genes`: Number of overlapping genes

5. **`adata.uns['X_cnmf_usage_names']`** (automatically added)
   - List of program names (e.g., ['cNMF_1', 'cNMF_2', ...])

## Method Parameters

### `transfer_to_adata()`

**Required Parameters:**
- `adata` (AnnData): The AnnData object to add results to
- `K` (int): Number of gene expression programs
- `density_threshold` (float): Density threshold for filtering outliers

**Optional Parameters:**
- `n_top_genes` (int, default=100): Number of top genes to return per program
- `norm_usage` (bool, default=True): Normalize usage scores to sum to 1
- `usage_key` (str, default='X_cnmf_usage'): Key for storing usage matrix
- `spectra_key` (str, default='cnmf_spectra'): Key for storing spectra scores
- `top_genes_key` (str, default='cnmf_top_genes'): Key for storing top genes

**Returns:**
- Modified AnnData object with cNMF results added

**Raises:**
- `ValueError`: If required parameters are missing or dimensions don't match
- `RuntimeError`: If cNMF results cannot be loaded

## Error Handling

The method includes comprehensive error checking:

```python
try:
    adata = cnmf_obj.transfer_to_adata(adata, K=8, density_threshold=0.1)
except ValueError as e:
    print(f"Parameter error: {e}")
except RuntimeError as e:
    print(f"Could not load cNMF results: {e}")
```

Common errors:
- **"K must be specified"**: You need to provide the K value used in your analysis
- **"density_threshold must be specified"**: Provide the threshold from consensus step
- **"Number of cells in cNMF results does not match"**: Make sure cNMF was run on the same dataset
- **"No common genes found"**: Gene names in cNMF results don't match AnnData gene names

## Compatibility

- **Python**: 3.7+
- **cNMF**: 1.4+
- **AnnData**: 0.7+
- **Scanpy**: 1.7+ (optional, for visualization)

## Example Workflow

Here's a complete example from data to visualization:

```python
import scanpy as sc
import pandas as pd
from cnmf_extensions import cNMF

# 1. Prepare data
adata = sc.read_10x_mtx('./data/filtered_feature_bc_matrix/')
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# 2. Initialize and run cNMF
cnmf_obj = cNMF(output_dir='./cnmf_output', name='example')

# Prepare: save counts for cNMF
cnmf_obj.prepare(counts_fn=adata, 
                 components=list(range(5,16)),  # Test K from 5 to 15
                 n_iter=100,
                 seed=14)

# Factorize: run NMF
cnmf_obj.factorize(worker_i=0, total_workers=4)

# Combine results
cnmf_obj.combine()

# Run consensus
cnmf_obj.consensus(k=10, density_threshold=0.10)

# 3. Transfer results to AnnData
adata = cnmf_obj.transfer_to_adata(adata, K=10, density_threshold=0.10)

# 4. Visualize
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Plot individual programs
for i in range(10):
    sc.pl.umap(adata, 
               color=adata.obsm['X_cnmf_usage'][:, i],
               title=f'cNMF Program {i+1}',
               save=f'_program_{i+1}.png')

# 5. Analyze top genes
top_genes = adata.uns['cnmf_top_genes']
print("Top genes for Program 1:")
print(top_genes.iloc[:, 0])
```

## Differences from Base cNMF

This extended class:
1. ✅ Inherits all original cNMF functionality
2. ✅ Adds `transfer_to_adata()` method
3. ✅ Automatically handles gene name matching
4. ✅ Provides detailed error messages
5. ✅ Includes comprehensive documentation
6. ✅ No breaking changes to existing code

You can use it as a drop-in replacement for the standard cNMF class!

## Contributing

If you find bugs or have feature requests, please open an issue or submit a pull request.

## License

This extension follows the same license as the original cNMF project.

## Citation

If you use this extension, please cite the original cNMF paper:

> Kotliar, D. et al. Identifying gene expression programs of cell-type identity and cellular activity with single-cell RNA-Seq. eLife (2019). https://doi.org/10.7554/eLife.43803

## References

- [cNMF GitHub Repository](https://github.com/dylkot/cNMF)
- [AnnData Documentation](https://anndata.readthedocs.io/)
- [Scanpy Documentation](https://scanpy.readthedocs.io/)
