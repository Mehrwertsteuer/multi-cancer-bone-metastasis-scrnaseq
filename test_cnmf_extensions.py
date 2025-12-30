"""
Integration test for cnmf_extensions module.

This test creates mock data to demonstrate the transfer_to_adata functionality
without requiring actual cNMF analysis results.
"""

import numpy as np
import pandas as pd
import tempfile
import os
from pathlib import Path


def create_mock_cnmf_results(output_dir, name, K=5, n_cells=100, n_genes=50):
    """
    Create mock cNMF result files for testing purposes.
    
    This mimics the structure of real cNMF output files.
    """
    # Create output directory
    results_dir = Path(output_dir) / name
    results_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate mock data
    np.random.seed(42)
    
    # Usage matrix (cells x programs)
    usage = pd.DataFrame(
        np.random.rand(n_cells, K),
        index=[f'cell_{i}' for i in range(n_cells)],
        columns=[f'program_{i}' for i in range(K)]
    )
    
    # Normalize usage to sum to 1 per cell
    usage = usage.div(usage.sum(axis=1), axis=0)
    
    # Gene spectra scores (programs x genes)
    genes = [f'gene_{i}' for i in range(n_genes)]
    spectra_scores = pd.DataFrame(
        np.random.randn(K, n_genes) * 2,  # z-scores
        index=[f'program_{i}' for i in range(K)],
        columns=genes
    )
    
    # Gene spectra TPM (programs x genes)
    spectra_tpm = pd.DataFrame(
        np.random.rand(K, n_genes) * 100,
        index=[f'program_{i}' for i in range(K)],
        columns=genes
    )
    
    # Top genes (sorted by score for each program)
    top_genes_data = {}
    for i in range(K):
        program_scores = spectra_scores.iloc[i].sort_values(ascending=False)
        top_genes_data[f'program_{i}'] = program_scores.index[:20].tolist()
    
    top_genes = pd.DataFrame(top_genes_data)
    
    # Save files in the format cNMF expects
    # cNMF uses k_%d and dt_%s format where:
    # - k_%d is the number of programs (e.g., k_5)
    # - dt_%s is the density threshold with dots replaced by underscores (e.g., 0.1 becomes dt_0_1)
    k_str = f'k_{K}'
    density_str = 'dt_0_1'  # density threshold 0.1 with dots replaced by underscores
    
    # Save files
    usage_file = results_dir / f'{name}.usages.{k_str}.{density_str}.consensus.txt'
    spectra_scores_file = results_dir / f'{name}.gene_spectra_score.{k_str}.{density_str}.txt'
    spectra_tpm_file = results_dir / f'{name}.gene_spectra_tpm.{k_str}.{density_str}.txt'
    top_genes_file = results_dir / f'{name}.top_genes.{k_str}.{density_str}.txt'
    
    usage.to_csv(usage_file, sep='\t')
    spectra_scores.T.to_csv(spectra_scores_file, sep='\t')  # Transpose to genes x programs
    spectra_tpm.T.to_csv(spectra_tpm_file, sep='\t')
    top_genes.to_csv(top_genes_file, sep='\t', index=False)
    
    return {
        'usage': usage,
        'spectra_scores': spectra_scores,
        'genes': genes,
        'cells': usage.index.tolist()
    }


def create_mock_adata(n_cells=100, n_genes=50):
    """Create a mock AnnData object for testing."""
    try:
        import anndata
        import scanpy as sc
    except ImportError:
        print("AnnData or Scanpy not installed. Skipping AnnData creation.")
        return None
    
    np.random.seed(42)
    
    # Create expression matrix
    X = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(float)
    
    # Create obs (cell metadata)
    obs = pd.DataFrame({
        'cell_type': np.random.choice(['Type_A', 'Type_B', 'Type_C'], n_cells)
    }, index=[f'cell_{i}' for i in range(n_cells)])
    
    # Create var (gene metadata)
    var = pd.DataFrame({
        'gene_name': [f'gene_{i}' for i in range(n_genes)]
    }, index=[f'gene_{i}' for i in range(n_genes)])
    
    # Create AnnData object
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    
    return adata


def test_transfer_to_adata():
    """Test the transfer_to_adata method with mock data."""
    print("="*80)
    print("Testing cNMF transfer_to_adata method")
    print("="*80)
    
    # Import required modules
    try:
        from cnmf_extensions import cNMF
        import anndata
    except ImportError as e:
        print(f"Error: Required modules not available: {e}")
        print("Please install: pip install cnmf anndata scanpy")
        return False
    
    # Create temporary directory for mock results
    with tempfile.TemporaryDirectory() as tmpdir:
        print(f"\n1. Creating mock cNMF results in {tmpdir}...")
        
        # Parameters
        K = 5
        n_cells = 100
        n_genes = 50
        name = 'test_analysis'
        density_threshold = 0.1
        
        # Create mock cNMF results
        mock_results = create_mock_cnmf_results(
            tmpdir, name, K=K, n_cells=n_cells, n_genes=n_genes
        )
        print(f"   ✓ Created mock results with {n_cells} cells, {n_genes} genes, K={K}")
        
        # Create mock AnnData object
        print("\n2. Creating mock AnnData object...")
        adata = create_mock_adata(n_cells=n_cells, n_genes=n_genes)
        if adata is None:
            return False
        print(f"   ✓ Created AnnData with shape {adata.shape}")
        
        # Initialize cNMF object
        print("\n3. Initializing cNMF object...")
        cnmf_obj = cNMF(output_dir=tmpdir, name=name)
        print("   ✓ cNMF object initialized")
        
        # Test transfer_to_adata
        print("\n4. Testing transfer_to_adata method...")
        try:
            adata_with_results = cnmf_obj.transfer_to_adata(
                adata,
                K=K,
                density_threshold=density_threshold,
                n_top_genes=10
            )
            print("   ✓ Successfully transferred results to AnnData")
        except Exception as e:
            print(f"   ✗ Error during transfer: {e}")
            import traceback
            traceback.print_exc()
            return False
        
        # Verify results
        print("\n5. Verifying results...")
        
        # Check usage matrix
        if 'X_cnmf_usage' in adata_with_results.obsm:
            usage_shape = adata_with_results.obsm['X_cnmf_usage'].shape
            expected_shape = (n_cells, K)
            if usage_shape == expected_shape:
                print(f"   ✓ Usage matrix has correct shape: {usage_shape}")
            else:
                print(f"   ✗ Usage matrix shape mismatch: {usage_shape} != {expected_shape}")
                return False
        else:
            print("   ✗ Usage matrix not found in adata.obsm")
            return False
        
        # Check spectra
        if 'cnmf_spectra' in adata_with_results.varm:
            spectra_shape = adata_with_results.varm['cnmf_spectra'].shape
            expected_shape = (n_genes, K)
            if spectra_shape == expected_shape:
                print(f"   ✓ Spectra matrix has correct shape: {spectra_shape}")
            else:
                print(f"   ✗ Spectra shape mismatch: {spectra_shape} != {expected_shape}")
                return False
        else:
            print("   ✗ Spectra matrix not found in adata.varm")
            return False
        
        # Check top genes
        if 'cnmf_top_genes' in adata_with_results.uns:
            top_genes = adata_with_results.uns['cnmf_top_genes']
            print(f"   ✓ Top genes DataFrame has shape: {top_genes.shape}")
        else:
            print("   ✗ Top genes not found in adata.uns")
            return False
        
        # Check metadata
        if 'cnmf_top_genes_metadata' in adata_with_results.uns:
            metadata = adata_with_results.uns['cnmf_top_genes_metadata']
            print(f"   ✓ Metadata stored: K={metadata['K']}, "
                  f"density_threshold={metadata['density_threshold']}")
        else:
            print("   ✗ Metadata not found in adata.uns")
            return False
        
        print("\n" + "="*80)
        print("All tests passed! ✓")
        print("="*80)
        
        # Display summary
        print("\nSummary of what was added to AnnData:")
        print(f"  - adata.obsm['X_cnmf_usage']: {adata_with_results.obsm['X_cnmf_usage'].shape}")
        print(f"  - adata.varm['cnmf_spectra']: {adata_with_results.varm['cnmf_spectra'].shape}")
        print(f"  - adata.uns['cnmf_top_genes']: DataFrame with {top_genes.shape[0]} genes")
        print(f"  - adata.uns['cnmf_top_genes_metadata']: Analysis parameters")
        
        return True


if __name__ == '__main__':
    success = test_transfer_to_adata()
    if success:
        print("\n✓ Integration test completed successfully!")
    else:
        print("\n✗ Integration test failed!")
        exit(1)
