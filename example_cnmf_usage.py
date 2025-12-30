"""
Example script demonstrating how to use the cnmf_extensions module
to add the transfer_to_adata method to cNMF objects.

This script shows:
1. How to import the extended cNMF class
2. How to use transfer_to_adata with your cNMF analysis
3. How to access the results in your AnnData object
"""

import scanpy as sc
from cnmf_extensions import cNMF

def main():
    """
    Example workflow for using cNMF with the transfer_to_adata method.
    
    Note: This is a template. You'll need to:
    - Replace paths with your actual data paths
    - Adjust K and density_threshold to match your analysis
    - Complete the cNMF analysis steps (prepare, factorize, combine, consensus)
    """
    
    print("="*80)
    print("cNMF transfer_to_adata Example")
    print("="*80)
    
    # Step 1: Initialize cNMF object
    print("\n1. Initializing cNMF object...")
    cnmf_obj = cNMF(
        output_dir='./cnmf_output',  # Directory where cNMF saves results
        name='my_analysis'            # Name of your analysis
    )
    print("   ✓ cNMF object created")
    
    # Step 2: Run cNMF analysis (if not already done)
    print("\n2. Running cNMF analysis...")
    print("   Note: In a real workflow, you would run:")
    print("   - cnmf_obj.prepare(counts_fn=..., components=[5,6,7,8,9,10], ...)")
    print("   - cnmf_obj.factorize(worker_i=0, total_workers=1)")
    print("   - cnmf_obj.combine()")
    print("   - cnmf_obj.consensus(k=8, density_threshold=0.1)")
    print("   (Skipping for this example - assuming analysis is complete)")
    
    # Step 3: Load your AnnData object
    print("\n3. Loading AnnData object...")
    print("   In a real workflow:")
    print("   adata = sc.read_h5ad('./data/integrated/47.integrated.h5ad')")
    print("   (Skipping for this example)")
    
    # Step 4: Transfer cNMF results to AnnData
    print("\n4. Transferring cNMF results to AnnData...")
    print("   Usage:")
    print("   adata = cnmf_obj.transfer_to_adata(")
    print("       adata,")
    print("       K=8,                    # Number of gene expression programs")
    print("       density_threshold=0.1   # Density threshold from consensus step")
    print("   )")
    print("   ")
    print("   This will add:")
    print("   - adata.obsm['X_cnmf_usage']: Cell × Programs usage matrix")
    print("   - adata.varm['cnmf_spectra']: Gene × Programs spectra scores")
    print("   - adata.uns['cnmf_top_genes']: Top marker genes per program")
    print("   - adata.uns['cnmf_top_genes_metadata']: Analysis metadata")
    
    # Step 5: Access results
    print("\n5. Accessing results...")
    print("   # Get usage scores for all cells")
    print("   usage_matrix = adata.obsm['X_cnmf_usage']")
    print("   print(f'Usage matrix shape: {usage_matrix.shape}')")
    print("   ")
    print("   # Get gene spectra scores")
    print("   spectra = adata.varm['cnmf_spectra']")
    print("   print(f'Spectra shape: {spectra.shape}')")
    print("   ")
    print("   # Get top genes for each program")
    print("   top_genes = adata.uns['cnmf_top_genes']")
    print("   print('Top 5 genes for Program 1:')")
    print("   print(top_genes.iloc[:5, 0])")
    
    # Step 6: Visualization
    print("\n6. Visualization examples...")
    print("   # Compute UMAP")
    print("   sc.pp.neighbors(adata)")
    print("   sc.tl.umap(adata)")
    print("   ")
    print("   # Plot program usage on UMAP")
    print("   sc.pl.umap(adata, color=['cNMF_1', 'cNMF_2', 'cNMF_3'])")
    print("   ")
    print("   # Or plot directly from obsm")
    print("   import matplotlib.pyplot as plt")
    print("   fig, ax = plt.subplots()")
    print("   sc.pl.umap(adata, color=adata.obsm['X_cnmf_usage'][:, 0], ax=ax)")
    print("   plt.show()")
    
    print("\n" + "="*80)
    print("Example complete!")
    print("="*80)
    print("\nFor a working example, make sure you have:")
    print("1. Completed cNMF analysis with saved results")
    print("2. An AnnData object with the same cells used in cNMF")
    print("3. Matching gene names between cNMF and AnnData")
    print("\nSee CNMF_EXTENSIONS_README.md for more details.")


if __name__ == '__main__':
    # Quick verification that the module works
    print("Verifying cnmf_extensions module...")
    try:
        from cnmf_extensions import cNMF
        print("✓ Successfully imported cNMF from cnmf_extensions")
        
        if hasattr(cNMF, 'transfer_to_adata'):
            print("✓ transfer_to_adata method is available")
        else:
            print("✗ transfer_to_adata method NOT found")
            
        print("\nRunning example...\n")
        main()
        
    except ImportError as e:
        print(f"✗ Error importing cnmf_extensions: {e}")
        print("\nMake sure:")
        print("1. cnmf_extensions.py is in the same directory or in your Python path")
        print("2. cNMF library is installed: pip install cnmf")
