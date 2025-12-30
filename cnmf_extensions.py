"""
cNMF Extensions Module

This module provides extended functionality for the cNMF (consensus Non-negative Matrix Factorization) 
library, specifically adding methods to integrate cNMF results with AnnData objects commonly used 
in single-cell RNA-seq analysis.

Usage:
    from cnmf_extensions import cNMF
    
    # Use the extended cNMF class as normal
    cnmf_obj = cNMF(output_dir='./output', name='sample')
    # ... run cNMF analysis ...
    
    # Transfer results to AnnData object
    import scanpy as sc
    adata = sc.read_h5ad('data.h5ad')
    adata = cnmf_obj.transfer_to_adata(adata, K=8, density_threshold=0.1)
"""

import numpy as np
import pandas as pd
from cnmf import cNMF as _cNMF_base


class cNMF(_cNMF_base):
    """
    Extended cNMF class with additional methods for AnnData integration.
    
    This class inherits all methods from the original cNMF class and adds
    the transfer_to_adata method for seamless integration with scanpy/AnnData workflows.
    """
    
    def transfer_to_adata(self, adata, K=None, density_threshold=None, 
                          n_top_genes=100, norm_usage=True, 
                          usage_key='X_cnmf_usage', 
                          spectra_key='cnmf_spectra',
                          top_genes_key='cnmf_top_genes'):
        """
        Transfer cNMF results to an AnnData object.
        
        This method loads the cNMF results using the specified parameters and adds them to
        the AnnData object's .obsm (for usage matrix), .varm (for gene spectra), and .uns
        (for top genes and metadata) attributes.
        
        Parameters
        ----------
        adata : AnnData
            The AnnData object to which cNMF results will be added
            
        K : int, optional
            Number of programs (must be within the k values specified in previous steps).
            If None, will try to use the last K value that was computed.
            
        density_threshold : float, optional
            Threshold for filtering outlier spectra (must be within the values specified 
            in consensus step). If None, will use a default threshold.
            
        n_top_genes : int, optional (default=100)
            Number of top genes per program to return
            
        norm_usage : bool, optional (default=True)
            If True, normalize cNMF usages to sum to 1
            
        usage_key : str, optional (default='X_cnmf_usage')
            Key name for storing usage matrix in adata.obsm
            
        spectra_key : str, optional (default='cnmf_spectra')
            Key name for storing gene spectra scores in adata.varm
            
        top_genes_key : str, optional (default='cnmf_top_genes')
            Key name for storing top genes DataFrame in adata.uns
            
        Returns
        -------
        adata : AnnData
            The input AnnData object with cNMF results added
            
        Raises
        ------
        ValueError
            If K or density_threshold are not specified and cannot be inferred
        RuntimeError
            If cNMF results cannot be loaded
            
        Examples
        --------
        >>> from cnmf_extensions import cNMF
        >>> import scanpy as sc
        >>> 
        >>> # Initialize and run cNMF
        >>> cnmf_obj = cNMF(output_dir='./output', name='sample')
        >>> # ... run cNMF analysis steps ...
        >>> 
        >>> # Load AnnData and transfer results
        >>> adata = sc.read_h5ad('data.h5ad')
        >>> adata = cnmf_obj.transfer_to_adata(adata, K=8, density_threshold=0.1)
        >>> 
        >>> # Access the results
        >>> print(adata.obsm['X_cnmf_usage'])  # Cell x Programs usage matrix
        >>> print(adata.varm['cnmf_spectra'])  # Gene x Programs spectra scores
        >>> print(adata.uns['cnmf_top_genes'])  # Top genes per program
        """
        
        # Validate inputs
        if K is None:
            raise ValueError(
                "K (number of programs) must be specified. "
                "Please provide the K value used in your cNMF analysis."
            )
        
        if density_threshold is None:
            raise ValueError(
                "density_threshold must be specified. "
                "Please provide the density_threshold value used in your cNMF consensus step."
            )
        
        try:
            # Load cNMF results
            print(f"Loading cNMF results with K={K}, density_threshold={density_threshold}...")
            usage, spectra_scores, spectra_tpm, top_genes = self.load_results(
                K=K,
                density_threshold=density_threshold,
                n_top_genes=n_top_genes,
                norm_usage=norm_usage
            )
            
            # Verify dimensions match
            if usage.shape[0] != adata.n_obs:
                raise ValueError(
                    f"Number of cells in cNMF results ({usage.shape[0]}) does not match "
                    f"number of observations in AnnData ({adata.n_obs}). "
                    "Please ensure the cNMF analysis was run on the same dataset."
                )
            
            # Store usage matrix in .obsm
            print(f"Adding usage matrix to adata.obsm['{usage_key}'] with shape {usage.shape}")
            adata.obsm[usage_key] = usage.values if hasattr(usage, 'values') else usage
            
            # Store gene spectra in .varm (need to align genes)
            # Find common genes between cNMF results and AnnData
            cnmf_genes = spectra_scores.columns
            adata_genes = adata.var_names
            common_genes = adata_genes.intersection(cnmf_genes)
            
            if len(common_genes) == 0:
                raise ValueError(
                    "No common genes found between cNMF results and AnnData. "
                    "Please check that gene names match."
                )
            
            print(f"Found {len(common_genes)} common genes out of {len(adata_genes)} genes in AnnData")
            
            # Create aligned spectra matrix
            spectra_aligned = np.zeros((len(adata_genes), K))
            for i, gene in enumerate(adata_genes):
                if gene in common_genes:
                    spectra_aligned[i, :] = spectra_scores[gene].values
            
            print(f"Adding spectra scores to adata.varm['{spectra_key}'] with shape {spectra_aligned.shape}")
            adata.varm[spectra_key] = spectra_aligned
            
            # Store top genes and metadata in .uns
            adata.uns[top_genes_key] = top_genes
            adata.uns[f'{top_genes_key}_metadata'] = {
                'K': K,
                'density_threshold': density_threshold,
                'n_top_genes': n_top_genes,
                'norm_usage': norm_usage,
                'n_common_genes': len(common_genes)
            }
            
            # Add program names to the usage matrix for easier interpretation
            program_names = [f'cNMF_{i+1}' for i in range(K)]
            adata.uns[f'{usage_key}_names'] = program_names
            
            print(f"Successfully transferred cNMF results to AnnData object")
            print(f"  - Usage matrix: adata.obsm['{usage_key}'] ({usage.shape[0]} cells x {K} programs)")
            print(f"  - Spectra scores: adata.varm['{spectra_key}'] ({len(adata_genes)} genes x {K} programs)")
            print(f"  - Top genes: adata.uns['{top_genes_key}']")
            print(f"  - Metadata: adata.uns['{top_genes_key}_metadata']")
            
            return adata
            
        except FileNotFoundError as e:
            raise RuntimeError(
                f"Could not load cNMF results. Error: {e}\n"
                "Please ensure that cNMF analysis has been completed and results are saved "
                "in the expected output directory."
            )
        except Exception as e:
            raise RuntimeError(
                f"Error transferring cNMF results to AnnData: {e}\n"
                "Please check that your cNMF analysis completed successfully and "
                "the parameters match those used in the analysis."
            )


# For backwards compatibility and easier imports
__all__ = ['cNMF']
