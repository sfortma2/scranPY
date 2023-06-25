# scranPY <img src="https://github.com/sfortma2/scranPY/assets/56206488/1e6ae6f9-60df-48ec-8bbe-1f07ae6e1560" width="39.59" height="75">


A python implementation of r-scran::computeSumFactors: normalization by deconvolution for single-cell RNA-sequencing.
 
 
 

# Installation
```
pip install git+https://github.com/sfortma2/scranPY.git
```

# Basic Usage

```ruby
compute_size_factors(adata=AnnData, clusters=None, parallelize=True, algorithm='CVXPY', sizes=np.arange(21, 102, 5), 
   max_size=3000, min_mean=None, plotting=True, lower_bound=0.1, normalize_counts=False, log1p=False, layer='scranPY', 
   verbose=True, save_plots_dir=None, stopwatch=True):
  """
  Args:
    1. adata: An AnnData file (unnormalized, non-log transformed counts in active adata.X matrix).
    2. clusters: An observation in adata.obs containing cluster annotations (None or Str; default: None).
    3. parallelize: Whether to parallelize the computation (bool; default: True).
    4. algorithm: The algorithm to use for solving linear equations (str: 'CVXPY' or 'QR'; default: 'CVXPY').
        'QR' uses numpy.linalg.qr to compute qr factorization. This is analogous to the current implementation of r-scran::computeSumFactors but it is slow b/c requires dense matrices. 
        'CVXPY' uses the cvxpy module to construct a quadratic optimization function for solving linear equations. This is the faster, recommended option.
    5. sizes: A NumPy array of pool sizes to use for normalization (default: np.arange(21, 102, 5)).
    6. max_size: The maximum size of a cluster before being split into smaller chunks for computations (default: 3000).
    7. min_mean: The minimum mean gene expression level to consider for reference cells. If None, will automatically determine the appropriate min_mean (default: None).
    8. plotting: Whether to plot the results (bool; default: True).
    9. lower_bound: The lower bound for 'constraints' in the CVXPY algorithm. This is a hyperparameter that can increase the scaling of the smallest returned size factors. (range: 0 to 0.5; default: 0.1).
        i.e. constraints = [G @ x >= H, x >= lower_bound]
    10. normalize_counts: Whether to normalize the active adata.X matrix by dividing the matrix by the returned size factors. (bool; default: False).
    11. log1p: Whether to np.log1p transform the active adata.X matrix after normalization. Only relevant if normalize_counts=True. (bool; default: False).
    12. layer: The layer to store normalized counts (only relevant if normalize_counts=True). (str or None; default: 'scranPY').
    13. verbose: Whether to print progress while solving linear equations. (bool; default: True).
    14. save_plots_dir: The directory to save the plots. (str; default: None).
    15. stopwatch: Whether to time the function. (bool; default: True).

  Returns:
    A NumPy array of size factors and stores size factors in passed adata as adata.obs['size_factors'].
  """
```

```ruby
import scranPY
scranPY.compute_sum_factors(adata, clusters='clusters', parallelize=True, algorithm='CVXPY', max_size=3000, stopwatch=True, 
    plotting=True, lower_bound=0.2, normalize_counts=False, log1p=False, layer=None, verbose=False, save_plots_dir=None)
```



# Original Credit

L Lun, A.T., Bach, K. and Marioni, J.C., 2016. Pooling across cells to normalize single-cell RNA sequencing data with many zero counts. Genome biology, 17(1), pp.1-14.


