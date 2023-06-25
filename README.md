# scranPY ![image](https://github.com/sfortma2/scranPY/assets/56206488/a3146dfe-b299-4aae-9155-97f74c3a2da2=250x250)


A python implementation of r-scran::computeSumFactors: normalization by deconvolution for single-cell RNA-sequencing.
 
 
 

# Installation
```
pip install git+https://github.com/sfortma2/scranPY.git
```

# Basic Usage
```ruby
import scranPY
scranPY.compute_sum_factors(adata, clusters='clusters', parallelize=True, algorithm='CVXPY', max_size=3000, stopwatch=True, 
    plotting=True, lower_bound=0.2, normalize_counts=False, log1p=False, layer=None, verbose=False, save_plots_dir=None)
```


# Original Credit

L Lun, A.T., Bach, K. and Marioni, J.C., 2016. Pooling across cells to normalize single-cell RNA sequencing data with many zero counts. Genome biology, 17(1), pp.1-14.


