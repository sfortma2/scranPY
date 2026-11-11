"""Frozen R reference for the QR path; not full-package equivalence."""
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scranPY import compute_sum_factors

def test_qr_pooling_against_r_reference():
    root=Path(__file__).parent/'fixtures'
    counts=np.loadtxt(root/'scran_counts.csv',delimiter=',').T
    expected=np.loadtxt(root/'scran_factors.csv',delimiter=',')
    a=ad.AnnData(counts.copy());a.obs['cluster']=pd.Categorical(np.repeat(['a','b'],120))
    compute_sum_factors(a,clusters='cluster',parallelize=False,algorithm='QR',min_mean=.1,plotting=False,stopwatch=False)
    np.testing.assert_allclose(a.obs['size_factors'],expected,atol=1e-8,rtol=1e-8)
    np.testing.assert_array_equal(a.X,counts)
