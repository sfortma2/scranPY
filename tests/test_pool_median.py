"""Even-gene pools must use the actual two middle order statistics."""
import numpy as np
from scranPY.main import forge_system


def test_pool_median_even_and_odd_gene_counts():
    rng=np.random.default_rng(42)
    for genes in (4,6,10,19,100):
        for _ in range(40):
            values=rng.uniform(.1,10,genes)
            _,_,result=forge_system(genes,1,values[None,:],[0,0],1,np.ones(genes))
            np.testing.assert_allclose(result[0],np.median(values),atol=1e-14,rtol=0)
