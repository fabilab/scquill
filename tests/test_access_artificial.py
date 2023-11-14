"""Test artificial data sets"""
import scquill


def test_300x500_init(shared_datadir):
    q = scquill.Approximation(
        filename=shared_datadir / 'artificial_300x500_approx.h5',
    )

def test_300x500_to_adata_average(shared_datadir):
    q = scquill.Approximation(
        filename=shared_datadir / 'artificial_300x500_approx.h5',
    )
    adata = q.to_adata(neighborhood=False)

    print(adata.obs_names)

    assert frozenset(adata.obs_names) == frozenset([f'celltype_{i+1}' for i in range(10)])

def test_300x500_to_adata_neighborhood(shared_datadir):
    q = scquill.Approximation(
        filename=shared_datadir / 'artificial_300x500_approx.h5',
    )
    adata = q.to_adata(neighborhood=True)
    assert frozenset(adata.obs_names) == frozenset([f'neighborhood_{i+1}' for i in range(30)])
    assert 'convex_hulls' in adata.uns.keys()
    assert list(adata.obsm.keys()) == ['X_ncells', 'X_umap']
