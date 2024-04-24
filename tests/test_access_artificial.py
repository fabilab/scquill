"""Test artificial data sets"""
import scquill


def test_300x500_init(shared_datadir):
    q = scquill.Approximation.read_h5(
        filename=shared_datadir / "artificial_300x500_approx.h5",
    )


def test_300x500_to_anndata_average(shared_datadir):
    q = scquill.Approximation.read_h5(
        filename=shared_datadir / "artificial_300x500_approx.h5",
    )
    adata = q.to_anndata(neighborhood=False)
    assert frozenset(adata.obs_names) == frozenset(
        [f"celltype_{i+1}" for i in range(10)]
    )


def test_300x500_to_anndata_neighborhood(shared_datadir):
    q = scquill.Approximation.read_h5(
        filename=shared_datadir / "artificial_300x500_approx.h5",
    )
    adata = q.to_anndata(neighborhood=True)
    assert frozenset(adata.obs_names) == frozenset(
        [str(i) for i in range(30)]
    )
    assert "convex_hulls" in adata.uns.keys()
    assert list(adata.obsm.keys()) == ["X_ncells", "X_umap"]


def test_300x500_to_anndata_coarse_grain_1level(shared_datadir):
    q = scquill.Approximation.read_h5(
        filename=shared_datadir / "artificial_300x500_disease_approx.h5",
    )
    adata = q.to_anndata(
        neighborhood=False,
        groupby=["celltype", "disease"],
    )
    assert list(adata.obs_names)[0] == 'celltype_1\tFalse'
    assert list(adata.obs_names)[1] == 'celltype_1\tTrue'

    adata_cg = scquill.coarse_grain_anndata(
        adata,
        "celltype",
    )
    assert frozenset(adata.obs_names) == frozenset(
        [f"celltype_{i+1}" for i in range(10)]
    )
