"""Test artificial data sets"""
import numpy as np
import scquill


def test_300x500_init(shared_datadir):
    q = scquill.Compressor(
        filename=shared_datadir / 'artificial_300x500.h5ad',
        output_filename='/tmp/artificial_300x500_approx.h5',
        celltype_column='celltype',
        include_neighborhood=True,
    )

def test_300x500_run(shared_datadir):
    q = scquill.Compressor(
        filename=shared_datadir / 'artificial_300x500.h5ad',
        output_filename='/tmp/artificial_300x500_approx.h5',
        celltype_column='celltype',
        include_neighborhood=True,
    )
    q()

    assert list(q.approximation.keys()) == ['gene_expression']
    approximation = q.approximation['gene_expression']
    assert list(approximation.keys()) == [
        'var_names',
        'obs',
        'obs_names',
        'Xave',
        'Xfrac',
        'neighborhood',
    ]
    assert list(approximation['neighborhood'].keys()) == [
        'obs_names',
        'cell_count',
        'coords_centroid',
        'convex_hull',
        'Xave',
        'Xfrac',
    ]
    assert approximation['obs'].columns.tolist() == ['celltype', 'cell_count']


def test_300x500_adjust(shared_datadir):
    q = scquill.Compressor(
        filename=shared_datadir / 'artificial_300x500.h5ad',
        output_filename='/tmp/artificial_300x500_disease_approx.h5',
        celltype_column='celltype',
        include_neighborhood=True,
    )
    q.prepare()
    q.adata.obs['disease'] = np.arange(2).astype(bool)[np.random.randint(2, size=q.adata.n_obs)]
    q.additional_groupby_columns = ['disease']
    q.compress()

    assert list(q.approximation.keys()) == ['gene_expression']
    approximation = q.approximation['gene_expression']
    assert list(approximation.keys()) == [
        'var_names',
        'obs',
        'obs_names',
        'Xave',
        'Xfrac',
        'neighborhood',
    ]
    assert list(approximation['neighborhood'].keys()) == [
        'obs_names',
        'cell_count',
        'coords_centroid',
        'convex_hull',
        'Xave',
        'Xfrac',
    ]
    assert approximation['obs'].columns.tolist() == ['celltype', 'disease', 'cell_count']
