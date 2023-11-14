"""Test artificial data sets"""
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
    assert list(q.approximation.keys()) == [
        'features',
        'ncells',
        'avg',
        'neighborhood',
        'frac',
    ]
    assert sorted(list(q.approximation['neighborhood'].keys())) == [
        'avg',
        'convex_hull',
        'coords_centroid',
        'frac',
        'ncells',
    ]

