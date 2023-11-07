'''
Utility functions for the compression
'''
import gc
import numpy as np
import pandas as pd
import scanpy as sc


def approximate_dataset(
    adata,
    celltype_column,
    celltype_order,
    measurement_type="gene_expression",
):
    """Compress atlas for one tissue after data is clean, normalised, and reannotated."""
    # Celltype averages
    resd = _compress_celltype(
        adata,
        celltype_column,
        celltype_order,
        measurement_type,
    )

    # Local neighborhoods
    neid = _compress_neighborhoods(
        resd['ncells'],
        adata,
        celltype_column,
        measurement_type=measurement_type,
    )

    features = adata.var_names
    result = {
        'features': features,
        'celltype': {
            'ncells': resd['ncells'],
            'avg': resd['avg'],
            'neighborhood': neid,
        },
    }
    if measurement_type == "gene_expression":
        result['celltype']['frac'] = resd['frac']
        result['celltype']['neighborhood']['frac'] = neid['frac']

    return result


def _compress_celltype(
    adata,
    celltype_column,
    celltype_order,
    measurement_type,
):
    """Compress at the cell type level"""
    features = adata.var_names
    avg = pd.DataFrame(
            np.zeros((len(features), len(celltype_order)), np.float32),
            index=features,
            columns=celltype_order,
            )
    if measurement_type == "gene_expression":
        frac = pd.DataFrame(
                np.zeros((len(features), len(celltype_order)), np.float32),
                index=features,
                columns=celltype_order,
                )
    ncells = pd.Series(
            np.zeros(len(celltype_order), np.int64), index=celltype_order,
            )

    for celltype in celltype_order:
        idx = adata.obs[celltype_column] == celltype
        
        # Number of cells
        ncells[celltype] = idx.sum()

        # Average across cell type
        Xidx = adata[idx].X
        avg[celltype] = np.asarray(Xidx.mean(axis=0))[0]
        if measurement_type == "gene_expression":
            frac[celltype] = np.asarray((Xidx > 0).mean(axis=0))[0]

    res = {
        'avg': avg,
        'ncells': ncells,
    }
    if measurement_type == "gene_expression":
        res['frac'] = frac
    return res


def _compress_neighborhoods(
    ncells,
    adata,
    celltype_column,
    max_cells_per_type=300,
    measurement_type='gene_expression',
    avg_neighborhoods=3,
):
    """Compress local neighborhood of a single cell type."""
    # Try something easy first, like k-means
    from sklearn.cluster import KMeans
    from scipy.spatial import ConvexHull

    features = adata.var_names

    celltypes = list(ncells.keys())

    # Subsample with some regard for cell typing
    cell_ids = []
    for celltype, ncell in ncells.items():
        cell_ids_ct = adata.obs_names[adata.obs[celltype_column] == celltype]
        if ncell > max_cells_per_type:
            idx_rand = np.random.choice(range(ncell), size=max_cells_per_type, replace=False)
            cell_ids_ct = cell_ids_ct[idx_rand]
        cell_ids.extend(list(cell_ids_ct))
    adata = adata[cell_ids].copy()

    ##############################################
    # USE AN EXISTING EMBEDDING OR MAKE A NEW ONE
    emb_keys = ['umap', 'tsne']
    for emb_key in emb_keys:
        if f'X_{emb_key}' in adata.obsm:
            break
    else:
        emb_key = 'umap'

        # Log
        sc.pp.log1p(adata)

        # Select features
        sc.pp.highly_variable_genes(adata)
        adata.raw = adata
        adata = adata[:, adata.var.highly_variable]

        # Create embedding, a proxy for cell states broadly
        sc.tl.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        points = adata.obsm[f'X_{emb_key}']

        # Back to all features for storage
        adata = adata.raw.to_adata()
        adata.obsm[f'X_{emb_key}'] = points

        # Back to cptt or equivalent for storage
        adata.X.data = np.expm1(adata.X.data)
    ##############################################

    points = adata.obsm[f'X_{emb_key}']

    # Do a global clustering, ensuring at least 3 cells
    # for each cluster so you can make convex hulls
    for n_clusters in range(avg_neighborhoods * len(celltypes), 1, -1):
        kmeans = KMeans(
            n_clusters=n_clusters,
            random_state=0,
            n_init='auto',
        ).fit(points) 
        labels = kmeans.labels_

        # Book keep how many cells of each time are in each cluster
        tmp = adata.obs[[celltype_column]].copy()
        tmp['kmeans'] = labels
        tmp['c'] = 1.0
        ncells_per_label = (
                tmp.groupby(['kmeans', celltype_column])
                   .size()
                   .unstack(fill_value=0)
                   .loc[:, celltypes])
        del tmp

        if ncells_per_label.sum(axis=1).min() >= 3:
            break
    else:
        raise ValueError("Cannot cluster neighborhoods")

    n_neis = kmeans.n_clusters
    nei_avg = pd.DataFrame(
            np.zeros((len(features), n_neis), np.float32),
            index=features,
            )
    nei_coords = pd.DataFrame(
            np.zeros((2, n_neis), np.float32),
            index=['x', 'y'],
            )
    convex_hulls = []
    if measurement_type == "gene_expression":
        nei_frac = pd.DataFrame(
                np.zeros((len(features), n_neis), np.float32),
                index=features,
                )
    for i in range(kmeans.n_clusters):
        idx = kmeans.labels_ == i

        # Add the average expression
        nei_avg.iloc[:, i] = np.asarray(adata.X[idx].mean(axis=0))[0]
        # Add the fraction expressing
        if measurement_type == "gene_expression":
            nei_frac.iloc[:, i] = np.asarray((adata.X[idx] > 0).mean(axis=0))[0]

        # Add the coordinates of the center
        points_i = points[idx]
        nei_coords.iloc[:, i] = points_i.mean(axis=0)

        # Add the convex hull
        hull = ConvexHull(points_i)
        convex_hulls.append(points_i[hull.vertices])

    # Clean up
    del adata
    gc.collect()

    nei_avg.columns = ncells_per_label.index
    nei_coords.columns = ncells_per_label.index
    if measurement_type == "gene_expression":
        nei_frac.columns = ncells_per_label.index

    neid = {
        'ncells': ncells_per_label,
        'avg': nei_avg,
        'coords_centroid': nei_coords,
        'convex_hull': convex_hulls,
    }
    if measurement_type == "gene_expression":
        neid['frac'] = nei_frac

    return neid
