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
    additional_groupby_columns=tuple(),
    measurement_type="gene_expression",
):
    """Compress atlas for one tissue after data is clean, normalised, and reannotated."""
    # Celltype averages
    resd = _compress_average(
        adata,
        celltype_column,
        additional_groupby_columns,
        measurement_type,
    )

    # Local neighborhoods
    neid = _compress_neighborhoods(
        resd['ncells'],
        adata,
        celltype_column,
        additional_groupby_columns,
        measurement_type,
    )

    features = adata.var_names
    result = {
        'features': features,
        'ncells': resd['ncells'],
        'avg': resd['avg'],
        'obs': resd['obs'],
        'neighborhood': neid,
    }
    if measurement_type == "gene_expression":
        result['frac'] = resd['frac']
        result['neighborhood']['frac'] = neid['frac']

    compressed_atlas = {
        measurement_type: result,
    }

    return compressed_atlas


def _compress_average(
    adata,
    celltype_column,
    additional_groupby_columns,
    measurement_type,
):
    """Compress at the cell type level"""
    features = adata.var_names
    nfeatures = len(features)

    groupby_columns = [celltype_column] + list(additional_groupby_columns)
    tmp = adata.obs[groupby_columns].copy()
    tmp['c'] = 1
    ncells = tmp.groupby(groupby_columns).sum()['c']
    # NOTE: make sure it's a multi-index for consistency
    if len(groupby_columns) == 1:
        ncells.index = pd.MultiIndex.from_arrays(
            [ncells.index.values],
            names=(ncells.index.name,),
        )
    ngroups = len(ncells)

    # Metadata for groups
    obs = ncells.index.to_frame()
    obs.index = ncells.index

    avg = pd.DataFrame(
            np.zeros((nfeatures, ngroups), np.float32),
            index=features,
            columns=ncells.index,
            )
    if measurement_type == "gene_expression":
        frac = pd.DataFrame(
                np.zeros((nfeatures, ngroups), np.float32),
                index=features,
                columns=ncells.index,
                )

    for groupid in ncells.index:
        # Reconstruct indices of focal cells
        idx = adata.obs[celltype_column] == groupid[0]
        if len(groupby_columns) > 1:
            for col, groupid_i in zip(groupby_columns[1:], groupid[1:]):
                idx &= adata.obs[col] == groupid_i

        # Average across group
        Xidx = adata[idx].X
        avg[groupid] = np.asarray(Xidx.mean(axis=0))[0]
        if measurement_type == "gene_expression":
            frac[groupid] = np.asarray((Xidx > 0).mean(axis=0))[0]

    res = {
        'avg': avg,
        'obs': obs,
        'ncells': ncells,
    }
    if measurement_type == "gene_expression":
        res['frac'] = frac
    return res


def _compress_neighborhoods(
    ncells,
    adata,
    celltype_column,
    additional_groupby_columns,
    measurement_type='gene_expression',
    max_cells_per_type=300,
    avg_neighborhoods=3,
):
    """Compress local neighborhood of a single cell type."""
    # Try something easy first, like k-means
    from sklearn.cluster import KMeans
    from scipy.spatial import ConvexHull

    features = adata.var_names

    groupby_columns = list(ncells.index.names)

    # Subsample with some regard for cell typing
    cell_ids = []
    for groupid, ncell in ncells.items():
        # Reconstruct indices of focal cells
        idx = adata.obs[celltype_column] == groupid[0]
        if len(groupby_columns) > 1:
            for col, groupid_i in zip(groupby_columns[1:], groupid[1:]):
                idx &= adata.obs[col] == groupid_i

        cell_ids_ct = adata.obs_names[idx]
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
    for n_clusters in range(avg_neighborhoods * len(ncells), 1, -1):
        kmeans = KMeans(
            n_clusters=n_clusters,
            random_state=0,
            n_init='auto',
        ).fit(points) 
        labels = kmeans.labels_

        # Book keep how many cells of each time are in each cluster
        tmp = adata.obs[groupby_columns].copy()
        tmp['kmeans'] = labels
        tmp['c'] = 1.0
        ncells_per_label = (
                tmp.groupby(['kmeans'] + groupby_columns)
                   .size()
                   .unstack(0, fill_value=0)
                   .T
                   )
        del tmp

        # Ensure the order is the same as the averages
        if len(groupby_columns) == 1:
            ncells_per_label = ncells_per_label.loc[:, ncells.index.get_level_values(0)]
            ncells_per_label.columns = ncells.index
        else:
            ncells_per_label = ncells_per_label.loc[:, ncells.index]

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
