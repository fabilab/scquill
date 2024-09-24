"""Compress single cell data using embedding neighborhoods."""

import gc
import numpy as np
import pandas as pd
import scanpy as sc


def _compress_neighborhood(
    adata,
    obs_bio,
    celltype_column,
    additional_groupby_columns,
    measurement_type="gene_expression",
    max_cells_per_type=300,
    avg_neighborhoods=3,
):
    """Compress local neighborhood of a single cell type."""
    # Try something easy first, like k-means
    from sklearn.cluster import KMeans
    from scipy.spatial import ConvexHull

    features = adata.var_names

    groupby_columns = [celltype_column] + list(additional_groupby_columns)

    # Subsample with some regard for cell typing
    # NOTE: this can be debated, but is also designed to get rid of "outlier cells"
    # which are the opposite of an approximation in a sense
    cell_ids = []
    idx = np.zeros(adata.n_obs, bool)
    for _, row in obs_bio.iterrows():
        # Reconstruct indices of focal cells
        idx[:] = True
        for col in groupby_columns:
            idx &= adata.obs[col] == row[col]

        cell_ids_ct = adata.obs_names[idx]
        ncell = row["cell_count"]
        if ncell > max_cells_per_type:
            idx_rand = np.random.choice(
                range(ncell), size=max_cells_per_type, replace=False
            )
            cell_ids_ct = cell_ids_ct[idx_rand]
        cell_ids.extend(list(cell_ids_ct))
    adata = adata[cell_ids].copy()

    ##############################################
    # USE AN EXISTING EMBEDDING OR MAKE A NEW ONE
    emb_keys = ["umap", "tsne"]
    for emb_key in emb_keys:
        if f"X_{emb_key}" in adata.obsm:
            break
    else:
        emb_key = "umap"

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
        points = adata.obsm[f"X_{emb_key}"]

        # Back to all features for storage
        adata = adata.raw.to_adata()
        adata.obsm[f"X_{emb_key}"] = points

        # Back to cptt or equivalent for storage
        adata.X.data = np.expm1(adata.X.data)
    ##############################################

    points = adata.obsm[f"X_{emb_key}"]

    # Do a global clustering, ensuring at least 3 cells
    # for each cluster so you can make convex hulls
    if len(groupby_columns) > 1:
        row_order = pd.MultiIndex.from_frame(obs_bio[groupby_columns])
    else:
        row_order = pd.Index(obs_bio[groupby_columns[0]])
    for n_clusters in range(avg_neighborhoods * len(obs_bio), 1, -1):
        kmeans = KMeans(
            n_clusters=n_clusters,
            random_state=0,
            n_init="auto",
        ).fit(points)
        labels = kmeans.labels_

        # Book keep how many cells of each time are in each cluster
        tmp = adata.obs[groupby_columns].copy()
        tmp["kmeans"] = labels
        tmp["c"] = 1.0
        ncells_per_label = (
            tmp.groupby(["kmeans"] + groupby_columns, observed=False)
            .size()
            .unstack(0, fill_value=0)
            .T
        )
        del tmp

        # Ensure the order is the same as the averages
        ncells_per_label = ncells_per_label.loc[:, row_order]

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
        index=["x", "y"],
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
        nei_avg.iloc[:, i] = np.asarray(adata.X[idx].mean(axis=0))[0].astype(np.float32)
        # Add the fraction expressing
        if measurement_type == "gene_expression":
            nei_frac.iloc[:, i] = np.asarray((adata.X[idx] > 0).mean(axis=0))[0].astype(
                np.float32
            )

        # Add the coordinates of the center
        points_i = points[idx]
        nei_coords.iloc[:, i] = points_i.mean(axis=0)

        # Add the convex hull
        hull = ConvexHull(points_i)
        convex_hulls.append(points_i[hull.vertices])

    # TODO: Approximate vector field for velocity in this embedding if available.
    # If not, we cannot recount the matrix outselves but if they layers are present,
    # we could recompute it and splash it onto our embedding, THEN recycle the arrow
    # flow algo from mpl to approximate

    # Clean up
    del adata
    gc.collect()

    nei_avg.columns = ncells_per_label.index
    nei_coords.columns = ncells_per_label.index
    if measurement_type == "gene_expression":
        nei_frac.columns = ncells_per_label.index

    neid = {
        "kind": "neighborhood",
        "obs_names": nei_avg.index.values,
        "cell_count": ncells_per_label.values,
        "coords_centroid": nei_coords.values,
        "convex_hull": convex_hulls,
        "Xave": nei_avg.values,
    }
    if measurement_type == "gene_expression":
        neid["Xfrac"] = nei_frac.values

    return neid
