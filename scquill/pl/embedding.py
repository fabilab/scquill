import numpy as np
import matplotlib.pyplot as plt


def _axis_heuristic(nplots):
    nrows = nplots // 4 + int(bool(nplots % 4))
    ncols = 1 + nplots % 4
    return (nrows, ncols)


def embedding(
    accessor,
    features,
    groupby='celltype',
    measurement_type=None,
    key='umap',
    axs=None,
    cmap='viridis_r',
    max_dot_size=100,
    **kwargs,
):
    """Show a neighborhood embedding."""
    if len(features) == 0:
        return

    if axs is None:
        nrows, ncols = _axis_heuristic(len(features))
        fig, axs = plt.subplots(nrows, ncols)
        if max(nrows, ncols) == 1:
            axs = [axs]
        elif min(nrows, ncols) > 1:
            axs = axs.ravel()

    adata = accessor.to_adata(
        groupby=groupby,
        neighborhood=True,
        measurement_type=measurement_type,
    )

    hue = adata[:, features].X.copy()
    # TODO: detect log?
    if 'fraction' in adata.layers:
        size = adata[:, features].layers['fraction'].copy()
    else:
        size = hue.copy()

    # Fractions are already between 0 and 1. Scale them remembering that
    # ax.scatter has the funky square thing
    size = max_dot_size * size**2

    centers = adata.obsm[f'X_{key}']
    
    # Make an embedding for each feature
    for i, ax in enumerate(axs):
        # Scatter the centers
        huei = hue[:, i]
        sizei = size[:, i]
        ax.scatter(
            *centers.T,
            c=huei,
            s=sizei,
            cmap=cmap,
        )

        # Draw the text and convex hulls
        hulls = adata.uns['convex_hulls']
        huemin = huei.min()
        huemax = huei.max() * 1.01 + 1e-10
        for j, hull in enumerate(hulls):
            # Hull
            hueij = huei[j]
            hueij = (hueij - huemin) / (huemax - huemin)
            norm = hueij
            hueij = plt.cm.get_cmap(cmap)(hueij)
            facecolor = tuple(list(hueij)[:3] + [0.4])
            edgecolor = tuple(list(hueij)[:3] + [0.9])
            poly = plt.Polygon(
                hull,
                facecolor=facecolor,
                edgecolor=edgecolor,
            )
            ax.add_patch(poly)

            # Text
            ax.text(
                *centers[j],
                str(j),
                ha='center',
                va='center',
                color='black' if norm < 0.6 else 'white',
            )

        feature = features[i]
        ax.set_title(feature)
    fig.tight_layout()

    return axs
