import numpy as np
import matplotlib.pyplot as plt


def dotplot(
    accessor,
    features,
    groupby='celltype',
    neighborhood=False,
    measurement_type=None,
    ax=None,
    cmap='viridis_r',
    max_dot_size=100,
    **kwargs,
):
    """Make a dot plot of features x groups."""
    if len(features) == 0:
        return

    if ax is None:
        ax = plt.gca()

    adata = accessor.to_adata(
        groupby=groupby,
        neighborhood=neighborhood,
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

    # Create x and y coordinates in a mesh
    ngroups, nfeatures = hue.shape
    xs, ys = np.meshgrid(np.arange(nfeatures), np.arange(ngroups))

    ax.scatter(
        xs.ravel(),
        ys.ravel(),
        s=size.ravel(),
        c=hue.ravel(),
        cmap=cmap,
    )
    ax.set_xlim(-0.5, nfeatures - 0.5)
    ax.set_xticks(np.arange(nfeatures))
    ax.set_yticks(np.arange(ngroups))
    ax.set_xticklabels(features, rotation=90)
    ax.set_yticklabels(adata.obs_names)
    ax.figure.tight_layout()

    return ax
