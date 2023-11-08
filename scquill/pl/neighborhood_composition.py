import numpy as np
import matplotlib.pyplot as plt


def neighborhood_composition(
    accessor,
    groupby='celltype',
    measurement_type=None,
    ax=None,
    cmap='viridis_r',
    **kwargs,
):
    """Composition of all neighborhoods."""

    if ax is None:
        ax = plt.gca()

    adata = accessor.to_adata(
        groupby=groupby,
        neighborhood=True,
        measurement_type=measurement_type,
    )

    Xncells = adata.obsm['X_ncells']
    celltypes = adata.uns['celltypes']
    
    nnei, nct = Xncells.shape

    ax.set_xlim(-0.5, nct - 0.5)
    ax.set_ylim(-0.5, nnei - 0.5)
    xs, ys = np.meshgrid(np.arange(nct), np.arange(nnei))
    hs = 1.0 * Xncells / Xncells.max()

    facecolors = plt.cm.get_cmap(cmap)(hs.ravel())
    ax.bar(
        x=xs.ravel(),
        height=hs.ravel(),
        bottom=(ys - hs / 2).ravel(),
        color=facecolors,
        edgecolor=tuple([0.1] * 3 + [0.5]),
    )
    ax.set_xticks(np.arange(nct))
    ax.set_xticklabels(celltypes, rotation=90)
    ax.set_yticks(np.arange(nnei))
    ax.set_ylabel('Neighborhood')
    ax.figure.tight_layout()

    return ax
