import numpy as np
import matplotlib.pyplot as plt

from .dotplot import dotplot
from .neighborhood_composition import neighborhood_composition


def neighborhoodplot(
    accessor,
    features,
    groupby='celltype',
    measurement_type=None,
    axs=None,
    cmap='viridis_r',
    max_dot_size=100,
    **kwargs,
):
    """Composition of all neighborhoods."""

    if axs is None:
        adata = accessor.to_adata(
            groupby=groupby,
            neighborhood=True,
            measurement_type=measurement_type,
        )
        nnei = adata.n_obs
        nct = len(adata.uns['celltypes'])
        fig, axs = plt.subplots(
            1, 2,
            gridspec_kw=dict(
                width_ratios=[0.3 + len(features), nct],
                ),
            sharey=True,
            figsize=(0.5 + 0.3 * (len(features) + nct), 0.5 + 0.3 * nnei),
            )


    dotplot(
        accessor,
        features,
        groupby=groupby,
        measurement_type=measurement_type,
        neighborhood=True,
        ax=axs[0],
        max_dot_size=max_dot_size,
        cmap=cmap,
        **kwargs,
    )
    neighborhood_composition(
        accessor,
        groupby=groupby,
        measurement_type=measurement_type,
        ax=axs[1],
        cmap=cmap,
        **kwargs,
    )
    axs[0].set_ylabel(axs[1].get_ylabel())
    axs[1].set_ylabel('')
    fig.tight_layout()

    return axs
