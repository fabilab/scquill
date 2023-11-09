from datetime import date
import numpy as np
import pandas as pd
import cellxgene_census
import scquill

if False:
    print('Access CENSUS database to retrieve dataset ids')
    with cellxgene_census.open_soma() as census:
        dataset_ids = census["census_data"]["homo_sapiens"].obs.read(
            column_names = ["dataset_id"]
        ).concat()
        dataset_ids = dataset_ids['dataset_id'].unique().to_numpy(zero_copy_only=False)
    print('Access completed')
    
    today = date.today().strftime("%Y-%m-%d")
    with open(f'data/cellxgene_census_dataset_id_{today}.csv', 'wt') as f:
        f.write('\n'.join(dataset_ids))
    print('Done')

print('Load dataset ids')
today = date.today().strftime("%Y-%m-%d")
with open(f'data/cellxgene_census_dataset_id_{today}.csv', 'rt') as f:
    dataset_ids = np.array(
        f.read().split('\n'),
    )

for dataset_id in dataset_ids:
    print(f'Access CENSUS database to retrieve data for dataset: {dataset_id}')
    with cellxgene_census.open_soma() as census:
        adata = cellxgene_census.get_anndata(
            census = census,
            organism = "Homo sapiens",
            obs_value_filter = f"dataset_id == '{dataset_id}'",
            column_names = {
                "obs": ["assay", "cell_type", "tissue", "tissue_general", "suspension_type", "disease"],
            },
        )
    print('Access completed')
    
    print(f'Build approximation for dataset: {dataset_id}')
    q = scquill.Compressor(
        adata=adata,
        celltype_column='cell_type',
        output_filename=f'data/approximations/{dataset_id}.h5',
    )
    q()
    #break

