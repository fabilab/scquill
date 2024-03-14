import os
import pathlib
from datetime import date
import argparse
import multiprocess as mp
import numpy as np
import pandas as pd
import cellxgene_census
import scquill


def func(dataset_id):
    import scquill
    import cellxgene_census

    output_fn = f'data/approximations/{dataset_id}.h5'
    if pathlib.Path(output_fn).exists():
        print(f'Already completed, skipping: {dataset_id}', flush=True)
        return dataset_id

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
    print(f'Access completed: {dataset_id}')

    print(f'Dataset size: {adata.shape}', flush=True)
    
    print(f'Build approximation for dataset: {dataset_id}')
    q = scquill.Compressor(
        adata=adata,
        celltype_column='cell_type',
        output_filename=output_fn,
    )
    q()
    return dataset_id


def pool_callback(dataset_id):
    print(f'Finished pool job for dataset: {dataset_id}', flush=True)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--threads', type=int, default=1,
        help='Number of processes to use in parallel',
    )
    parser.add_argument(
        '--test', action='store_true',
        help='Test only the first dataset, on 1 core',
    )
    parser.add_argument(
        '--dry', action='store_true',
        help='Dry run, don\'t call anything except the dataset_id generation',
    )
    args = parser.parse_args()

    today = date.today().strftime("%Y-%m-%d")
    datasets_fn = pathlib.Path(f'data/cellxgene_census_dataset_id_{today}.csv')

    if not pathlib.Path('data').exists():
        os.mkdir('data')
    if not pathlib.Path('data/approximations').exists():
        os.mkdir('data/approximations')
    
    if not datasets_fn.exists():
        print('Access CENSUS database to retrieve dataset ids')
        with cellxgene_census.open_soma() as census:
            dataset_ids = census["census_data"]["homo_sapiens"].obs.read(
                column_names = ["dataset_id"]
            ).concat()
            # They are already sorted by size, now invert it
            dataset_id_count = dataset_ids['dataset_id'].to_pandas().value_counts()[::-1]
        print('Access completed')
        
        today = date.today().strftime("%Y-%m-%d")
        dataset_id_count.to_csv(datasets_fn, sep=',')
        print('Done')
    
    print('Load dataset ids')
    dataset_id_count = pd.read_csv(datasets_fn, index_col=0)['count']
    dataset_ids = dataset_id_count.index.values

    if args.dry:
        print(f"Dry run with {args.threads} processes and test={args.test}")
    elif args.test:
        for dataset_id in dataset_ids:
            func(dataset_id)
            break
    elif args.threads == 1:
        ndata = len(dataset_ids)
        for i, dataset_id in enumerate(dataset_ids):
            ncells = dataset_id_count.loc[dataset_id]
            print(f'Start dataset {i+1} / {ndata}: {ncells} cells', flush=True)
            func(dataset_id)
    else:
        with mp.Pool(args.threads) as pool:
            # Get promise for good coding style
            async_res = pool.map_async(func, dataset_ids, callback=pool_callback)
            # Wait for the end on this process
            async_res.wait()
