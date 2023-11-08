[![PyPI version](https://badge.fury.io/py/scquill.svg)](https://badge.fury.io/py/scquill)

<img src="https://raw.githubusercontent.com/fabilab/scquill/main/logo.png" width="150" height="150">

# scquill
Approximate any single cell data set, saving >99% of memory and runtime.


## Approximating a single cell data set
```python
import scquill

q = scquill.Compressor(
    filename='myscdata.h5ad',
    output_filename='myapprox.h5',
    celltype_column="cell_annotation",
)

q()
```

## Exploring an approximation
To load an approximation:
```
import scquill

app = scquill.ApproxAnnData(
    filename='myapprox.h5',
)
```

To show a neighborhood plot:
```
scquill.pl.neighborhoodplot(app, ['gene1', 'gene2', 'gene3'])
```
<img src="https://raw.githubusercontent.com/fabilab/scquill/main/neighborhoodplot.png" width="350">

To show embeddings of cell neighborhoods, similar to single-cell UMAPs:

```
scquill.pl.embedding(app, ['gene1', 'gene2', 'gene3'])
```
<img src="https://raw.githubusercontent.com/fabilab/scquill/main/embeddings.png" width="650">

**MORE TO COME**

## Authors
Fabio Zanini @[fabilab](https://fabilab.org)
