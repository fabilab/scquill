<img src="https://raw.githubusercontent.com/fabilab/scquill/main/logo.png" width="150" height="150">

# scquill
Approximate any single cell data set, saving >99% of memory and runtime.


## Approximating a single cell data set
```python
import scquill

q = Quill(
    filename='myscdata.h5ad',
    output_filename='myapprox.h5',
    celltype_column="cell_annotation",
)

q()
```

**Steps:**
- Load dataset if necessary
- Preprocess
- Compress
- Store to output file
- (TODO): Provide an interface to explore approximations.

## Authors
Fabio Zanini @[fabilab](https://fabilab.org)
