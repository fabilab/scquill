from typing import Sequence, Union
import h5py
import hdf5plugin
import numpy as np
import pandas as pd


def get_feature_sequences(
    app,
    measurement_type : Union[None, str] = None,
    features : Union[None, Sequence[str]] = None,
):
    """Get the sequences of some features if available."""

    if app.filename is not None:

        # FIXME: outsource this to a separate function in scquill.io.read
        if measurement_type is None:
            with h5py.File(app.filename) as h5_data:
                measurement_type = app._infer_measurement_type(h5_data)

        with h5py.File(app.filename) as h5_data:
            if measurement_type not in h5_data["measurements"]:
                raise KeyError("Measurement type not found: {measurement_type}")

            me = h5_data["measurements"][measurement_type]

            # Get the feature names
            var_names = me["var_names"].asstr()[:]

            if "feature_sequences" not in me:
                raise KeyError("No feature sequences found.")

            # Get the pointer to the sequences
            seq_pointer = me["feature_sequences"]["sequences"].asstr()

            # Construct a mapping from feature names to indices
            var_series = pd.Series(np.arange(len(var_names)), index=var_names)

            if features is None:
                seqs = seq_pointer[:]
            else:
                # HDF5 requires sorted numerical indices
                fea_idx = var_series.loc[features].copy().to_frame(name="idx")
                fea_idx = fea_idx.sort_values("idx")
                fea_idx["idx_sorted"] = np.arange(len(fea_idx))

                # Extract the sequences
                seqs = seq_pointer[fea_idx["idx"].values]

                # Reorder the sequences
                seqs = seqs[fea_idx.loc[features, "idx_sorted"].values]

    else:
        raise NotImplementedError
                
    return seqs
