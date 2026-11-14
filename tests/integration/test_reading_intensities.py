"""Ways to read nested intensities give the same values."""

import pandas as pd
import pytest

import qpx


def test_arrow_backed_frame_has_arrow_dtypes(dataset_dir):
    with qpx.open_dataset(str(dataset_dir)) as ds:
        frame = ds.feature.to_df(dtype_backend="pyarrow")
        assert isinstance(frame["intensities"].dtype, pd.ArrowDtype)
        assert len(frame) == ds.feature.count()
        with pytest.raises(ValueError, match="dtype_backend"):
            ds.feature.to_df(dtype_backend="numpy")


def test_intensities_wide_matches_the_long_table(dataset_dir):
    with qpx.open_dataset(str(dataset_dir)) as ds:
        long = ds.feature.peptide_intensities().to_df()
        wide = ds.feature.intensities_wide().to_df()
        labels = sorted(long["label"].unique())
        assert set(labels) <= set(wide.columns)
        assert wide[labels].sum().sum() == pytest.approx(long["intensity"].sum())
        only = ds.feature.intensities_wide(labels=labels[:1]).to_df()
        assert labels[0] in only.columns
        assert set(labels[1:]).isdisjoint(only.columns)
