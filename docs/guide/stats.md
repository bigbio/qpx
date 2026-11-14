# Statistics and Quality Control

Compute statistics and quality metrics of a QPX dataset with the Python API. Every example here runs on a converted dataset in `./output/`.

## Getting Started

```python
import qpx

ds = qpx.Dataset("./output/")

print(ds.available_structures)
print(f"PSMs: {ds.psm.count()}")
print(f"Features: {ds.feature.count()}")
print(f"Protein group rows: {ds.pg.count()}")
```

A view opens lazily: `count()`, `filter()` and `select()` run in DuckDB without loading the file, and `to_df()` loads only what the query asks for.

## Built-in Summaries

```python
import qpx

ds = qpx.Dataset("./output/")

# One row: unique peptides and proteins, runs, features (targets only)
print(ds.qc_view.metrics().to_df())

# One row per run: n_proteins, n_peptides, n_features
print(ds.identification_summary.summary().to_df())

# One row per modification: n_psms, n_peptides
print(ds.modification_view.frequency().to_df())

# Through the sample table
print(ds.sample_summary.peptides_per_sample().to_df())
print(ds.sample_summary.proteins_per_sample().to_df())
```

## Your Own Queries

`ds.sql()` runs DuckDB SQL over the views, which are registered by name:

```python
import qpx

ds = qpx.Dataset("./output/")

# Charge states of target features
charges = ds.sql("""
    SELECT charge, count(*) AS n
    FROM feature
    WHERE NOT coalesce(is_decoy, false)
    GROUP BY charge
    ORDER BY charge
""").to_df()
print(charges)

# Peptides and peptidoforms
counts = ds.sql("""
    SELECT count(DISTINCT sequence) AS peptides, count(DISTINCT peptidoform) AS peptidoforms
    FROM feature
""").fetchone()
print(counts)
```

!!! warning "`unique` is a reserved word"
    Write `"unique"` in SQL when you select that feature column.

The same results are available through the query builder:

```python
import qpx

ds = qpx.Dataset("./output/")

print(ds.feature.filter("charge = 2").count())
print(ds.feature.filter("charge = 2").select("sequence", "charge").limit(5).to_df())
print(ds.psm.targets_only().select("peptidoform").to_df()["peptidoform"].nunique())
```

## Missing Values

`design_matrix()` pivots protein quantities to one row per sample and one column per protein. A missing quantity is `NaN`.

```python
import qpx

ds = qpx.Dataset("./output/")

matrix = ds.design_matrix()
print(f"{matrix.shape[0]} samples, {matrix.shape[1]} proteins")
print(f"Missing overall: {matrix.isna().to_numpy().mean():.1%}")
print(matrix.isna().mean(axis=1).sort_values(ascending=False).head())
```

## A Quality Report

```python
import qpx


def quality_report(path, min_peptides_per_run=1000):
    """Flag runs with few peptides, against the dataset's own median."""
    ds = qpx.Dataset(path)
    runs = ds.identification_summary.summary().to_df()
    median = runs["n_peptides"].median()
    low = runs[(runs["n_peptides"] < min_peptides_per_run) | (runs["n_peptides"] < median / 2)]
    return {
        "metrics": ds.qc_view.metrics().to_df().iloc[0].to_dict(),
        "median_peptides_per_run": float(median),
        "low_runs": low["run_file_name"].tolist(),
    }


print(quality_report("./output/"))
```

The thresholds are examples: choose them for your instrument and gradient. Do not take a fixed PSM count as a sign of quality across experiments.

## Validating Alongside Statistics

`validate()` returns one result per view. Use it before reporting numbers:

```python
import qpx

ds = qpx.Dataset("./output/")

for structure, result in ds.validate().items():
    errors = [issue for issue in result.issues if issue.severity == "error"]
    print(structure, "valid" if not errors else f"{len(errors)} error(s)")
```

## Comparing Datasets

Ids and identities are per analysis ([Row Identity](../spec/identity.md#scope)), so compare natural measures. To compare several conversions, loop over their folders:

```python
import pandas as pd
import qpx

rows = []
for path in ["./output/"]:
    ds = qpx.Dataset(path)
    rows.append({"dataset": path, **ds.qc_view.metrics().to_df().iloc[0].to_dict()})
print(pd.DataFrame(rows))
```

## Related Documentation

- [Visualization Guide](visualize.md)
- [Validate](validate.md)
- [API Views](../spec/views.md)
