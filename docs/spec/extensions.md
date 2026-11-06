# Extensions

Three places hold what a view has no field for. Pick by who else should read the value.

| Place | Holds | Readers |
|-------|-------|---------|
| `cv_params` | A string per row, named below | Anyone can list it; the type is the reader's concern |
| `x_<namespace>_<name>` column | A string column of your tool, typed by the reader | Queryable like any column, and the validator accepts it |
| `uns` of an AnnData view, or a Parquet footer key | Results about the whole file | Your own tooling |

A column named otherwise makes `qpxc validate` warn (`undeclared_column`), and a writer refuses or drops it with a warning.

## `cv_params` names

An entry is `{cv_name, cv_value}`. Use the snake_case name of a source column or a PSI-MS term, and put a tool's private meaning behind its name: `diann_quantity_quality`. [Scores & CV Terms](scores.md) lists the shared names, and [Converter Coverage](converter-coverage.md) those each converter writes. A value used for filtering or joining belongs in a column, where the type and the statistics help.

## Extension columns

`x_` columns are declared to a writer, and are nullable strings:

```python
with FeatureWriter(path, extra_columns=["x_mytool_cluster"]) as writer:
    writer.write_batch(records)
```

The `sample` view also takes any column: the SDRF's `characteristics[...]` and `factor value[...]` arrive there as `organism_part` and `factor_<name>`.

## Targeted runs: PRM, SRM and heavy standards

QPX has no Skyline converter and no field for a spiked amount, and a heavy channel needs a sample, or the `pg` check fails. This convention validates:

- Give the standard a channel and its own row in the `sample` view, `sample_accession` `SIS-1`, with `x_` or `factor_` columns for the amount and unit.
- On the feature row, the heavy peptidoform's `cv_params` holds `spiked_amount`, `spiked_amount_unit` and `paired_peptidoform`, the endogenous peptide it stands for.
- A transition's area is a row of the [fragment view](fragment.md).

## Curve fits: thermal profiling and dose response

The raw values per temperature or dose fit the views: one channel or run per point, with the point in `factor_<name>` of the `sample` view. A fit's result (Tm, EC50) describes a protein across the design, so it goes in the `.de.h5ad` or `.ae.h5ad` view:

```python
adata.uns["fit_results"] = pd.DataFrame(
    {
        "protein": ["P12345"],
        "condition": ["treated"],
        "model": ["sigmoid"],
        "parameter": ["Tm"],
        "value": [52.3],
        "stderr": [0.4],
        "r2": [0.98],
        "n_points": [10],
    }
)
```

This is a recommendation, which `qpxc validate` does not read.
