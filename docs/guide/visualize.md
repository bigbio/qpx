# Visualization Guide

The summary views have a `plot()` method that returns a [plotly](https://plotly.com/python/) `Figure`. Install the extra with `pip install "qpx[plotting]"`. Examples run on a converted dataset in `./output/`.

## Built-in Plots

```python
from pathlib import Path

import qpx

ds = qpx.Dataset("./output/")
Path("./plots").mkdir(exist_ok=True)

ds.identification_summary.plot().write_html("./plots/identifications.html")
ds.run_summary.plot().write_html("./plots/runs.html")
ds.modification_view.plot(top_n=10).write_html("./plots/modifications.html")
ds.qc_view.plot().write_html("./plots/qc.html")
ds.sample_summary.plot("peptides_per_sample").write_html("./plots/peptides_per_sample.html")
```

| View | Plot | `plot()` shows |
| --- | --- | --- |
| `ds.identification_summary` | Grouped bar | Proteins and peptides per run |
| `ds.run_summary` | Grouped bar | Peptides, proteins and features per run |
| `ds.modification_view` | Bar | The `top_n` modifications by PSM count |
| `ds.qc_view` | Bar | Unique peptides and proteins, runs, features |
| `ds.sample_summary` | Bar | `metric`: `peptides_per_run`, `peptides_per_sample` or `proteins_per_sample` |

A plotly figure is an interactive HTML page. For SVG, PDF or PNG, install [kaleido](https://github.com/plotly/Kaleido) and call `fig.write_image("plot.svg")`. A plotly figure has no `savefig`.

## Your Own Plots

Each view returns plain data, so any plotting library works. With matplotlib:

### Intensity Distributions

```python
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import qpx

ds = qpx.Dataset("./output/")
Path("./plots").mkdir(exist_ok=True)

matrix = ds.design_matrix()
values = [np.log2(row.dropna()) for _, row in matrix.iterrows()]

fig, ax = plt.subplots(figsize=(10, 5))
ax.boxplot(values, tick_labels=matrix.index.tolist())
ax.set_ylabel("log2 protein intensity")
ax.tick_params(axis="x", rotation=45)
fig.savefig("./plots/intensity_boxplot.svg", bbox_inches="tight")
```

### Peptides per Protein

```python
from pathlib import Path

import matplotlib.pyplot as plt
import qpx

ds = qpx.Dataset("./output/")
Path("./plots").mkdir(exist_ok=True)

counts = ds.sql("""
    SELECT anchor_protein, count(DISTINCT sequence) AS peptides
    FROM feature
    WHERE anchor_protein IS NOT NULL AND NOT coalesce(is_decoy, false)
    GROUP BY anchor_protein
""").to_df()

fig, ax = plt.subplots()
ax.hist(counts["peptides"], bins=range(1, 40))
ax.set_xlabel("Peptides per protein")
ax.set_ylabel("Proteins")
fig.savefig("./plots/peptides_per_protein.svg", bbox_inches="tight")
```

### Missing Values

```python
from pathlib import Path

import matplotlib.pyplot as plt
import qpx

ds = qpx.Dataset("./output/")
Path("./plots").mkdir(exist_ok=True)

missing = ds.design_matrix().isna().mean(axis=1) * 100

fig, ax = plt.subplots()
ax.bar(missing.index.astype(str), missing.to_numpy())
ax.set_ylabel("Missing proteins (%)")
ax.tick_params(axis="x", rotation=45)
fig.savefig("./plots/missing_values.svg", bbox_inches="tight")
```

### Volcano Plots

A differential expression file (`.de.h5ad`) holds one table per contrast in `uns["de_results"]`, as the [Differential Expression spec](../spec/differential.md) describes. Read it with `anndata` and plot `logfoldchanges` against `-log10(pvals_adj)`.

## Related Documentation

- [Statistics and Quality Control](stats.md)
- [API Views](../spec/views.md)
