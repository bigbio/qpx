# Project Metadata

QPX keeps project metadata in tables next to the data views. Examples run on a converted dataset in `./output/`.

| Table | Holds | Source |
|---|---|---|
| `dataset` | Project accession, title, software, file checksums | Converter options, PRIDE |
| `sample` | Organism, tissue, disease and other properties, per sample | SDRF |
| `run` | Instrument, fraction, enzyme, the samples of each run | SDRF |
| `provenance` | Processing steps and parameters | Converter |
| `ontology` | CV terms of the columns | Converter |

## Dataset

```python
import qpx

ds = qpx.Dataset("./output/")

meta = ds.dataset_meta.to_df().iloc[0]
print(meta["project_accession"], meta["project_title"])
print(meta["software_name"], meta["software_version"])
print(meta["analysis_id"])
```

The specification version is in each file's footer, not in a column ([Versioning](../spec/versioning.md#version-in-files)).

## Samples

```python
import qpx

ds = qpx.Dataset("./output/")

samples = ds.sample.to_df()
print(samples[["sample_accession", "organism", "organism_part"]].head())
print(samples["organism"].value_counts())
```

A property with several values is one string joined with `"; "` ([SDRF Mapping](../spec/sdrf-mapping.md#multi-valued-metadata)). `ds.sample.by_organism("Homo sapiens")` selects samples by organism.

## Runs

A run lists its samples, one entry per channel:

```python
import qpx

ds = qpx.Dataset("./output/")

runs = ds.run.to_df()
print(runs[["run_file_name", "instrument", "fraction"]].head())

channels = ds.sql("""
    SELECT r.run_file_name, s.sample_accession, s.label
    FROM run r, UNNEST(r.samples) AS t(s)
""").to_df()
print(channels.head())
```

## Provenance

```python
import qpx

ds = qpx.Dataset("./output/")

for _, step in ds.provenance.to_df().iterrows():
    print(step["step_order"], step["step_name"], step["tool_name"], step["tool_version"])
```

## SDRF

Converters take the SDRF with `--sdrf-file`, write `sample` and `run` from it and keep a copy as `{prefix}.sdrf.tsv`. To read it directly, use the [sdrf-pipelines](https://github.com/bigbio/sdrf-pipelines) reader or pandas:

```python
from pathlib import Path

import pandas as pd

sdrf_path = next(Path("./output/").glob("*.sdrf.tsv"))
sdrf = pd.read_csv(sdrf_path, sep="\t")
print(sdrf.columns.tolist())
```

## PRIDE Metadata

`--enrich-pride` on a convert command, or `Dataset.enrich_from_pride()`, fills the dataset table from PRIDE:

```bash
qpxc convert maxquant \
    --evidence-file evidence.txt \
    --sdrf-file metadata.sdrf.tsv \
    --project-accession PXD014414 \
    --enrich-pride \
    --output-folder ./output
```

## Exporting Metadata

```python
from pathlib import Path

import qpx

ds = qpx.Dataset("./output/")
Path("./metadata").mkdir(exist_ok=True)

ds.sample.to_df().to_csv("./metadata/samples.tsv", sep="\t", index=False)
ds.dataset_meta.to_df().drop(columns=["file_checksums"], errors="ignore").to_csv(
    "./metadata/dataset.tsv", sep="\t", index=False
)
```

## Related Documentation

- [Dataset Specification](../spec/dataset.md)
- [Sample Specification](../spec/sample.md)
- [Run Specification](../spec/run.md)
- [Provenance Specification](../spec/provenance.md)
- [SDRF Mapping](../spec/sdrf-mapping.md)
