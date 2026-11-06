# YAML Schema Reference

QPX data views have formal schema definitions in YAML format. These schemas define the structure and types of each field using Arrow-native type names and serve as the canonical, language-independent schema definition for the QPX format.

!!! note "Expression views use AnnData"
    The expression views (Absolute Expression, Differential Expression) use [AnnData](https://anndata.readthedocs.io/) (`.h5ad`) as their primary format and do not have YAML schemas. See [AnnData Concepts](anndata.md) for details.

!!! note "API views have no schemas"
    [API views](views.md) (e.g., peptide summaries, protein summaries) are computed on demand from the primary data views. They are programmable and do not have formal schemas.

!!! info "YAML vs Parquet"
    YAML schemas (`.yaml`) define the **logical data model** -- the field names, Arrow types, nullability, and documentation for each view. The actual data is serialized in **Apache Parquet** format for efficient columnar storage and querying. The YAML schemas are the single source of truth for the QPX data model, while Parquet handles the physical storage layer.

## How to use these schemas

- **Validation**: Use these schemas to validate that a QPX file conforms to the expected structure.
- **Code generation**: The QPX loader converts YAML schemas to PyArrow schemas at import time.
- **Documentation**: Each field includes a `doc` attribute describing its purpose and semantics.
- **Shared types**: Common struct types (e.g., `score`, `modification`, `intensity`) are defined once in [`types.yaml`](schemas/types.yaml) and reused across schemas.

## Schemas

### Data Views

=== "PSM"

    ```yaml title="psm.yaml"
    --8<-- "docs/spec/schemas/psm.yaml"
    ```

=== "Feature"

    ```yaml title="feature.yaml"
    --8<-- "docs/spec/schemas/feature.yaml"
    ```

=== "Protein Group"

    ```yaml title="pg.yaml"
    --8<-- "docs/spec/schemas/pg.yaml"
    ```

=== "Mass Spectra"

    ```yaml title="mz.yaml"
    --8<-- "docs/spec/schemas/mz.yaml"
    ```

=== "Extracted-Ion Chromatograms"

    ```yaml title="xic.yaml"
    --8<-- "docs/spec/schemas/xic.yaml"
    ```

### Metadata Views

=== "Dataset"

    ```yaml title="dataset.yaml"
    --8<-- "docs/spec/schemas/dataset.yaml"
    ```

=== "Sample"

    ```yaml title="sample.yaml"
    --8<-- "docs/spec/schemas/sample.yaml"
    ```

=== "Run"

    ```yaml title="run.yaml"
    --8<-- "docs/spec/schemas/run.yaml"
    ```

=== "Ontology"

    ```yaml title="ontology.yaml"
    --8<-- "docs/spec/schemas/ontology.yaml"
    ```

=== "Provenance"

    ```yaml title="provenance.yaml"
    --8<-- "docs/spec/schemas/provenance.yaml"
    ```

### Shared Types

```yaml title="types.yaml"
--8<-- "docs/spec/schemas/types.yaml"
```

!!! tip "Viewing schemas programmatically"
    You can load and inspect any schema in Python using the QPX loader:

    ```python
    from qpx.core.data.loader import load_schema

    schema = load_schema("feature")
    for field in schema.get_arrow_schema():
        print(f"{field.name}: {field.type} (nullable={field.nullable})")
    ```
