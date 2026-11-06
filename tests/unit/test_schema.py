"""Tests for QPX schema layer: FieldDef, ViewSchema, and all schema classes."""

import pyarrow as pa

from qpx.core.data import (
    DatasetSchema,
    FeatureSchema,
    FieldDef,
    MzSchema,
    OntologySchema,
    PgSchema,
    ProvenanceSchema,
    PsmSchema,
    RunSchema,
    SampleSchema,
    ViewSchema,
    XicSchema,
)
from tests.conftest import _placeholder

ALL_SCHEMAS = [
    FeatureSchema,
    PsmSchema,
    PgSchema,
    MzSchema,
    SampleSchema,
    RunSchema,
    DatasetSchema,
    OntologySchema,
    ProvenanceSchema,
    XicSchema,
]


def test_field_def():
    """FieldDef: nullable, optional, types, doc metadata, complex nested types."""
    # Required field not nullable
    f = FieldDef("my_col", pa.string(), nullable=False, doc="test")
    arrow = f.to_arrow_field()
    assert arrow.name == "my_col"
    assert arrow.type == pa.string()
    assert arrow.nullable is False

    # Nullable field
    f = FieldDef("score", pa.float64(), nullable=True, doc="test nullable")
    assert f.to_arrow_field().nullable is True

    # Optional field
    f = FieldDef("opt_col", pa.int32(), nullable=True, optional=True, doc="optional")
    assert f.to_arrow_field().nullable is True
    assert f.optional is True

    # Default field is nullable
    f = FieldDef("default_col", pa.string())
    assert f.to_arrow_field().nullable is True

    # Preserves arrow type
    f = FieldDef("scan", pa.list_(pa.int32()), nullable=False)
    assert f.to_arrow_field().type == pa.list_(pa.int32())

    # Doc in metadata
    f = FieldDef("documented_col", pa.string(), nullable=False, doc="My documentation")
    arrow = f.to_arrow_field()
    assert arrow.metadata is not None
    assert arrow.metadata[b"doc"] == b"My documentation"

    # No doc, no metadata
    f = FieldDef("no_doc_col", pa.string(), nullable=False)
    assert f.to_arrow_field().metadata is None

    # Complex nested type
    score_type = pa.list_(
        pa.struct(
            [
                pa.field("score_name", pa.string()),
                pa.field("score_value", pa.float64()),
                pa.field("higher_better", pa.bool_(), nullable=True),
            ]
        )
    )
    f = FieldDef("scores", score_type, nullable=True)
    assert f.to_arrow_field().type == score_type


def test_view_schema_basic():
    """ViewSchema: arrow schema, iter_fields, field names match, is_optional."""
    schema = FeatureSchema.get_arrow_schema()
    assert isinstance(schema, pa.Schema)

    fields = list(FeatureSchema._iter_fields())
    assert len(fields) > 0
    for name, field in fields:
        assert isinstance(name, str)
        assert isinstance(field, FieldDef)

    # Field names match
    field_names_from_schema = set(schema.names)
    field_names_from_iter = {name for name, _ in FeatureSchema._iter_fields()}
    assert field_names_from_schema == field_names_from_iter

    # is_optional
    assert FeatureSchema._is_optional("pg_global_qvalue") is True
    assert FeatureSchema._is_optional("gg_accessions") is True
    assert FeatureSchema._is_optional("sequence") is False
    assert FeatureSchema._is_optional("nonexistent_field") is False


def test_de_novo_fields_are_nullable():
    """De novo workflows have no target-decoy search or protein mapping."""
    assert PsmSchema.get_arrow_schema().field("is_decoy").nullable is False
    assert FeatureSchema.get_arrow_schema().field("is_decoy").nullable is False
    assert FeatureSchema.get_arrow_schema().field("anchor_protein").nullable is True


def test_schema_validation():
    """Validation: valid table, missing required, type mismatch, optional absent."""
    schema = FeatureSchema.get_arrow_schema()

    # Valid table returns no errors. Non-nullable columns must carry real
    # values (strict validation flags nulls in required columns as errors).
    arrays = {f.name: (pa.nulls(1, type=f.type) if f.nullable else pa.array([_placeholder(f.type)], type=f.type)) for f in schema}
    table = pa.table(arrays, schema=schema)
    assert FeatureSchema.validate(table) == []

    # Missing required column
    fields_to_keep = [f for f in schema if f.name != "sequence"]
    arrays = {f.name: pa.nulls(1, type=f.type) for f in fields_to_keep}
    table = pa.table(arrays, schema=pa.schema(fields_to_keep))
    errors = FeatureSchema.validate(table)
    assert len(errors) >= 1
    assert any("sequence" in e for e in errors)

    # Type mismatch
    wrong_fields = [pa.field(f.name, pa.string()) if f.name == "charge" else f for f in schema]
    arrays = {f.name: pa.nulls(1, type=f.type) for f in wrong_fields}
    arrays["charge"] = pa.array(["wrong_type"], type=pa.string())
    table = pa.table(arrays, schema=pa.schema(wrong_fields))
    errors = FeatureSchema.validate(table)
    assert any("charge" in e and "int16" in e for e in errors)

    # Optional column absent is fine
    fields_to_keep = [f for f in schema if f.name != "pg_global_qvalue"]
    arrays = {f.name: pa.nulls(1, type=f.type) for f in fields_to_keep}
    table = pa.table(arrays, schema=pa.schema(fields_to_keep))
    errors = FeatureSchema.validate(table)
    assert not any("pg_global_qvalue" in e for e in errors)

    # Required but nullable columns must be present even when their values are null.
    fields_to_keep = [f for f in schema if f.name != "anchor_protein"]
    arrays = {f.name: pa.nulls(1, type=f.type) for f in fields_to_keep}
    table = pa.table(arrays, schema=pa.schema(fields_to_keep))
    errors = FeatureSchema.validate(table)
    assert any("anchor_protein" in e for e in errors)


def test_run_with_only_required_columns_is_valid():
    """Run columns the spec does not require may be absent (bigbio/qpx#356)."""
    arrow = RunSchema.get_arrow_schema()
    keep = ["run_accession", "run_file_name", "samples"]
    sample = {"sample_accession": "s1", "label": "LFQ", "biological_replicate": None, "technical_replicate": None}
    table = pa.Table.from_pylist(
        [{"run_accession": "r1", "run_file_name": "r1", "samples": [sample]}],
        schema=pa.schema([arrow.field(name) for name in keep]),
    )
    for strict in (False, True):
        assert not RunSchema.validate_full(table, strict=strict).issues


def test_mz_and_psm_with_optional_columns_absent():
    """Mz and PSM columns may be absent in both validation modes (bigbio/qpx#356)."""
    cases = [
        (
            MzSchema,
            {
                "run_file_name",
                "scan",
                "inverse_ion_mobility",
                "drift_time",
                "compensation_voltage",
                "ion_injection_time",
                "precursors",
                "ion_mobility_array",
                "ion_mobility_array_unit",
                "cv_params",
            },
        ),
        (
            PsmSchema,
            {
                "mz_array",
                "intensity_array",
                "charge_array",
                "ion_type_array",
                "ion_mobility_array",
                "ion_mobility",
                "cross_links",
            },
        ),
    ]
    for schema_class, dropped in cases:
        fields = [f for f in schema_class.get_arrow_schema() if f.name not in dropped]
        arrays = {
            f.name: (pa.nulls(1, type=f.type) if f.nullable else pa.array([_placeholder(f.type)], type=f.type)) for f in fields
        }
        table = pa.table(arrays, schema=pa.schema(fields))
        for strict in (False, True):
            issues = schema_class.validate_full(table, strict=strict).issues
            assert not any(issue.check == "missing_column" for issue in issues)


def test_all_schemas_have_valid_primary_keys():
    """Every schema has a _primary_key tuple whose fields exist in the schema."""
    for schema_inst in ALL_SCHEMAS:
        assert isinstance(schema_inst._primary_key, tuple)
        assert len(schema_inst._primary_key) > 0
        schema = schema_inst.get_arrow_schema()
        for pk_field in schema_inst._primary_key:
            assert pk_field in schema.names, f"PK field '{pk_field}' not in {schema_inst.__name__}"


def test_yaml_loader_and_nested_types():
    """YAML loader returns ViewSchema; nested structs resolve; sample extension works."""
    from qpx.core.data.loader import load_schema

    schema = load_schema("feature")
    assert isinstance(schema, ViewSchema)

    # Nested struct: modifications = list<struct<name, accession, positions>>
    arrow_schema = FeatureSchema.get_arrow_schema()
    mod_type = arrow_schema.field("modifications").type
    assert isinstance(mod_type, pa.ListType)
    inner = mod_type.value_type
    assert isinstance(inner, pa.StructType)
    assert inner.get_field_index("name") >= 0
    assert inner.get_field_index("positions") >= 0

    # Nested struct: intensities = list<struct<label, intensity>>
    int_type = arrow_schema.field("intensities").type
    assert isinstance(int_type, pa.ListType)
    inner = int_type.value_type
    assert isinstance(inner, pa.StructType)
    assert inner.get_field_index("label") >= 0
    assert inner.get_field_index("intensity") >= 0

    # Sample extra columns
    extended = SampleSchema.get_extended_arrow_schema(["bmi", "biological_replicate"])
    assert "bmi" in extended.names
    assert extended.field("bmi").type == pa.string()
    assert extended.field("organism").type == pa.string()


def test_schema_mirrors_are_byte_identical():
    """The runtime schemas (qpx/core/data/schemas) and the published spec mirror
    (docs/spec/schemas) MUST stay byte-identical. They have silently drifted
    before (bigbio/qpx#220 review); this guards against re-drift so the docs
    always describe exactly what the library enforces."""
    from pathlib import Path

    repo = Path(__file__).resolve().parents[2]
    runtime = repo / "qpx" / "core" / "data" / "schemas"
    mirror = repo / "docs" / "spec" / "schemas"
    runtime_files = {p.name for p in runtime.glob("*.yaml")}
    mirror_files = {p.name for p in mirror.glob("*.yaml")}
    assert runtime_files == mirror_files, (
        f"schema file sets differ: only in runtime={runtime_files - mirror_files}, only in mirror={mirror_files - runtime_files}"
    )
    drifted = [name for name in sorted(runtime_files) if (runtime / name).read_bytes() != (mirror / name).read_bytes()]
    assert not drifted, (
        f"schema mirror drift in {drifted}: qpx/core/data/schemas and docs/spec/schemas must be byte-identical — sync them."
    )


def test_run_schema_acquisition_fields():
    """WS4: run schema exposes acquisition_method and additional_terms."""
    schema = RunSchema.get_arrow_schema()

    # acquisition_method: nullable string
    assert "acquisition_method" in schema.names
    acq = schema.field("acquisition_method")
    assert acq.type == pa.string()
    assert acq.nullable is True

    # additional_terms: nullable list<ontology_term struct>
    assert "additional_terms" in schema.names
    terms = schema.field("additional_terms")
    assert terms.nullable is True
    assert isinstance(terms.type, pa.ListType)
    inner = terms.type.value_type
    assert isinstance(inner, pa.StructType)
    assert inner.get_field_index("accession") >= 0
    assert inner.get_field_index("name") >= 0


def test_feature_schema_peptide_qvalue():
    """WS4: feature schema exposes a nullable peptide_qvalue float."""
    schema = FeatureSchema.get_arrow_schema()
    assert "peptide_qvalue" in schema.names
    qval = schema.field("peptide_qvalue")
    assert qval.type == pa.float64()
    assert qval.nullable is True


def test_missing_nullable_column_says_it_must_be_present():
    kept = [field for field in FeatureSchema.get_arrow_schema() if field.name != "modifications"]
    table = pa.table({field.name: pa.nulls(1, type=field.type) for field in kept}, schema=pa.schema(kept))
    messages = [issue.message for issue in FeatureSchema.validate_full(table).issues if issue.check == "missing_column"]
    assert messages == ["Missing column: modifications (its values may be null, but the column must be present)"]


def test_run_samples_without_the_later_nullable_members_stay_valid():
    """A 1.1 or 1.2 run view has no samples[].channel_role or mixture; a required member stays required."""
    written = pa.struct(
        [
            ("sample_accession", pa.string()),
            ("label", pa.string()),
            ("biological_replicate", pa.int32()),
            ("technical_replicate", pa.int32()),
        ]
    )
    expected = RunSchema.get_arrow_schema().field("samples").type

    def type_issues(samples_type):
        table = pa.table({"run_accession": ["r"], "run_file_name": ["r"], "samples": pa.nulls(1, type=pa.list_(samples_type))})
        return [issue for issue in RunSchema.validate_full(table).issues if issue.check == "type_mismatch"]

    assert expected.value_type.num_fields > written.num_fields
    assert type_issues(written) == []
    assert type_issues(pa.struct(list(written)[1:])) != []
    assert type_issues(pa.struct([*written, ("unknown", pa.string()), ("other", pa.string()), ("more", pa.string())])) != []
