"""SDRF converter -- produces sample.parquet and run.parquet.

Reads an SDRF TSV file, loads it into DuckDB, and uses SQL transforms to
produce rows that conform to ``SampleSchema`` and ``RunSchema``.  The
results are written through ``SampleWriter`` and ``RunWriter``.
"""

from __future__ import annotations

import logging
import re
from functools import cached_property
from pathlib import Path

import pandas as pd

from qpx.converters.base import BaseConverter
from qpx.converters.channel_labels import is_ms1_label, normalize_label
from qpx.core.files import run_file_stem
from qpx.core.ontology import PublicOntology
from qpx.core.scores import field_ontology_entries
from qpx.core.sdrf import (
    get_name_from_complex_sdrf_value,
    parse_modification_parameters,
    parse_sdrf_kv,
    validate_sdrf_data_files,
)
from qpx.writers.run import RunWriter
from qpx.writers.sample import SampleWriter

logger = logging.getLogger(__name__)


def convert_sdrf_views(sdrf_path: str, output_folder: Path, prefix: str, ontology_entries: list[dict]) -> bool:
    """Write the sample and run views and add the run ontology terms to ``ontology_entries``.

    An SDRF the converter cannot handle (value, key, OS or type error) is logged
    and its partial files removed; returns whether both views were written.
    """
    try:
        with SdrfConverter() as sdrf_conv:
            sdrf_conv.convert(
                sdrf_path=sdrf_path,
                sample_output=str(output_folder / f"{prefix}.sample.parquet"),
                run_output=str(output_folder / f"{prefix}.run.parquet"),
            )
            ontology_entries.extend(sdrf_conv.run_ontology_entries())
    except (ValueError, KeyError, OSError, TypeError) as exc:
        logger.warning("SDRF conversion skipped (incomplete SDRF?): %s", exc)
        for suffix in (".sample.parquet", ".run.parquet"):
            corrupt = output_folder / f"{prefix}{suffix}"
            if corrupt.exists():
                corrupt.unlink()
                logger.debug("Removed corrupt %s", corrupt)
        return False
    return True


# Standard SDRF column names (lower-cased)
_COL_SOURCE_NAME = "source name"
_COL_DATA_FILE = "comment[data file]"
_COL_LABEL = "comment[label]"
_COL_FRACTION = "comment[fraction identifier]"
_COL_TECH_REP = "comment[technical replicate]"
_COL_BIO_REP = "characteristics[biological replicate]"
_COL_ORGANISM = "characteristics[organism]"
_COL_ORGANISM_PART = "characteristics[organism part]"
_COL_DISEASE = "characteristics[disease]"
_COL_CELL_LINE = "characteristics[cell line]"
_COL_CELL_TYPE = "characteristics[cell type]"
_COL_SEX = "characteristics[sex]"
_COL_AGE = "characteristics[age]"
_COL_DEV_STAGE = "characteristics[developmental stage]"
_COL_ANCESTRY = "characteristics[ancestry category]"
_COL_INDIVIDUAL = "characteristics[individual]"
_COL_INSTRUMENT = "comment[instrument]"
_COL_ENZYME = "comment[cleavage agent details]"
_COL_DISSOC = "comment[dissociation method]"
_COL_ACQUISITION = "comment[proteomics data acquisition method]"
_COL_GRADIENT = "comment[gradient length]"
_COL_CHANNEL_ROLE = "comment[channel role]"
_COL_REFERENCE_CHANNEL = "comment[reference channel]"
_COL_MIXTURE = "comment[mixture]"
_LABEL_MOD = re.compile(r"(?i)^(?:label:|dimethyl|mtraq)")
_DURATION = re.compile(r"(?i)(\d+(?:\.\d+)?)\s*(h|hr|hours?|m|min|mins|minutes?|s|sec|seconds?)?")
_SECONDS_PER_UNIT = {"h": 3600.0, "m": 60.0, "s": 1.0}
_COL_MOD_PARAM = "comment[modification parameter"  # prefix, may have ]s suffix
_KNOWN_RUN_COMMENT_COLUMNS = {
    _COL_DATA_FILE,
    _COL_LABEL,
    _COL_FRACTION,
    _COL_TECH_REP,
    _COL_INSTRUMENT,
    _COL_ENZYME,
    _COL_DISSOC,
    _COL_ACQUISITION,
    _COL_GRADIENT,
    _COL_CHANNEL_ROLE,
    _COL_REFERENCE_CHANNEL,
    _COL_MIXTURE,
}

# Sample fields the schema requires; a blank value becomes "not available" rather than losing the views
_REQUIRED_SAMPLE_FIELDS = {
    "organism": _COL_ORGANISM,
    "organism_part": _COL_ORGANISM_PART,
}

# Run fields whose SDRF values are ontology terms
_RUN_TERM_COLUMNS = {
    "instrument": _COL_INSTRUMENT,
    "enzymes": _COL_ENZYME,
    "dissociation_method": _COL_DISSOC,
    "acquisition_method": _COL_ACQUISITION,
}


def _extract_nt(value: str) -> str:
    """Extract the ``NT=`` name from a complex SDRF value."""
    if pd.isna(value):
        return ""
    return get_name_from_complex_sdrf_value(str(value))


def _extract_ac(value: str) -> str:
    """Extract the ``AC=`` accession from a complex SDRF value ('' when absent)."""
    return parse_sdrf_kv(value).get("AC", "")


def _collect_column_values(df: pd.DataFrame, substr: str) -> list[str]:
    """Collect unique values from all columns whose name contains *substr*."""
    cols = [c for c in df.columns if substr in c]
    if not cols:
        return []
    return pd.unique(df[cols].values.flatten()).tolist()


def _join_values(vals: list) -> str | None:
    """Join a list of values into a single '; '-separated string, or None."""
    clean = [str(v) for v in vals if pd.notna(v)]
    return "; ".join(clean) if clean else None


# Sentinel values that indicate an SDRF field has no meaningful content
_SENTINELS = {"not applicable", "not available", "none", "na", "n/a", ""}


def _is_sentinel(value: str | None) -> bool:
    """Return True if the value is an SDRF sentinel (not applicable, etc.)."""
    if value is None:
        return True
    return str(value).strip().lower() in _SENTINELS


def _all_sentinels(sdrf_df: pd.DataFrame, col_pattern: str) -> bool:
    """Return True if ALL unique values matching *col_pattern* are sentinels.

    Used to skip optional sample fields when every sample has "not applicable"
    or similar placeholder values.
    """
    vals = _collect_column_values(sdrf_df, col_pattern)
    if not vals:
        return True
    return all(_is_sentinel(v) for v in vals)


def _normalize_char_name(sdrf_col: str) -> str:
    """Convert an SDRF ``characteristics[X]`` or ``factor value[X]`` column name to snake_case.

    Example: ``characteristics[biological replicate]`` -> ``biological_replicate``
    """
    # Extract the name inside brackets
    m = re.match(r"(?:characteristics|factor value)\[(.+)\]", sdrf_col)
    name = m.group(1) if m else sdrf_col
    # Convert to snake_case
    name = re.sub(r"[^a-z0-9]+", "_", name.lower().strip())
    return name.strip("_")


def _merge_characteristic_values(group: pd.DataFrame, columns: list[str]) -> str | None:
    """Merge unique values, retaining sentinels only when no real value exists."""
    vals = pd.unique(group[columns].values.flatten()).tolist()
    meaningful = [v for v in vals if pd.notna(v) and not _is_sentinel(v)]
    return _join_values(meaningful or vals)


def _add_unique_run_label(
    seen_labels: dict[str, tuple[str, str]],
    run_file: str,
    label: str,
    sample_acc: str,
    label_raw: object,
) -> None:
    """Register a run label, rejecting duplicate canonical join keys."""
    if label in seen_labels:
        previous_sample, previous_label = seen_labels[label]
        raise ValueError(
            f"SDRF run {run_file!r} maps canonical label {label!r} more than once: "
            f"{previous_sample!r} ({previous_label!r}) and {sample_acc!r} ({label_raw!r})"
        )
    seen_labels[label] = (sample_acc, str(label_raw))


def _sample_characteristic_maps(sdrf_df: pd.DataFrame, skip_empty_properties: bool) -> tuple[dict, dict]:
    """Resolve standard and extra sample characteristic columns."""
    # Optional fields — only written when SDRF has meaningful values
    _optional_fields = {
        "disease": _COL_DISEASE,
        "cell_line": _COL_CELL_LINE,
        "cell_type": _COL_CELL_TYPE,
        "sex": _COL_SEX,
        "age": _COL_AGE,
        "developmental_stage": _COL_DEV_STAGE,
        "ancestry": _COL_ANCESTRY,
        "individual": _COL_INDIVIDUAL,
    }

    # Determine which optional fields have real data across the full SDRF
    active_optional: dict[str, str] = {}
    for field_name, col_pattern in _optional_fields.items():
        if skip_empty_properties and _all_sentinels(sdrf_df, col_pattern):
            continue
        active_optional[field_name] = col_pattern

    # Combine into the full char_map for per-sample extraction
    char_map = {**_REQUIRED_SAMPLE_FIELDS, **active_optional}

    # Discover extra characteristics columns not in the standard set
    # biological replicate belongs in run.parquet, not sample.parquet
    known_prefixes = (
        set(_REQUIRED_SAMPLE_FIELDS.values())
        | set(_optional_fields.values())
        | {
            "characteristics[sample description",
            _COL_BIO_REP,
        }
    )
    extra_char_cols = [
        c for c in sdrf_df.columns if c.startswith("characteristics[") and not any(c.startswith(k) for k in known_prefixes)
    ]
    # Group all SDRF columns that map to the same output characteristic.
    extra_col_map: dict[str, list[str]] = {}
    for col in extra_char_cols:
        extra_col_map.setdefault(_normalize_char_name(col), []).append(col)
    for col in sdrf_df.columns:
        if col.startswith("factor value["):
            extra_col_map.setdefault(f"factor_{_normalize_char_name(col)}", []).append(col)
    return char_map, extra_col_map


def _run_samples(group: pd.DataFrame, run_file: str) -> list[dict]:
    """Build run samples and their replicate values."""
    # Build per-run samples list
    samples = []
    seen_labels: dict[str, tuple[str, str]] = {}
    for row in group.to_dict("records"):
        label_raw = row.get(_COL_LABEL, "")
        # Normalize to the canonical channel label (e.g. the SDRF ontology
        # form "label free sample" -> "LFQ") so run.samples[].label is the
        # same join key that feature/pg intensities[].label carry.
        label = normalize_label(_extract_nt(label_raw)) if pd.notna(label_raw) else ""

        sample_acc = str(row.get(_COL_SOURCE_NAME, ""))
        _add_unique_run_label(seen_labels, run_file, label, sample_acc, label_raw)

        bio_rep = None
        if _COL_BIO_REP in row and pd.notna(row[_COL_BIO_REP]):
            try:
                bio_rep = int(row[_COL_BIO_REP])
            except (ValueError, TypeError):
                bio_rep = None

        tech_rep = None
        if _COL_TECH_REP in row and pd.notna(row[_COL_TECH_REP]):
            try:
                tech_rep = int(row[_COL_TECH_REP])
            except (ValueError, TypeError):
                tech_rep = None

        samples.append(
            {
                "sample_accession": sample_acc,
                "label": label,
                "biological_replicate": bio_rep,
                "technical_replicate": tech_rep,
                "channel_role": _channel_role(row, label),
                "mixture": str(row[_COL_MIXTURE]) if pd.notna(row.get(_COL_MIXTURE)) else None,
            }
        )

    return samples


def _gradient_seconds(group: pd.DataFrame) -> float | None:
    """``comment[gradient length]`` in seconds: ``120 min``, ``2 h``, ``7200 s``; a bare number is minutes."""
    if _COL_GRADIENT not in group.columns or group[_COL_GRADIENT].dropna().empty:
        return None
    match = _DURATION.fullmatch(str(group[_COL_GRADIENT].dropna().iloc[0]).strip())
    if match is None:
        return None
    return float(match.group(1)) * _SECONDS_PER_UNIT[(match.group(2) or "m")[0].lower()]


def _channel_role(row: dict, label: str) -> str | None:
    """A channel's declared role; a label matching ``comment[reference channel]`` is the reference."""
    role = row.get(_COL_CHANNEL_ROLE)
    if pd.notna(role) and role:
        return str(role)
    reference = row.get(_COL_REFERENCE_CHANNEL)
    if pd.isna(reference) or not reference:
        return None
    text = str(reference).strip()
    if text.lower() in ("yes", "true") or normalize_label(_extract_nt(text)) == label:
        return "reference"
    return "sample"


def _term_pairs(group: pd.DataFrame, column: str) -> list[tuple[str, str]]:
    """Return the ``(name, accession)`` pairs of a run's SDRF term column, in file order."""
    values = [str(value) for value in _collect_column_values(group, column) if pd.notna(value)]
    pairs = ((_extract_nt(value), _extract_ac(value)) for value in values)
    return [pair for pair in pairs if pair[0]]


class SdrfConverter(BaseConverter):
    """Convert an SDRF TSV into ``sample.parquet`` and ``run.parquet``.

    Usage::

        with SdrfConverter() as conv:
            conv.convert(
                sdrf_path="sdrf.tsv",
                sample_output="sample.parquet",
                run_output="run.parquet",
            )
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # Run fields that hold values, and the SDRF column each came from, for the ontology view
        self._run_field_sources: dict[str, str] = {}
        self._term_cache: dict[tuple[str, str], tuple[str, str] | None] = {}

    def convert(
        self,
        sdrf_path: str,
        sample_output: str | None = None,
        run_output: str | None = None,
        creator: str = "qpx",
    ) -> None:
        """Run the full SDRF conversion.

        Args:
            sdrf_path: Path to the SDRF TSV file.
            sample_output: Optional output path for ``sample.parquet``.
            run_output: Optional output path for ``run.parquet``.
            creator: Creator tag stored in Parquet metadata.
        """
        if sample_output is None and run_output is None:
            raise ValueError("At least one of sample_output or run_output is required")
        sdrf_df = self._read_sdrf(sdrf_path)
        if run_output is not None:
            validate_sdrf_data_files(sdrf_df)
        if sample_output is not None:
            self._write_samples(sdrf_df, sample_output, creator)
        if run_output is not None:
            self._write_runs(sdrf_df, run_output, creator)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _read_sdrf(self, sdrf_path: str) -> pd.DataFrame:
        """Load the SDRF file and normalise column names to lower-case."""
        df = pd.read_csv(sdrf_path, sep="\t", header=0)
        df.columns = df.columns.str.lower()
        return df

    # ------------------------------------------------------------------
    # Sample output
    # ------------------------------------------------------------------

    def _write_samples(
        self,
        sdrf_df: pd.DataFrame,
        output_path: str,
        creator: str,
        skip_empty_properties: bool = True,
    ) -> None:
        """Build one row per unique sample and write ``sample.parquet``.

        All ``characteristics[X]`` columns become individual string columns.
        Multi-valued fields are joined with ``"; "``.
        Repeated extra characteristics are merged; sentinel placeholders are
        ignored when meaningful values are present for that sample.

        Args:
            sdrf_df: Loaded SDRF DataFrame.
            output_path: Destination path for sample.parquet.
            creator: Creator tag stored in Parquet metadata.
            skip_empty_properties: When True (default), skip optional fields
                whose values are all sentinels ("not applicable", "not available",
                "none") across every sample. Set to False to write all fields.
        """

        char_map, extra_col_map = _sample_characteristic_maps(sdrf_df, skip_empty_properties)
        extra_column_names = sorted(extra_col_map)

        sample_groups = sdrf_df.groupby(_COL_SOURCE_NAME)

        records: list[dict] = []
        for sample_name, group in sample_groups:
            rec: dict = {"sample_accession": str(sample_name)}

            # Standard fields -> scalar strings
            for field_name, col_pattern in char_map.items():
                vals = _collect_column_values(group, col_pattern)
                rec[field_name] = _join_values(vals)
            blank = [field for field in _REQUIRED_SAMPLE_FIELDS if not rec[field]]
            for field_name in blank:
                self.logger.warning("Sample %s has no %s in the SDRF; writing 'not available'", sample_name, field_name)
                rec[field_name] = "not available"

            # Free-text description
            desc_col = "characteristics[sample description"
            if not (skip_empty_properties and _all_sentinels(sdrf_df, desc_col)):
                desc_vals = _collect_column_values(group, desc_col)
                rec["sample_description"] = _join_values(desc_vals)

            # Extra characteristics -> individual columns
            for norm_name, sdrf_cols in extra_col_map.items():
                rec[norm_name] = _merge_characteristic_values(group, sdrf_cols)

            records.append(rec)

        with SampleWriter(
            output_path,
            creator=creator,
            extra_columns=extra_column_names,
        ) as writer:
            writer.write_batch(records)

        self.logger.info(f"Wrote {len(records)} samples ({len(extra_column_names)} extra columns) to {output_path}")

    # ------------------------------------------------------------------
    # Run output
    # ------------------------------------------------------------------

    def _write_runs(
        self,
        sdrf_df: pd.DataFrame,
        output_path: str,
        creator: str,
    ) -> None:
        """Build one row per run and write ``run.parquet``."""

        # Strip file extension from data file column
        sdrf_df = sdrf_df.copy()
        sdrf_df["_run_file"] = sdrf_df[_COL_DATA_FILE].map(run_file_stem)
        sdrf_df["_file_name"] = sdrf_df[_COL_DATA_FILE].astype(str).str.strip()

        # Group by run file
        run_groups = sdrf_df.groupby("_run_file")

        extra_comment_cols: dict[str, list[str]] = {}
        for column in sdrf_df.columns:
            match = re.fullmatch(r"comment\[(.+)\]", column)
            if not match or column in _KNOWN_RUN_COMMENT_COLUMNS or column.startswith(_COL_MOD_PARAM):
                continue
            name = _normalize_char_name(f"characteristics[{match.group(1)}]")
            if name:
                extra_comment_cols.setdefault(f"x_sdrf_{name}", []).append(column)
        extra_column_names = sorted(extra_comment_cols)

        records: list[dict] = []
        for run_file, group in run_groups:
            run_file = str(run_file)

            samples = _run_samples(group, run_file)

            # Assay name / run accession
            run_accession = run_file

            acquisition_fields = self._run_acquisition_fields(group)

            # Original file name with extension (e.g. "S1_Frontal_1.raw")
            file_name_vals = group["_file_name"].dropna().unique()
            file_name = str(file_name_vals[0]) if len(file_name_vals) > 0 else None

            rec = {
                "run_accession": run_accession,
                "run_file_name": run_file,
                "file_name": file_name,
                "samples": samples,
                **acquisition_fields,
                "modification_parameters": parse_modification_parameters(group),
                "label_definitions": self._label_definitions(group),
                "gradient_length": _gradient_seconds(group),
            }
            for output_column, source_columns in extra_comment_cols.items():
                rec[output_column] = _merge_characteristic_values(group, source_columns)
            records.append(rec)

        with RunWriter(output_path, creator=creator, extra_columns=extra_column_names) as writer:
            writer.write_batch(records)

        self.logger.info(f"Wrote {len(records)} runs ({len(extra_column_names)} extra columns) to {output_path}")

    def _run_acquisition_fields(self, group: pd.DataFrame) -> dict:
        """Extract acquisition fields as plain strings, with their resolved terms in ``additional_terms``."""
        # Fraction
        fraction = None
        if _COL_FRACTION in group.columns:
            frac_vals = group[_COL_FRACTION].dropna().unique()
            if len(frac_vals) > 0:
                fraction = str(frac_vals[0])

        pairs = {field: _term_pairs(group, column) for field, column in _RUN_TERM_COLUMNS.items()}
        for field, values in pairs.items():
            if values:
                self._run_field_sources[field] = _RUN_TERM_COLUMNS[field]

        def first(field: str) -> str | None:
            return pairs[field][0][0] if pairs[field] else None

        return {
            "fraction": fraction,
            "instrument": first("instrument"),
            "enzymes": [name for name, _ in pairs["enzymes"]] or None,
            "dissociation_method": first("dissociation_method"),
            "acquisition_method": first("acquisition_method"),
            "additional_terms": self._additional_terms(pairs) or None,
        }

    # ------------------------------------------------------------------
    # Ontology enrichment — resolve run metadata to CV accessions
    # ------------------------------------------------------------------

    @cached_property
    def _psi_ms(self) -> PublicOntology:
        """The PSI-MS vocabulary, loaded once per conversion."""
        return PublicOntology("psi_ms", auto_update=False)

    def _resolve_term(self, name: str, accession: str) -> tuple[str, str] | None:
        """Resolve an SDRF term to ``(accession, name)``, trying the ``AC=`` accession before the ``NT=`` name.

        A PSI-MS accession must exist in the vocabulary. An accession of another
        ontology cannot be checked here, so it is kept with the SDRF's name.
        """
        key = (name, accession)
        if key not in self._term_cache:
            term = (self._psi_ms.lookup(accession) if accession else None) or self._psi_ms.lookup_name(name)
            if term is not None:
                self._term_cache[key] = (term["accession"], term["name"])
            else:
                self._term_cache[key] = (accession, name) if accession and not accession.startswith("MS:") else None
        return self._term_cache[key]

    def _additional_terms(self, pairs: dict[str, list[tuple[str, str]]]) -> list[dict]:
        """Build the run's ontology terms; ``properties`` names the plain field each one resolves."""
        terms = []
        for field, values in pairs.items():
            for name, accession in values:
                resolved = self._resolve_term(name, accession)
                if resolved is not None:
                    terms.append(
                        {
                            "accession": resolved[0],
                            "name": resolved[1],
                            "properties": [{"key": "field", "value": field}, {"key": "sdrf_value", "value": name}],
                        }
                    )
        return terms

    def run_ontology_entries(self) -> list[dict]:
        """Return ontology rows for the run fields this SDRF filled; their values resolve in ``additional_terms``."""
        return field_ontology_entries(view="run", resolved_mappings=self._run_field_sources, tool_name="SDRF")

    # ------------------------------------------------------------------
    # Modification parsing
    # ------------------------------------------------------------------

    @classmethod
    def _label_definitions(cls, group: pd.DataFrame) -> list[dict] | None:
        """The label modifications declared on the rows of each MS1 label of a run."""
        definitions = []
        for raw_label, rows in group.groupby(_COL_LABEL, sort=True):
            label = normalize_label(_extract_nt(raw_label))
            modifications = [mod for mod in parse_modification_parameters(rows) or [] if _LABEL_MOD.match(mod["name"] or "")]
            if is_ms1_label(label) and modifications:
                definitions.append({"label": label, "modifications": modifications})
        return definitions or None
