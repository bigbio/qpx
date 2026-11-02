"""SDRF converter -- produces sample.parquet and run.parquet.

Reads an SDRF TSV file, loads it into DuckDB, and uses SQL transforms to
produce rows that conform to ``SampleSchema`` and ``RunSchema``.  The
results are written through ``SampleWriter`` and ``RunWriter``.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import pandas as pd

from qpx.converters.base import BaseConverter
from qpx.converters.channel_labels import normalize_label
from qpx.core.files import run_file_stem
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
_COL_MOD_PARAM = "comment[modification parameter"  # prefix, may have ]s suffix


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
    """Convert an SDRF ``characteristics[X]`` column name to snake_case.

    Example: ``characteristics[biological replicate]`` -> ``biological_replicate``
    """
    # Extract the name inside brackets
    m = re.match(r"characteristics\[(.+)\]", sdrf_col)
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
    # Required fields — always written
    _required_fields = {
        "organism": _COL_ORGANISM,
        "organism_part": _COL_ORGANISM_PART,
    }

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
    char_map = {**_required_fields, **active_optional}

    # Discover extra characteristics columns not in the standard set
    # biological replicate belongs in run.parquet, not sample.parquet
    known_prefixes = (
        set(_required_fields.values())
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
            }
        )

    return samples


def _run_term(group: pd.DataFrame, column: str, terms: set[tuple[str, str]]) -> str | None:
    """Extract one run term and collect its ontology accession."""
    name = None
    values = _collect_column_values(group, column)
    if values:
        name = _extract_nt(str(values[0]))
        if name:
            terms.add((name, _extract_ac(str(values[0]))))
    return name


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
        # Track unique run metadata values for ontology enrichment
        self._run_instruments: set[tuple[str, str]] = set()
        self._run_enzymes: set[tuple[str, str]] = set()
        self._run_dissociation: set[tuple[str, str]] = set()

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
            compression=self._compression,
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

        # Parse modification parameters (same across all rows usually)
        mod_params = parse_modification_parameters(sdrf_df)

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
                "modification_parameters": mod_params or None,
            }
            records.append(rec)

        with RunWriter(output_path, creator=creator, compression=self._compression) as writer:
            writer.write_batch(records)

        self.logger.info(f"Wrote {len(records)} runs to {output_path}")

    def _run_acquisition_fields(self, group: pd.DataFrame) -> dict:
        """Extract acquisition fields and collect their ontology terms."""
        # Fraction
        fraction = None
        if _COL_FRACTION in group.columns:
            frac_vals = group[_COL_FRACTION].dropna().unique()
            if len(frac_vals) > 0:
                fraction = str(frac_vals[0])

        # Instrument (plain string; CV mapping in ontology.parquet)
        instrument = _run_term(group, _COL_INSTRUMENT, self._run_instruments)

        # Enzymes (list of plain strings)
        enzymes = None
        enz_vals = _collect_column_values(group, _COL_ENZYME)
        if enz_vals:
            enzyme_pairs = [(_extract_nt(str(v)), _extract_ac(str(v))) for v in enz_vals if pd.notna(v)]
            enzymes = [name for name, _ in enzyme_pairs]
            self._run_enzymes.update(pair for pair in enzyme_pairs if pair[0])

        # Dissociation method (plain string)
        dissoc = _run_term(group, _COL_DISSOC, self._run_dissociation)

        # Acquisition method (plain string; e.g. DDA / DIA)
        acquisition_method = None
        acq_vals = _collect_column_values(group, _COL_ACQUISITION)
        if acq_vals:
            acquisition_method = _extract_nt(str(acq_vals[0])) or None

        return {
            "fraction": fraction,
            "instrument": instrument,
            "enzymes": enzymes,
            "dissociation_method": dissoc,
            "acquisition_method": acquisition_method,
        }

    # ------------------------------------------------------------------
    # Ontology enrichment — resolve run metadata to CV accessions
    # ------------------------------------------------------------------

    def run_ontology_entries(self) -> list[dict]:
        """Generate ontology.parquet entries for run metadata terms.

        Resolves instrument, enzyme, and dissociation method terms to
        PSI-MS CV accessions using ``PublicOntology``. The SDRF ``AC=``
        accession takes precedence over the ``NT=`` name.
        """
        from qpx.core.ontology import PublicOntology

        ms = PublicOntology("psi_ms", auto_update=False)

        entries: list[dict] = []
        seen: set[str] = set()

        # A name with an accession sorts before the same name without one, so the accession is tried first.
        run_terms = self._run_instruments | self._run_enzymes | self._run_dissociation
        for name, accession in sorted(run_terms, key=lambda term: (term[0], not term[1], term[1])):
            if not name or name in seen:
                continue

            term = ms.lookup(accession) if accession else None
            if term is None:
                term = ms.lookup_name(name)
            if term is None:
                continue
            seen.add(name)

            entries.append(
                {
                    "field_name": name,
                    "ontology_name": term["name"],
                    "ontology_accession": term["accession"],
                    "ontology_source": term["source"],
                    "ontology_version": ms.version if ms.version != "unknown" else None,
                    "view": "run",
                    "description": term["definition"],
                }
            )

        return entries
