"""Convert MSNet parquet data into a unified QPX dataset.

Reads MSNet projects, transforms PSM data to QPX schema, extracts
Sample/Run metadata from SDRF, and writes a partitioned QPX dataset
(by organism + run_file_name).

Usage:
    python scripts/convert_msnet.py --input /path/to/msnet --output /path/to/output --projects PXD014877-Sulfolobus_solfataricus
    python scripts/convert_msnet.py --input /path/to/msnet --output /path/to/output --all
"""

from __future__ import annotations

import argparse
import csv
import datetime
import logging
import re
import threading
import uuid
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pyarrow as pa
import pyarrow.dataset as pds
import pyarrow.parquet as pq

from qpx._version import __version__
from qpx.core.scores import field_ontology_entries
from qpx.core.sql import escape_path, sql_build
from qpx.writers.base import parquet_write_options
from qpx.writers.ontology import OntologyWriter
from qpx.writers.partitioned import (
    _PARTITIONED_MAX_ROWS_PER_GROUP,
    _PARTITIONED_MIN_ROWS_PER_GROUP,
    _stamp_footer_metadata,
)
from qpx.writers.provenance import ProvenanceWriter

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s [%(threadName)s] %(message)s",
    datefmt="%H:%M:%S",
)
_log = logging.getLogger(__name__)

DEFAULT_MSNET_DIR = "msnet"
DEFAULT_OUTPUT_DIR = "msnet_qpx"


# ---------------------------------------------------------------------------
# SDRF parsing
# ---------------------------------------------------------------------------


def _strip_raw_ext(filename: str) -> str:
    """Remove .raw / .mzML / .d etc. extensions from a filename."""
    return re.sub(r"\.(raw|mzML|d|wiff|RAW)$", "", filename, flags=re.IGNORECASE)


def parse_sdrf(sdrf_path: Path) -> list[dict]:
    """Parse an SDRF TSV and return a list of row dicts."""
    with open(sdrf_path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def _parse_mod_kv(raw: str) -> dict | None:
    """Parse a single ``comment[modification parameters]`` cell.

    Example: ``NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed``
    Returns ``{accession, name, fixed, position, target_amino_acid}``
    matching the QPX ``modification_param`` struct, or *None* if empty.
    """
    if not raw or not raw.strip():
        return None
    parts = {}
    for token in raw.split(";"):
        token = token.strip()
        if "=" in token:
            k, v = token.split("=", 1)
            parts[k.strip()] = v.strip()
    accession = parts.get("AC")
    if not accession:
        return None
    return {
        "accession": accession,
        "name": parts.get("NT"),
        "fixed": parts.get("MT", "").lower() == "fixed" if "MT" in parts else None,
        "position": parts.get("PP"),
        "target_amino_acid": parts.get("TA"),
    }


def extract_modification_params(sdrf_path: Path) -> list[dict]:
    """Extract unique modification parameters from all ``comment[modification parameters]`` columns.

    Uses raw TSV parsing (not DictReader) because the column name is
    repeated once per modification.
    """
    mods_seen: dict[str, dict] = {}  # keyed by accession to deduplicate
    with open(sdrf_path, encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        mod_indices = [i for i, h in enumerate(header) if h.strip().lower() == "comment[modification parameters]"]
        if not mod_indices:
            return []
        for row in reader:
            for idx in mod_indices:
                if idx < len(row):
                    mod = _parse_mod_kv(row[idx])
                    if mod and mod["accession"] not in mods_seen:
                        mods_seen[mod["accession"]] = mod
            break  # modifications are the same for all rows in a project
    return list(mods_seen.values())


def extract_organism(sdrf_rows: list[dict]) -> str:
    """Extract organism from SDRF rows (lowercase, underscored)."""
    for row in sdrf_rows:
        org = row.get("characteristics[organism]", "").strip()
        if org and org.lower() not in ("not available", "not applicable"):
            return org.lower().replace(" ", "_")
    return "unknown"


def extract_precursor_tolerance_da(sdrf_rows: list[dict]) -> float:
    """Extract precursor mass tolerance from SDRF and convert to Da.

    Parses ``comment[precursor mass tolerance]`` which looks like
    "10 ppm", "20 ppm", or "0.6 Da".  When the unit is ppm the value
    is converted to Da using a reference m/z of 1000 (mid-range peptide).

    Returns the tolerance in Da, or *DEFAULT_QC_TOLERANCE_DA* when the
    column is missing or cannot be parsed.
    """
    for row in sdrf_rows:
        raw = row.get("comment[precursor mass tolerance]", "").strip()
        if not raw:
            continue
        m = re.match(r"([\d.]+)\s*(ppm|da)", raw, re.IGNORECASE)
        if not m:
            continue
        value = float(m.group(1))
        unit = m.group(2).lower()
        if unit == "ppm":
            # Convert ppm to Da at reference m/z 1000
            return value * 1000.0 / 1e6
        return value
    return DEFAULT_QC_TOLERANCE_DA


def build_run_file_map(sdrf_rows: list[dict]) -> dict[str, dict]:
    """Map run_file_name (no ext) → {sample_accession, organism, instrument, ...}."""
    mapping = {}
    for row in sdrf_rows:
        data_file = row.get("comment[data file]", "").strip()
        if not data_file:
            continue
        run_file_name = _strip_raw_ext(data_file)
        source_name = row.get("source name", "").strip()
        organism = row.get("characteristics[organism]", "").strip()
        instrument_raw = row.get("comment[instrument]", "")
        # Parse instrument: "NT=Q Exactive HF-X;AC=MS:1002877"
        instrument = ""
        for part in instrument_raw.split(";"):
            if part.strip().startswith("NT="):
                instrument = part.strip()[3:]
        # Parse enzyme
        enzyme_raw = row.get("comment[cleavage agent details]", "")
        enzyme = ""
        for part in enzyme_raw.split(";"):
            if part.strip().startswith("NT="):
                enzyme = part.strip()[3:]
        # Parse dissociation
        dissoc_raw = row.get("comment[dissociation method]", "")
        dissoc = ""
        for part in dissoc_raw.split(";"):
            if part.strip().startswith("NT="):
                dissoc = part.strip()[3:]
        # Parse fraction
        fraction = row.get("comment[fraction identifier]", "").strip()
        # Parse label
        label = ""
        label_raw = row.get("comment[label]", "")
        for part in label_raw.split(";"):
            if part.strip().startswith("NT="):
                label = part.strip()[3:]
        # Parse assay name
        assay_name = row.get("assay name", "").strip()
        # Parse bio/tech replicate
        bio_rep = row.get("characteristics[biological replicate]", "").strip()
        tech_rep = row.get("comment[technical replicate]", "").strip()

        mapping[run_file_name] = {
            "sample_accession": source_name,
            "organism": organism,
            "instrument": instrument,
            "enzyme": enzyme,
            "dissociation_method": dissoc,
            "fraction": fraction,
            "label": label,
            "assay_name": assay_name,
            "data_file": data_file,
            "bio_replicate": int(bio_rep) if bio_rep.isdigit() else None,
            "tech_replicate": int(tech_rep) if tech_rep.isdigit() else None,
            "organism_part": row.get("characteristics[organism part]", "").strip(),
            "disease": row.get("characteristics[disease]", "").strip(),
            "cell_line": row.get("characteristics[cell line]", "").strip(),
            "cell_type": row.get("characteristics[cell type]", "").strip(),
            "sex": row.get("characteristics[sex]", "").strip(),
            "age": row.get("characteristics[age]", "").strip(),
            "developmental_stage": row.get("characteristics[developmental stage]", "").strip(),
            "ancestry": row.get("characteristics[ancestry category]", "").strip(),
        }
    return mapping


# ---------------------------------------------------------------------------
# Modification transform
# ---------------------------------------------------------------------------


def _transform_modifications(mod_col: pa.ChunkedArray) -> pa.Array:
    """Transform MSNet modifications to QPX modification format.

    MSNet: list<struct<name: string, positions: list<struct<position: string>>>>
    QPX:   list<struct<name: string, accession: string?, positions: list<struct<position: int32, amino_acid: string?, scores: list<score>?>>>>
    """
    py_mods = mod_col.to_pylist()
    transformed = []
    for row_mods in py_mods:
        if row_mods is None:
            transformed.append(None)
            continue
        new_mods = []
        for mod in row_mods:
            if mod is None:
                continue
            name = mod.get("name", "")
            positions = mod.get("positions", [])
            new_positions = []
            if positions:
                for pos in positions:
                    pos_str = pos.get("position", "0") if pos else "0"
                    try:
                        pos_int = int(pos_str)
                    except (ValueError, TypeError):
                        pos_int = 0
                    new_positions.append(
                        {
                            "position": pos_int,
                            "amino_acid": None,
                            "scores": None,
                        }
                    )
            new_mods.append(
                {
                    "name": name,
                    "accession": None,
                    "positions": new_positions,
                }
            )
        transformed.append(new_mods)

    # Build QPX modification Arrow type
    score_type = pa.struct(
        [
            pa.field("score_name", pa.string()),
            pa.field("score_value", pa.float64()),
            pa.field("higher_better", pa.bool_(), nullable=True),
        ]
    )
    pos_type = pa.struct(
        [
            pa.field("position", pa.int32()),
            pa.field("amino_acid", pa.string(), nullable=True),
            pa.field("scores", pa.list_(score_type), nullable=True),
        ]
    )
    mod_type = pa.struct(
        [
            pa.field("name", pa.string()),
            pa.field("accession", pa.string(), nullable=True),
            pa.field("positions", pa.list_(pos_type)),
        ]
    )
    return pa.array(transformed, type=pa.list_(mod_type))


# ---------------------------------------------------------------------------
# CV params transform
# ---------------------------------------------------------------------------


def _transform_cv_params(cv_col: pa.ChunkedArray) -> pa.Array:
    """Transform MSNet cv_params struct to QPX list<cv_param>.

    MSNet: struct<Instrument, Fragmentation Method, Collision Energy, Mass Analyzer>
    QPX:   list<struct<cv_name: string, cv_value: string>>
    """
    py_cv = cv_col.to_pylist()
    transformed = []
    for row in py_cv:
        if row is None:
            transformed.append(None)
            continue
        pairs = []
        for k, v in row.items():
            if v is not None:
                pairs.append({"cv_name": k, "cv_value": str(v)})
        transformed.append(pairs)

    cv_param_type = pa.struct(
        [
            pa.field("cv_name", pa.string()),
            pa.field("cv_value", pa.string()),
        ]
    )
    return pa.array(transformed, type=pa.list_(cv_param_type))


# ---------------------------------------------------------------------------
# Additional scores
# ---------------------------------------------------------------------------


def _build_additional_scores(table: pa.Table) -> pa.Array:
    """Pack global_qvalue and consensus_support into list<score>."""
    gq = table.column("global_qvalue").to_pylist()
    cs = table.column("consensus_support").to_pylist()

    scores_list = []
    for g, c in zip(gq, cs, strict=True):
        row_scores = []
        if g is not None:
            row_scores.append(
                {
                    "score_name": "global_qvalue",
                    "score_value": float(g),
                    "higher_better": None,
                }
            )
        if c is not None:
            row_scores.append(
                {
                    "score_name": "consensus_support",
                    "score_value": float(c),
                    "higher_better": None,
                }
            )
        scores_list.append(row_scores if row_scores else None)

    score_type = pa.struct(
        [
            pa.field("score_name", pa.string()),
            pa.field("score_value", pa.float64()),
            pa.field("higher_better", pa.bool_(), nullable=True),
        ]
    )
    return pa.array(scores_list, type=pa.list_(score_type))


# ---------------------------------------------------------------------------
# PSM conversion
# ---------------------------------------------------------------------------


def _col(table: pa.Table, *names: str) -> pa.ChunkedArray | None:
    """Get column by first matching name, or None if absent."""
    for name in names:
        if name in table.schema.names:
            return table.column(name)
    return None


# Columns required by convert_psm_table
_NEEDED_COLUMNS = {
    "sequence",
    "peptidoform",
    "precursor_charge",
    "posterior_error_probability",
    "retention_time",
    "reference_file_name",
    "scan",
    "modifications",
    "cv_params",
    "global_qvalue",
    "consensus_support",
    "mz_array",
    "intensity_array",
    "ion_type_array",
    "charge_array",
    # Schema 4 alternatives
    "observed_mz",
    "exp_mass_to_charge",
    "calculated_mz",
    "cal_mass_to_charge",
    # Optional
    "ion_mobility",
    "Luciphor_global_flr",
    "Luciphor_local_flr",
}


def _needed_column_names(schema: pa.Schema) -> list[str]:
    """Return column names (unique) that convert_psm_table needs."""
    return sorted({f.name for f in schema} & _NEEDED_COLUMNS)


def _deduplicate_scan(table: pa.Table) -> pa.ChunkedArray:
    """Handle duplicate 'scan' columns and string scan values.

    PXD000561 has two scan columns (string + int32).
    PXD014877-Dorea has scan as string type.
    """
    indices = [i for i, f in enumerate(table.schema) if f.name == "scan"]
    # Prefer int32 column if multiple exist
    for idx in indices:
        if pa.types.is_integer(table.schema.field(idx).type):
            return table.column(idx)
    # Fallback: first scan column (may be string, handled in caller)
    return table.column(indices[0])


def convert_psm_table(msnet_table: pa.Table) -> pa.Table:
    """Convert an MSNet parquet table to QPX PSM schema.

    Handles 4 known MSNet schema variants:
      1. Standard (40 projects): exp_mass_to_charge / cal_mass_to_charge
      2. PXD000561: duplicate 'scan' column
      3. PXD000612: extra Luciphor columns
      4. PXD040266: observed_mz / calculated_mz + ion_mobility
    """
    n = msnet_table.num_rows
    col_names = set(msnet_table.schema.names)

    # Direct columns
    sequence = msnet_table.column("sequence")
    peptidoform = msnet_table.column("peptidoform")
    charge = pa.compute.cast(msnet_table.column("precursor_charge"), pa.int16())
    pep = msnet_table.column("posterior_error_probability")
    rt = pa.compute.cast(msnet_table.column("retention_time"), pa.float32())
    run_file_name = msnet_table.column("reference_file_name")

    # Schema 4 uses observed_mz/calculated_mz; others use exp_/cal_mass_to_charge
    obs_raw = _col(msnet_table, "observed_mz", "exp_mass_to_charge")
    calc_raw = _col(msnet_table, "calculated_mz", "cal_mass_to_charge")
    if obs_raw is None or calc_raw is None:
        raise ValueError("Missing required m/z columns: need observed_mz/exp_mass_to_charge and calculated_mz/cal_mass_to_charge")
    observed_mz = pa.compute.cast(obs_raw, pa.float32())
    calculated_mz = pa.compute.cast(calc_raw, pa.float32())

    # Compute mass_error_ppm = 1e6 * (obs - calc) / calc
    obs_f64 = pa.compute.cast(observed_mz, pa.float64())
    calc_f64 = pa.compute.cast(calculated_mz, pa.float64())
    diff = pa.compute.subtract(obs_f64, calc_f64)
    ratio = pa.compute.divide(diff, calc_f64)
    mass_error_ppm = pa.compute.cast(
        pa.compute.multiply(ratio, pa.scalar(1e6, pa.float64())),
        pa.float32(),
    )

    # is_decoy = False for all (MSNet has targets only)
    is_decoy = pa.nulls(n, type=pa.bool_()).fill_null(False)

    # Scan: handle duplicate column, string type, int32 → list<int32>
    scan_col = _deduplicate_scan(msnet_table)
    scan_raw = scan_col.to_pylist()
    scan_list = []
    for s in scan_raw:
        if s is None:
            scan_list.append([0])
        elif isinstance(s, str):
            try:
                scan_list.append([int(s)])
            except ValueError:
                scan_list.append([0])
        else:
            scan_list.append([s])
    scan = pa.array(scan_list, type=pa.list_(pa.int32()))

    # Structured transforms
    modifications = _transform_modifications(msnet_table.column("modifications"))
    cv_params = _transform_cv_params(msnet_table.column("cv_params"))
    additional_scores = _build_additional_scores(msnet_table)

    # Extra scores from Schema 3 (Luciphor)
    if "Luciphor_global_flr" in col_names or "Luciphor_local_flr" in col_names:
        extra = additional_scores.to_pylist()
        gflr = _col(msnet_table, "Luciphor_global_flr")
        lflr = _col(msnet_table, "Luciphor_local_flr")
        gflr_list = gflr.to_pylist() if gflr is not None else [None] * n
        lflr_list = lflr.to_pylist() if lflr is not None else [None] * n
        for i in range(n):
            row = extra[i] or []
            if gflr_list[i] is not None:
                row.append(
                    {
                        "score_name": "Luciphor_global_flr",
                        "score_value": float(gflr_list[i]),
                        "higher_better": None,
                    }
                )
            if lflr_list[i] is not None:
                row.append(
                    {
                        "score_name": "Luciphor_local_flr",
                        "score_value": float(lflr_list[i]),
                        "higher_better": None,
                    }
                )
            extra[i] = row if row else None
        score_type = pa.struct(
            [
                pa.field("score_name", pa.string()),
                pa.field("score_value", pa.float64()),
                pa.field("higher_better", pa.bool_(), nullable=True),
            ]
        )
        additional_scores = pa.array(extra, type=pa.list_(score_type))

    # Spectral arrays
    mz_array = msnet_table.column("mz_array")
    intensity_array = msnet_table.column("intensity_array")
    ion_type_array = msnet_table.column("ion_type_array")
    charge_array = msnet_table.column("charge_array")

    # ion_mobility: Schema 4 (PXD040266) has it natively
    ion_mob_col = _col(msnet_table, "ion_mobility")
    ion_mobility = pa.compute.cast(ion_mob_col, pa.float32()) if ion_mob_col is not None else pa.nulls(n, type=pa.float32())

    # Nullable columns not present in MSNet — fill with nulls
    null_f32 = pa.nulls(n, type=pa.float32())
    null_i16 = pa.nulls(n, type=pa.int16())
    null_str_list = pa.nulls(n, type=pa.list_(pa.string()))
    null_f32_list = pa.nulls(n, type=pa.list_(pa.float32()))

    cross_link_type = pa.struct(
        [
            pa.field("xl_type", pa.string()),
            pa.field("partner_sequence", pa.string(), nullable=True),
            pa.field("partner_peptidoform", pa.string(), nullable=True),
            pa.field("donor_position", pa.int32()),
            pa.field("acceptor_position", pa.int32(), nullable=True),
            pa.field("linker_name", pa.string()),
            pa.field("linker_accession", pa.string(), nullable=True),
            pa.field("linker_mass", pa.float64()),
        ]
    )
    null_cross_links = pa.nulls(n, type=pa.list_(cross_link_type))

    # Build QPX PSM table (column order follows psm.yaml)
    qpx_table = pa.table(
        {
            "sequence": sequence,
            "peptidoform": peptidoform,
            "modifications": modifications,
            "charge": charge,
            "posterior_error_probability": pep,
            "is_decoy": is_decoy,
            "calculated_mz": calculated_mz,
            "observed_mz": observed_mz,
            "mass_error_ppm": mass_error_ppm,
            "additional_scores": additional_scores,
            "predicted_rt": null_f32,
            "run_file_name": run_file_name,
            "cv_params": cv_params,
            "scan": scan,
            "rt": rt,
            "ion_mobility": ion_mobility,
            "missed_cleavages": null_i16,
            "protein_accessions": null_str_list,
            "cross_links": null_cross_links,
            "mz_array": mz_array,
            "intensity_array": intensity_array,
            "charge_array": charge_array,
            "ion_type_array": ion_type_array,
            "ion_mobility_array": null_f32_list,
        }
    )

    return qpx_table


# ---------------------------------------------------------------------------
# Sample table from SDRF
# ---------------------------------------------------------------------------


def build_sample_table(run_file_map: dict[str, dict]) -> pa.Table:
    """Build QPX sample.parquet from SDRF run_file_map."""
    seen = {}
    for info in run_file_map.values():
        sa = info["sample_accession"]
        if sa and sa not in seen:
            seen[sa] = info

    rows = []
    for sa, info in seen.items():
        row = {
            "sample_accession": sa,
            "organism": info["organism"] or "unknown",
            "organism_part": info["organism_part"] or "not available",
        }
        # Optional fields
        for field, key in [
            ("disease", "disease"),
            ("cell_line", "cell_line"),
            ("cell_type", "cell_type"),
            ("sex", "sex"),
            ("age", "age"),
            ("developmental_stage", "developmental_stage"),
            ("ancestry", "ancestry"),
        ]:
            val = info.get(key, "")
            if val and val.lower() not in ("not applicable", "not available", ""):
                row[field] = val
        rows.append(row)

    if not rows:
        return None

    from qpx.core.data import SampleSchema

    # Determine which optional columns are present
    all_keys = set()
    for r in rows:
        all_keys.update(r.keys())
    schema = SampleSchema.get_arrow_schema()
    fields = []
    for f in schema:
        if f.name in all_keys:
            fields.append(f)

    # Build arrays
    arrays = {}
    for f in fields:
        vals = [r.get(f.name) for r in rows]
        arrays[f.name] = pa.array(vals, type=f.type)

    return pa.table(arrays)


# ---------------------------------------------------------------------------
# Run table from SDRF
# ---------------------------------------------------------------------------


def build_run_table(run_file_map: dict[str, dict]) -> pa.Table:
    """Build QPX run.parquet from SDRF run_file_map."""
    rows = []
    for _key, info in run_file_map.items():
        # Use stored run_file_name (keys may be namespaced to avoid collisions)
        rfn = info.get("_run_file_name_key", _key)
        sample_channels = [
            {
                "sample_accession": info["sample_accession"],
                "label": info["label"] or "label free sample",
                "biological_replicate": info["bio_replicate"],
                "technical_replicate": info["tech_replicate"],
            }
        ]
        row = {
            "run_accession": info["assay_name"] or rfn,
            "run_file_name": rfn,
            "file_name": info["data_file"],
            "samples": sample_channels,
            "fraction": info["fraction"] or None,
            "instrument": info["instrument"] or None,
            "enzymes": [info["enzyme"]] if info["enzyme"] else None,
            "dissociation_method": info["dissociation_method"] or None,
            "modification_parameters": info.get("modification_parameters"),
        }
        rows.append(row)

    if not rows:
        return None

    from qpx.core.data import RunSchema

    schema = RunSchema.get_arrow_schema()
    # Build each column
    arrays = {}
    for f in schema:
        vals = [r.get(f.name) for r in rows]
        arrays[f.name] = pa.array(vals, type=f.type)

    return pa.table(arrays)


# ---------------------------------------------------------------------------
# Dataset table
# ---------------------------------------------------------------------------


def build_dataset_table(project_accessions: list[str]) -> pa.Table:
    """Build QPX dataset.parquet with project-level metadata."""
    rows = []
    for pa_str in project_accessions:
        rows.append(
            {
                "project_accession": pa_str,
                "project_title": None,
                "project_description": None,
                "pubmed_id": None,
                "software_name": "quantms/MSNet",
                "software_version": None,
                "creation_date": datetime.datetime.now().isoformat(),
            }
        )

    from qpx.core.data import DatasetSchema

    schema = DatasetSchema.get_arrow_schema()
    arrays = {}
    for f in schema:
        vals = [r.get(f.name) for r in rows]
        try:
            arrays[f.name] = pa.array(vals, type=f.type)
        except (pa.ArrowInvalid, pa.ArrowNotImplementedError):
            _log.debug("Skipping dataset column '%s' (unsupported type)", f.name)

    return pa.table(arrays)


# ---------------------------------------------------------------------------
# Main conversion
# ---------------------------------------------------------------------------


def find_parquet_files(proj_dir: Path) -> list[Path]:
    """Find all parquet files in psm/ or root."""
    psm_dir = proj_dir / "psm"
    if psm_dir.is_dir():
        files = list(psm_dir.glob("*.parquet"))
        if files:
            return files
    return list(proj_dir.glob("*.parquet"))


def find_sdrf(proj_dir: Path) -> Path | None:
    """Find SDRF file."""
    sdrf_dir = proj_dir / "sdrf"
    if sdrf_dir.is_dir():
        files = list(sdrf_dir.glob("*.sdrf.tsv"))
        if files:
            return files[0]
    return None


_write_lock = threading.Lock()

# ---------------------------------------------------------------------------
# QC: mass error check
# ---------------------------------------------------------------------------

DEFAULT_QC_TOLERANCE_DA = 0.5  # fallback when SDRF lacks precursor tolerance
QC_TOLERANCE_MULTIPLIER = 10  # flag PSMs exceeding 10x the search tolerance
MAX_OUTLIER_FRACTION = 0.20  # skip project if >20% PSMs exceed threshold


def _check_mass_error_qc(
    output_dir: Path,
    organism: str,
    run_file_names: list[str],
    tolerance_da: float = DEFAULT_QC_TOLERANCE_DA,
) -> dict:
    """Check fraction of PSMs with |mass_error_da| > tolerance.

    The tolerance is derived from the SDRF precursor mass tolerance
    multiplied by *QC_TOLERANCE_MULTIPLIER*.  This adapts the check
    to the instrument resolution of each project.

    Returns dict with keys: total_psms, outlier_count, outlier_pct,
    tolerance_da, passed.
    """
    import duckdb

    qc_limit = tolerance_da * QC_TOLERANCE_MULTIPLIER

    psm_dir = output_dir / "psm"
    safe_path = escape_path(str(psm_dir))

    if not run_file_names:
        return {"total_psms": 0, "outlier_count": 0, "outlier_pct": 0.0, "tolerance_da": tolerance_da, "passed": True}

    con = duckdb.connect()
    try:
        sql = sql_build(
            "SELECT count(*) AS total, "
            "count_if(abs(observed_mz - calculated_mz) > ?) AS outliers "
            "FROM parquet_scan('$path/**/*.parquet', hive_partitioning=true) "
            "WHERE organism = ? "
            "AND list_contains(?, run_file_name)",
            path=safe_path,
        )
        row = con.execute(sql, [qc_limit, organism, run_file_names]).fetchone()
        total = row[0] or 0
        outliers = row[1] or 0
    except Exception as e:
        _log.warning("QC check failed: %s", e)
        return {"total_psms": 0, "outlier_count": 0, "outlier_pct": 0.0, "tolerance_da": tolerance_da, "passed": True}
    finally:
        con.close()

    pct = (outliers / total * 100) if total > 0 else 0.0
    passed = (outliers / total) <= MAX_OUTLIER_FRACTION if total > 0 else True
    return {
        "total_psms": total,
        "outlier_count": outliers,
        "outlier_pct": round(pct, 2),
        "tolerance_da": tolerance_da,
        "passed": passed,
    }


def _write_partitioned_streaming(
    table: pa.Table,
    output_dir: Path,
    partition_cols: list[str],
):
    """Write an Arrow table as Hive-partitioned Parquet (thread-safe).

    Each call uses a unique basename_template so that multiple batches
    writing to the same partition do not overwrite each other.
    """
    # Encoding parity with the single-file writers: same parquet_write_options
    # path (tuned ZSTD level, selective BYTE_STREAM_SPLIT + dictionary, v2.6),
    # a sane row-group size (streamed batches otherwise fragment into ~80-row
    # groups and defeat ZSTD), and the qpx footer KV metadata so partitioned
    # PSM files keep qpx_version / file_type / uuid.
    table = _stamp_footer_metadata(table, "zstd", {b"file_type": b"psm"})
    part_schema = pa.schema([table.schema.field(c) for c in partition_cols])
    partitioning = pds.partitioning(part_schema, flavor="hive")
    opts = parquet_write_options(table.schema, "zstd")
    file_options = pds.ParquetFileFormat().make_write_options(**opts)
    basename = f"part-{uuid.uuid4().hex[:12]}-{{i}}.parquet"
    with _write_lock:
        pds.write_dataset(
            table,
            str(output_dir),
            format="parquet",
            partitioning=partitioning,
            existing_data_behavior="overwrite_or_ignore",
            basename_template=basename,
            file_options=file_options,
            max_partitions=4096,
            min_rows_per_group=_PARTITIONED_MIN_ROWS_PER_GROUP,
            max_rows_per_group=_PARTITIONED_MAX_ROWS_PER_GROUP,
        )


def convert_project(
    proj_name: str,
    output_dir: Path,
    msnet_dir: Path = Path(DEFAULT_MSNET_DIR),
) -> tuple[dict[str, dict] | None, dict | None]:
    """Convert a single MSNet project to QPX format (streaming by row group).

    Returns (run_file_map, qc_metrics) where qc_metrics contains mass-error
    quality metrics, or (None, None) on failure.
    """
    proj_dir = msnet_dir / proj_name
    if not proj_dir.is_dir():
        _log.error("Project directory not found: %s", proj_dir)
        return None, None

    _log.info("Converting %s...", proj_name)

    # Find SDRF
    sdrf_path = find_sdrf(proj_dir)
    if not sdrf_path:
        _log.error("No SDRF found for %s", proj_name)
        return None, None

    sdrf_rows = parse_sdrf(sdrf_path)
    organism = extract_organism(sdrf_rows)
    precursor_tol_da = extract_precursor_tolerance_da(sdrf_rows)
    mod_params = extract_modification_params(sdrf_path)
    run_file_map = build_run_file_map(sdrf_rows)
    # Attach modification params to each run (same for all runs in a project)
    for info in run_file_map.values():
        info["modification_parameters"] = mod_params or None
    _log.info(
        "  Organism: %s, SDRF runs: %d, precursor tol: %.4f Da, mods: %d",
        organism,
        len(run_file_map),
        precursor_tol_da,
        len(mod_params),
    )

    # Find and convert parquet files
    pq_files = find_parquet_files(proj_dir)
    if not pq_files:
        _log.error("No parquet files found for %s", proj_name)
        return None, None

    psm_dir = output_dir / "psm"
    total_rows = 0

    for pf in pq_files:
        pf_meta = pq.read_metadata(pf)
        n_rg = pf_meta.num_row_groups
        size_mb = pf.stat().st_size / (1024 * 1024)
        _log.info("  %s (%.0f MB, %d row groups)", pf.name, size_mb, n_rg)

        pf_reader = pq.ParquetFile(pf)
        col_names = _needed_column_names(pf_reader.schema_arrow)

        # Large files: smaller batch to avoid OOM; others: full row group
        batch_size = 100_000 if size_mb > 2000 else 1_048_576
        batch_num = 0

        for batch in pf_reader.iter_batches(batch_size=batch_size, columns=col_names):
            rg_table = pa.Table.from_batches([batch])
            qpx_batch = convert_psm_table(rg_table)

            # Add organism column
            n = qpx_batch.num_rows
            qpx_batch = qpx_batch.append_column("organism", pa.array([organism] * n, type=pa.string()))

            # Stream write immediately
            _write_partitioned_streaming(qpx_batch, psm_dir, ["organism", "run_file_name"])
            total_rows += n
            batch_num += 1

            if batch_num % 5 == 0 or size_mb > 2000:
                _log.info("    batch %d: %d rows (total %d)", batch_num, n, total_rows)

    _log.info("  Total PSM rows: %d", total_rows)

    # QC: check mass error outlier rate (adaptive to SDRF precursor tolerance)
    qc = _check_mass_error_qc(output_dir, organism, list(run_file_map.keys()), tolerance_da=precursor_tol_da)
    qc["project"] = proj_name
    qc["organism"] = organism
    qc_limit = precursor_tol_da * QC_TOLERANCE_MULTIPLIER
    if not qc["passed"]:
        _log.error(
            "QC FAILED %s: %.1f%% PSMs with |mass_error| > %.4f Da (precursor_tol=%.4f Da x%d, threshold: %s%%)",
            proj_name,
            qc["outlier_pct"],
            qc_limit,
            precursor_tol_da,
            QC_TOLERANCE_MULTIPLIER,
            MAX_OUTLIER_FRACTION * 100,
        )
        return None, qc

    _log.info(
        "  QC passed: %.1f%% outliers (%d/%d), limit=%.4f Da",
        qc["outlier_pct"],
        qc["outlier_count"],
        qc["total_psms"],
        qc_limit,
    )
    return run_file_map, qc


def _write_metadata(
    output_dir: Path,
    all_run_file_maps: dict[str, dict],
    project_accessions: list[str],
):
    """Write sample, run, and dataset parquet files."""
    if all_run_file_maps:
        sample_table = build_sample_table(all_run_file_maps)
        if sample_table is not None:
            pq.write_table(sample_table, str(output_dir / "msnet.sample.parquet"), compression="zstd")
            _log.info("Wrote sample.parquet: %d rows", sample_table.num_rows)

        run_table = build_run_table(all_run_file_maps)
        if run_table is not None:
            pq.write_table(run_table, str(output_dir / "msnet.run.parquet"), compression="zstd")
            _log.info("Wrote run.parquet: %d rows", run_table.num_rows)

    if project_accessions:
        ds_table = build_dataset_table(list(set(project_accessions)))
        if ds_table is not None:
            pq.write_table(ds_table, str(output_dir / "msnet.dataset.parquet"), compression="zstd")
            _log.info("Wrote dataset.parquet: %d projects", ds_table.num_rows)


def _write_provenance(output_dir: Path):
    """Write provenance.parquet with pipeline steps."""
    provenance_rows = [
        {
            "step_order": 1,
            "step_category": "database_search",
            "step_name": "sequence_search",
            "tool_name": "Multiple (project-dependent)",
            "tool_version": None,
            "tool_uri": None,
            "parameters": None,
            "config": None,
            "output_views": ["psm"],
        },
        {
            "step_order": 2,
            "step_category": "psm_rescoring",
            "step_name": "msnet_rescoring",
            "tool_name": "MSNet",
            "tool_version": None,
            "tool_uri": "https://github.com/bigbio/MSNet",
            "parameters": None,
            "config": None,
            "output_views": ["psm"],
        },
        {
            "step_order": 3,
            "step_category": "conversion",
            "step_name": "msnet_to_qpx",
            "tool_name": "QPX convert_msnet.py",
            "tool_version": __version__,
            "tool_uri": None,
            "parameters": [
                {"key": "partition_by", "value": "organism, run_file_name"},
                {"key": "compression", "value": "zstd"},
            ],
            "config": None,
            "output_views": ["psm"],
        },
    ]
    with ProvenanceWriter(output_dir / "msnet.provenance.parquet", creator="qpx", compression="zstd") as writer:
        writer.write_batch(provenance_rows)
    _log.info("Wrote provenance.parquet: %d steps", len(provenance_rows))


def _write_ontology(output_dir: Path):
    """Write ontology.parquet with field mappings."""
    entries = field_ontology_entries(
        view="psm",
        resolved_mappings={
            "sequence": "sequence",
            "peptidoform": "peptidoform",
            "charge": "precursor_charge",
            "posterior_error_probability": "posterior_error_probability",
            "observed_mz": "exp_mass_to_charge",
            "calculated_mz": "cal_mass_to_charge",
            "rt": "retention_time",
            "is_decoy": "is_decoy",
            "predicted_rt": "predicted_rt",
        },
        tool_name="MSNet",
    )
    with OntologyWriter(output_dir / "msnet.ontology.parquet", creator="qpx", compression="zstd") as writer:
        writer.write_batch(entries)
    _log.info("Wrote ontology.parquet: %d entries", len(entries))


def _write_qc_report(output_dir: Path, all_qc_metrics: list[dict]):
    """Write QC report TSV summarising mass-error checks per project."""
    if not all_qc_metrics:
        return
    qc_path = output_dir / "msnet.qc_report.tsv"
    with open(qc_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["project", "organism", "tolerance_da", "total_psms", "outlier_count", "outlier_pct", "status"],
            delimiter="\t",
        )
        writer.writeheader()
        for qc in sorted(all_qc_metrics, key=lambda x: x.get("project", "")):
            writer.writerow({k: qc.get(k) for k in writer.fieldnames})
    _log.info("Wrote QC report: %s (%d entries)", qc_path, len(all_qc_metrics))


# ---------------------------------------------------------------------------
# Parallel conversion orchestration
# ---------------------------------------------------------------------------


def _run_conversions(
    projects: list[str],
    msnet_dir: Path,
    output_dir: Path,
    proj_sizes: dict[str, int],
) -> tuple[int, int, int, list[str], dict[str, dict], list[dict]]:
    """Run project conversions (small in parallel, large sequentially).

    Returns (ok, qc_skipped, fail, project_accessions,
             all_run_file_maps, all_qc_metrics).
    """
    SIZE_THRESHOLD = 500 * 1024 * 1024  # 500 MB
    small_projs = [p for p in projects if proj_sizes[p] < SIZE_THRESHOLD]
    large_projs = [p for p in projects if proj_sizes[p] >= SIZE_THRESHOLD]
    _log.info(
        "Converting %d projects (%d small, %d large >500MB) to %s",
        len(projects),
        len(small_projs),
        len(large_projs),
        output_dir,
    )

    ok = 0
    fail = 0
    qc_skipped = 0
    project_accessions: list[str] = []
    all_run_file_maps: dict[str, dict] = {}
    all_qc_metrics: list[dict] = []
    lock = threading.Lock()
    counter = [0]

    def _worker(proj: str) -> tuple[str, dict[str, dict] | None, dict | None]:
        with lock:
            counter[0] += 1
            idx = counter[0]
        size_mb = proj_sizes[proj] / (1024 * 1024)
        _log.info("[%d/%d] %s (%.0f MB)", idx, len(projects), proj, size_mb)
        try:
            rfm, qc = convert_project(proj, output_dir, msnet_dir)
            return proj, rfm, qc
        except (OSError, ValueError, pa.ArrowException) as e:
            _log.error("FAILED %s: %s", proj, e, exc_info=True)
            return proj, None, None

    def _collect(proj: str, rfm: dict[str, dict] | None, qc: dict | None):
        nonlocal ok, fail, qc_skipped
        if qc is not None:
            qc["status"] = "passed" if qc["passed"] else "qc_failed"
            all_qc_metrics.append(qc)
        if rfm is not None:
            ok += 1
            acc = proj.split("-")[0] if "-" in proj else proj
            project_accessions.append(acc)
            # Namespace keys by project to avoid cross-project collisions
            for rfn, info in rfm.items():
                info["_run_file_name_key"] = rfn
                all_run_file_maps[f"{proj}::{rfn}"] = info
        elif qc is not None and not qc["passed"]:
            qc_skipped += 1
        else:
            fail += 1

    if small_projs:
        workers = min(8, len(small_projs))
        _log.info("Phase 1: %d small projects with %d workers", len(small_projs), workers)
        with ThreadPoolExecutor(max_workers=workers, thread_name_prefix="conv") as pool:
            futures = {pool.submit(_worker, p): p for p in small_projs}
            for future in as_completed(futures):
                proj, rfm, qc = future.result()
                _collect(proj, rfm, qc)

    if large_projs:
        _log.info("Phase 2: %d large projects sequentially", len(large_projs))
        for proj in large_projs:
            proj, rfm, qc = _worker(proj)
            _collect(proj, rfm, qc)

    return ok, qc_skipped, fail, project_accessions, all_run_file_maps, all_qc_metrics


def main():
    parser = argparse.ArgumentParser(description="Convert MSNet projects to QPX dataset")
    parser.add_argument("--projects", nargs="+", help="Project names to convert")
    parser.add_argument("--all", action="store_true", help="Convert all local projects")
    parser.add_argument("--input", type=str, default=DEFAULT_MSNET_DIR, help="MSNet data directory")
    parser.add_argument("--output", type=str, default=DEFAULT_OUTPUT_DIR, help="Output directory")
    args = parser.parse_args()

    msnet_dir = Path(args.input)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.all:
        projects = sorted([d.name for d in msnet_dir.iterdir() if d.is_dir() and find_sdrf(d) is not None])
    elif args.projects:
        projects = args.projects
    else:
        parser.error("Specify --projects or --all")

    proj_sizes = {}
    for proj in projects:
        proj_dir = msnet_dir / proj
        proj_sizes[proj] = sum(f.stat().st_size for f in proj_dir.rglob("*.parquet"))

    ok, qc_skipped, fail, project_accessions, all_run_file_maps, all_qc_metrics = _run_conversions(
        projects, msnet_dir, output_dir, proj_sizes
    )

    _log.info("Writing unified metadata files...")
    _write_metadata(output_dir, all_run_file_maps, project_accessions)
    _write_provenance(output_dir)
    _write_ontology(output_dir)
    _write_qc_report(output_dir, all_qc_metrics)

    _log.info("\nDone: %d converted, %d QC-skipped, %d failed", ok, qc_skipped, fail)


if __name__ == "__main__":
    main()
