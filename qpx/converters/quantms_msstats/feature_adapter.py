"""Convert QuantMS ``*_msstats_in.csv`` tables to QPX Features."""

from __future__ import annotations

import logging
import re
from collections.abc import Iterable

import pandas as pd

from qpx.converters.base import BaseConverter
from qpx.converters.channel_labels import (
    experiment_type_from_labels,
    normalize_label,
    read_sdrf_labels,
    resolve_channel_labels,
)
from qpx.converters.openms_consensus.feature_adapter import to_modifications, to_proforma
from qpx.core.files import run_file_stem
from qpx.core.sdrf import validate_sdrf_data_files
from qpx.core.sql import escape_path, sql_build, validate_identifier, validate_table
from qpx.writers.feature import FeatureWriter

logger = logging.getLogger(__name__)

FEATURE_IDENTITY_COMPOSITE = ("run_file_name", "peptidoform", "charge", "rt", "scan")

_RUN_ALIAS_COLUMNS = (
    "comment[data file]",
    "comment[file uri]",
    "comment[associated file uri]",
    "assay name",
)
_RUN_EXTENSION = re.compile(
    r"^(.*)\.(?:mzml(?:\.gz)?|mzxml|raw|d|wiff|mgf|dia)(?:$|[._\s])",
    re.IGNORECASE,
)
_DECOY_PREFIXES = ("decoy_", "rev_", "reverse_")


def _run_alias_key(value: object) -> str:
    """Return the normalized run key used to match MSstats and SDRF values."""
    basename = str(value).strip().replace("\\", "/").rsplit("/", maxsplit=1)[-1].lower()
    match = _RUN_EXTENSION.match(basename)
    if match:
        return match.group(1)
    return run_file_stem(basename).lower()


def _split_proteins(values: Iterable[str]) -> list[str]:
    """Return distinct protein accessions in deterministic source order."""
    accessions: list[str] = []
    seen: set[str] = set()
    for raw in values:
        for value in str(raw).split(";"):
            accession = value.strip()
            if accession and accession not in seen:
                seen.add(accession)
                accessions.append(accession)
    return accessions


def _is_decoy(accessions: list[str]) -> bool:
    """Return true only when every reported protein is a decoy accession."""
    return bool(accessions) and all(accession.lower().startswith(_DECOY_PREFIXES) for accession in accessions)


def _as_channel_position(value: object) -> int | None:
    """Parse a one-based integer MSstats channel position."""
    try:
        numeric = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    position = int(numeric)
    return position if numeric == position and position > 0 else None


def _parse_peptide(raw: str) -> tuple[str, str, list[dict] | None]:
    """Convert an OpenMS peptide string into QPX sequence/PTM fields."""
    import pyopenms as oms

    try:
        aa_sequence = getattr(oms, "AASequence").fromString(raw)
    except RuntimeError as exc:
        raise ValueError(f"Invalid QuantMS PeptideSequence value: {raw!r}") from exc
    sequence = aa_sequence.toUnmodifiedString()
    peptidoform = to_proforma(aa_sequence)
    if not sequence or not peptidoform:
        raise ValueError(f"Invalid QuantMS PeptideSequence value: {raw!r}")
    return sequence, peptidoform, to_modifications(aa_sequence) or None


class QuantmsMsstatsFeatureAdapter(BaseConverter):
    """Convert one QuantMS MSstats input table into ``feature.parquet``."""

    _COLUMN_ALIASES = {
        "protein": ("ProteinName",),
        "peptide": ("PeptideSequence",),
        "charge": ("Charge", "PrecursorCharge"),
        "intensity": ("Intensity",),
        "run": ("Run",),
        "reference": ("Reference",),
        "channel": ("Channel",),
        "rt": ("RetentionTime",),
    }

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._columns: dict[str, str | None] = {}
        self._resolved: dict[str, str] = {}
        self._modifications: dict[str, list[dict] | None] = {}

    def convert(
        self,
        msstats_path: str,
        sdrf_path: str,
        output_path: str,
        chunksize: int = 50_000,
        creator: str = "quantms-msstats",
    ) -> None:
        """Validate and convert a QuantMS MSstats table."""
        self._load_msstats(msstats_path)
        self._resolve_columns()
        aliases, valid_pairs = self._load_sdrf_context(sdrf_path)
        self._register_frame("_run_aliases", aliases)
        self._register_frame("_valid_run_labels", valid_pairs)
        channel_map = self._build_channel_map(sdrf_path)
        if channel_map is not None:
            self._register_frame("_channel_map", channel_map)
        self._register_peptides()
        self._create_normalized_table(channel_map is not None)
        self._validate_normalized_rows()
        self._validate_run_labels()
        self._validate_feature_groups()
        self._write_features(output_path, chunksize, creator)

    def _load_msstats(self, path: str) -> None:
        safe_path = escape_path(path)
        self._conn.execute(
            sql_build(
                """
                CREATE OR REPLACE VIEW _msstats_input AS
                SELECT * FROM read_csv_auto(
                    '$path',
                    header=true,
                    all_varchar=true,
                    delim=',',
                    quote='"',
                    escape='"'
                )
                """,
                path=safe_path,
            )
        )
        count = self._conn.execute("SELECT count(*) FROM _msstats_input").fetchone()[0]
        if count == 0:
            raise ValueError("QuantMS MSstats input contains no data rows")
        self.logger.info("Loaded %d QuantMS MSstats rows from %s", count, path)

    def _resolve_columns(self) -> None:
        available = {
            row[0]
            for row in self._conn.execute(
                "SELECT column_name FROM information_schema.columns WHERE table_name='_msstats_input'"
            ).fetchall()
        }
        required = {"protein", "peptide", "charge", "intensity"}
        for field, aliases in self._COLUMN_ALIASES.items():
            self._columns[field] = next((name for name in aliases if name in available), None)
            if field in required and self._columns[field] is None:
                raise ValueError(f"QuantMS MSstats input is missing required column: {' or '.join(aliases)}")
        if self._columns["run"] is None and self._columns["reference"] is None:
            raise ValueError("QuantMS MSstats input requires a Run or Reference column")

        self._resolved = {
            "anchor_protein": self._columns["protein"],
            "pg_accessions": self._columns["protein"],
            "sequence": self._columns["peptide"],
            "peptidoform": self._columns["peptide"],
            "charge": self._columns["charge"],
            "intensities": self._columns["intensity"],
        }
        if self._columns["rt"]:
            self._resolved["rt"] = self._columns["rt"]
        if self._columns["reference"]:
            self._resolved["scan"] = self._columns["reference"]
        if self._columns["channel"]:
            self._resolved["intensity_label"] = self._columns["channel"]

    @staticmethod
    def _load_sdrf_context(path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        sdrf = pd.read_csv(path, sep="\t", dtype=str)
        sdrf.columns = sdrf.columns.str.strip().str.lower()
        required = {"source name", "comment[data file]", "comment[label]"}
        missing = sorted(required - set(sdrf.columns))
        if missing:
            raise ValueError(f"SDRF input is missing required column(s): {', '.join(missing)}")
        validate_sdrf_data_files(sdrf)

        alias_to_run: dict[str, str] = {}
        valid_pairs: list[dict[str, str]] = []
        seen_pairs: dict[tuple[str, str], str] = {}
        for row in sdrf.to_dict("records"):
            data_file = str(row["comment[data file]"]).strip()
            run_file_name = run_file_stem(data_file)
            label = normalize_label(str(row["comment[label]"]).strip())
            sample = str(row["source name"]).strip()
            if not run_file_name or not label or not sample:
                raise ValueError("SDRF source name, data file, and label values must be non-empty")

            pair = (run_file_name, label)
            previous = seen_pairs.get(pair)
            if previous is not None:
                raise ValueError(
                    f"SDRF run {run_file_name!r} maps canonical label {label!r} more than once: {previous!r} and {sample!r}"
                )
            seen_pairs[pair] = sample
            valid_pairs.append({"run_file_name": run_file_name, "label": label})

            for column in _RUN_ALIAS_COLUMNS:
                value = row.get(column)
                if value is None or pd.isna(value) or not str(value).strip():
                    continue
                key = _run_alias_key(value)
                if not key:
                    continue
                existing = alias_to_run.get(key)
                if existing is not None and existing != run_file_name:
                    raise ValueError(f"SDRF run alias {key!r} maps to both {existing!r} and {run_file_name!r}")
                alias_to_run[key] = run_file_name

        if not alias_to_run:
            raise ValueError("SDRF input contains no usable run aliases")
        aliases = pd.DataFrame(
            sorted(alias_to_run.items()),
            columns=["alias_key", "run_file_name"],
        )
        return aliases, pd.DataFrame(valid_pairs)

    def _build_channel_map(self, sdrf_path: str) -> pd.DataFrame | None:
        raw_labels = read_sdrf_labels(sdrf_path)
        experiment_type = experiment_type_from_labels(raw_labels)
        channel_column = self._columns["channel"]
        raw_channels = self._raw_channels(channel_column)

        if experiment_type == "LFQ":
            if raw_channels:
                raise ValueError("LFQ SDRF conflicts with non-empty MSstats Channel values")
            return None
        if channel_column is None:
            raise ValueError(f"{experiment_type} QuantMS MSstats input requires a Channel column")
        if not raw_channels:
            raise ValueError(f"{experiment_type} QuantMS MSstats input contains no Channel values")

        records = self._channel_records(experiment_type, raw_labels, raw_channels)
        return pd.DataFrame(records)

    def _raw_channels(self, channel_column: str | None) -> list[str]:
        """Return distinct non-empty MSstats channel values."""
        if channel_column is None:
            return []
        column = validate_identifier(channel_column)
        return [
            str(row[0]).strip()
            for row in self._conn.execute(
                sql_build(
                    "SELECT DISTINCT trim(CAST($channel AS VARCHAR)) FROM _msstats_input "
                    "WHERE NULLIF(trim(CAST($channel AS VARCHAR)), '') IS NOT NULL",
                    channel=column,
                )
            ).fetchall()
        ]

    @staticmethod
    def _channel_records(
        experiment_type: str,
        raw_labels: set[str] | None,
        raw_channels: list[str],
    ) -> list[dict]:
        """Resolve source channel values to canonical QPX labels."""
        numeric_positions = [position for value in raw_channels if (position := _as_channel_position(value))]
        position_map = resolve_channel_labels(experiment_type, raw_labels, numeric_positions)
        declared = {normalize_label(value) for value in raw_labels or set()}
        order_by_label = {label: order for order, label in position_map.items()}

        records: list[dict] = []
        for raw in raw_channels:
            position = _as_channel_position(raw)
            label = position_map.get(position) if position is not None else normalize_label(raw)
            if not label or label not in declared:
                raise ValueError(f"MSstats Channel {raw!r} does not match an SDRF {experiment_type} label")
            records.append(
                {
                    "source_channel": raw,
                    "canonical_label": label,
                    "channel_order": order_by_label.get(label, position or 0),
                }
            )
        return records

    def _register_peptides(self) -> None:
        peptide_column = validate_identifier(self._columns["peptide"])
        rows = self._conn.execute(
            sql_build(
                "SELECT DISTINCT trim(CAST($peptide AS VARCHAR)) FROM _msstats_input",
                peptide=peptide_column,
            )
        ).fetchall()
        records: list[dict[str, str]] = []
        for (raw_value,) in rows:
            raw = str(raw_value or "").strip()
            sequence, peptidoform, modifications = _parse_peptide(raw)
            self._modifications[raw] = modifications
            records.append(
                {
                    "peptide_raw": raw,
                    "sequence": sequence,
                    "peptidoform": peptidoform,
                }
            )
        self._register_frame("_peptide_lookup", pd.DataFrame(records))

    def _create_normalized_table(self, is_isobaric: bool) -> None:
        protein = validate_identifier(self._columns["protein"])
        peptide = validate_identifier(self._columns["peptide"])
        charge = validate_identifier(self._columns["charge"])
        intensity = validate_identifier(self._columns["intensity"])
        run = validate_identifier(self._columns["run"]) if self._columns["run"] else None
        reference = validate_identifier(self._columns["reference"]) if self._columns["reference"] else None
        channel = validate_identifier(self._columns["channel"]) if self._columns["channel"] else None
        rt = validate_identifier(self._columns["rt"]) if self._columns["rt"] else None

        run_expr = f"m.{run}" if run else "NULL::VARCHAR"
        reference_expr = f"m.{reference}" if reference else "NULL::VARCHAR"
        run_key = self._run_key_expression(run_expr)
        reference_key = self._run_key_expression(reference_expr)
        rt_expr = f"TRY_CAST(m.{rt} AS FLOAT)" if rt else "NULL::FLOAT"
        scan_expr = (
            f"TRY_CAST(regexp_extract(COALESCE(m.{reference}, ''), '(?i)scan=([0-9]+)', 1) AS INTEGER)"
            if reference
            else "NULL::INTEGER"
        )
        if is_isobaric:
            channel_join = f"LEFT JOIN _channel_map cm ON trim(CAST(m.{channel} AS VARCHAR)) = cm.source_channel"
            label_expr = "cm.canonical_label"
            order_expr = "cm.channel_order"
        else:
            channel_join = ""
            label_expr = "'LFQ'"
            order_expr = "1"

        query = sql_build(
            """
            CREATE OR REPLACE TEMP TABLE _normalized AS
            SELECT trim(CAST(m.$protein AS VARCHAR)) AS protein_raw,
                   trim(CAST(m.$peptide AS VARCHAR)) AS peptide_raw,
                   p.sequence,
                   p.peptidoform,
                   TRY_CAST(m.$charge AS SMALLINT) AS charge,
                   COALESCE(rr.run_file_name, ru.run_file_name) AS run_file_name,
                   rr.run_file_name AS reference_file,
                   ru.run_file_name AS run_file,
                   TRY_CAST(m.$intensity AS FLOAT) AS intensity,
                   $rt_expr AS rt,
                   $scan_expr AS scan_number,
                   $reference_expr AS reference_raw,
                   $label_expr AS intensity_label,
                   $order_expr AS channel_order,
                   count(DISTINCT trim(CAST(m.$protein AS VARCHAR)))
                     OVER (PARTITION BY p.sequence) = 1
                     AND NOT contains(trim(CAST(m.$protein AS VARCHAR)), ';') AS is_unique
            FROM _msstats_input m
            LEFT JOIN _run_aliases rr ON rr.alias_key = $reference_key
            LEFT JOIN _run_aliases ru ON ru.alias_key = $run_key
            LEFT JOIN _peptide_lookup p
              ON p.peptide_raw = trim(CAST(m.$peptide AS VARCHAR))
            $channel_join
            """,
            protein=protein,
            peptide=peptide,
            charge=charge,
            intensity=intensity,
            rt_expr=rt_expr,
            scan_expr=scan_expr,
            reference_expr=reference_expr,
            label_expr=label_expr,
            order_expr=order_expr,
            reference_key=reference_key,
            run_key=run_key,
            channel_join=channel_join,
        )
        self._conn.execute(query)

    @staticmethod
    def _run_key_expression(column: str) -> str:
        basename = "regexp_extract(replace(lower(coalesce(" + column + ", '')), '\\\\', '/'), '([^/]+)$', 1)"
        return (
            "CASE WHEN regexp_matches(" + basename + ", '\\.(mzml(?:\\.gz)?|mzxml|raw|d|wiff|mgf|dia)(?:$|[._\\s])') "
            "THEN regexp_extract(" + basename + ", '^(.*)(?:\\.(?:mzml(?:\\.gz)?|mzxml|raw|d|wiff|mgf|dia))(?:$|[._\\s])', 1) "
            "ELSE regexp_replace(" + basename + ", '\\.(?:mzml(?:\\.gz)?|mzxml|raw|d|wiff|mgf|dia)$', '') END"
        )

    def _validate_normalized_rows(self) -> None:
        total, invalid, unresolved, conflicts = self._conn.execute(
            """
            SELECT count(*),
                   count(*) FILTER (
                       WHERE NULLIF(protein_raw, '') IS NULL
                          OR NULLIF(peptidoform, '') IS NULL
                          OR NULLIF(sequence, '') IS NULL
                          OR charge IS NULL OR charge <= 0
                          OR intensity IS NULL OR NOT isfinite(intensity)
                          OR NULLIF(intensity_label, '') IS NULL
                   ),
                   count(*) FILTER (WHERE run_file_name IS NULL),
                   count(*) FILTER (
                       WHERE reference_file IS NOT NULL AND run_file IS NOT NULL
                         AND reference_file <> run_file
                   )
            FROM _normalized
            """
        ).fetchone()
        if invalid:
            raise ValueError(f"{invalid} of {total} QuantMS MSstats rows contain invalid required values")
        if unresolved:
            raise ValueError(f"{unresolved} of {total} QuantMS MSstats rows cannot be mapped to an SDRF data file")
        if conflicts:
            raise ValueError(f"{conflicts} QuantMS MSstats rows map Run and Reference to different SDRF files")

    def _validate_run_labels(self) -> None:
        count, run_file_name, label = self._conn.execute(
            """
            SELECT count(*), min(n.run_file_name), min(n.intensity_label)
            FROM _normalized n
            LEFT JOIN _valid_run_labels v
              ON v.run_file_name = n.run_file_name AND v.label = n.intensity_label
            WHERE v.run_file_name IS NULL
            """
        ).fetchone()
        if count:
            raise ValueError(f"{count} QuantMS MSstats rows have no SDRF sample for run {run_file_name!r} and label {label!r}")

    def _validate_feature_groups(self) -> None:
        conflict_count = self._conn.execute(
            """
            SELECT count(*) FROM (
                SELECT run_file_name, peptidoform, charge, rt, scan_number, intensity_label
                FROM _normalized
                GROUP BY ALL
                HAVING count(DISTINCT intensity) > 1
            )
            """
        ).fetchone()[0]
        if conflict_count:
            raise ValueError(f"{conflict_count} Feature/channel group(s) contain conflicting intensity values")

        location_count = self._conn.execute(
            """
            SELECT count(*) FROM (
                SELECT run_file_name, peptidoform, charge, rt, scan_number
                FROM _normalized
                GROUP BY ALL
                HAVING count(DISTINCT reference_raw) FILTER (WHERE reference_raw IS NOT NULL) > 1
            )
            """
        ).fetchone()[0]
        if location_count:
            raise ValueError(
                f"{location_count} Feature group(s) collapse distinct Reference values; "
                "the MSstats input lacks a sufficient RT/scan locator"
            )

    def _write_features(self, output_path: str, chunksize: int, creator: str) -> None:
        query = """
            WITH per_channel AS (
                SELECT run_file_name, peptide_raw, sequence, peptidoform, charge, rt, scan_number,
                       intensity_label, min(channel_order) AS channel_order,
                       first(intensity ORDER BY intensity) AS intensity
                FROM _normalized
                GROUP BY ALL
            ), proteins AS (
                SELECT run_file_name, peptidoform, charge, rt, scan_number,
                       list(DISTINCT protein_raw ORDER BY protein_raw) AS protein_names,
                       bool_and(is_unique) AS is_unique
                FROM _normalized
                GROUP BY ALL
            )
            SELECT c.run_file_name, c.peptide_raw, c.sequence, c.peptidoform, c.charge,
                   c.rt, c.scan_number,
                   list(c.intensity_label ORDER BY c.channel_order, c.intensity_label) AS labels,
                   list(c.intensity ORDER BY c.channel_order, c.intensity_label) AS intensity_values,
                   p.protein_names, p.is_unique
            FROM per_channel c
            JOIN proteins p
              ON p.run_file_name = c.run_file_name
             AND p.peptidoform = c.peptidoform
             AND p.charge = c.charge
             AND p.rt IS NOT DISTINCT FROM c.rt
             AND p.scan_number IS NOT DISTINCT FROM c.scan_number
            GROUP BY c.run_file_name, c.peptide_raw, c.sequence, c.peptidoform,
                     c.charge, c.rt, c.scan_number, p.protein_names, p.is_unique
            ORDER BY c.run_file_name, c.peptidoform, c.charge, c.rt, c.scan_number
        """
        written = 0
        with FeatureWriter(
            output_path,
            creator=creator,
            compression=self._compression,
            batch_size=chunksize,
            identity_composite=FEATURE_IDENTITY_COMPOSITE,
        ) as writer:
            for batch in self._query_batched(query, chunksize):
                records = [self._feature_record(row) for row in batch.to_pylist()]
                writer.write_batch(records)
                written += len(records)
        if not written:
            raise ValueError("QuantMS MSstats conversion produced no Feature rows")
        self.logger.info("Wrote %d Features to %s", written, output_path)

    def _feature_record(self, row: dict) -> dict:
        accessions = _split_proteins(row["protein_names"] or [])
        pg_accessions = [
            {"accession": accession, "start": None, "end": None, "pre": None, "post": None} for accession in accessions
        ]
        scan_number = row["scan_number"]
        return {
            "sequence": row["sequence"],
            "peptidoform": row["peptidoform"],
            "modifications": self._modifications[row["peptide_raw"]],
            "charge": row["charge"],
            "is_decoy": _is_decoy(accessions),
            "calculated_mz": None,
            "observed_mz": None,
            "run_file_name": row["run_file_name"],
            "scan": [scan_number] if scan_number is not None else [],
            "rt": row["rt"],
            "intensities": [
                {"label": label, "intensity": intensity}
                for label, intensity in zip(row["labels"], row["intensity_values"], strict=True)
            ],
            "pg_accessions": pg_accessions,
            "anchor_protein": accessions[0] if accessions else None,
            "unique": row["is_unique"],
            "id_run_file_name": row["run_file_name"],
        }

    def _register_frame(self, name: str, frame: pd.DataFrame) -> None:
        table = validate_table(name)
        self._conn.execute(sql_build("DROP TABLE IF EXISTS $table", table=table))
        self._conn.from_df(frame).create(table)
