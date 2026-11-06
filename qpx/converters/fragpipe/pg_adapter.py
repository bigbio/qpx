"""FragPipe PG adapter -- combined_protein.tsv to pg.parquet.

Loads a FragPipe ``combined_protein.tsv`` into DuckDB, transforms
protein-group rows into ``PgSchema``, and writes through ``PgWriter``.

Key schema changes:
    - ``is_decoy`` (int, 0/1) -> ``is_decoy`` (bool)
    - Intensities struct: ``{sample_accession, channel, intensity}`` ->
      ``{label, intensity}``
"""

from __future__ import annotations

import logging

import pandas as pd

from qpx.converters.base import BaseConverter, resolve_columns
from qpx.converters.channel_labels import experiment_runs_from_sdrf, experiment_type_from_labels, read_sdrf_labels
from qpx.converters.fragpipe.constants import is_decoy_accession, protein_accessions
from qpx.converters.fragpipe.isobaric import TMT_REPORT_INTENSITY, channel_intensities, require_channel_columns
from qpx.converters.mappings import get_field_mappings
from qpx.converters.utils import is_contaminant_accession, safe_float
from qpx.core.files import run_file_stem
from qpx.core.sql import escape_path, sql_build
from qpx.writers.pg import PgWriter

logger = logging.getLogger(__name__)

# Derive field map from central YAML mappings
_PG_MAP = get_field_mappings("fragpipe", "pg")

# Additional quantities of an experiment by column suffix, named as in pg.md.
_ADDITIONAL_INTENSITIES = (
    ("MaxLFQ Intensity", "maxlfq"),
    ("MaxLFQ Unique Intensity", "maxlfq_unique"),
    ("Unique Intensity", "unique_intensity"),
)
_SPECTRAL_COUNTS = (
    ("Spectral Count", "spectral_count"),
    ("Unique Spectral Count", "unique_spectral_count"),
    ("Total Spectral Count", "total_spectral_count"),
)
# Per-experiment column suffixes, longest first so that "<exp> MaxLFQ Intensity"
# is not read as an experiment "<exp> MaxLFQ".
_EXPERIMENT_SUFFIXES = sorted(
    {suffix for suffix, _ in _ADDITIONAL_INTENSITIES + _SPECTRAL_COUNTS}
    | {"Intensity", "Total Intensity", "Razor Intensity", "Razor Spectral Count", "Total Peptides", "Unique Peptides"},
    key=len,
    reverse=True,
)


def _experiment_quantities(row, experiment: str) -> tuple[float | None, list[dict]] | None:
    """Primary intensity and additional quantities of one experiment; ``None`` without evidence in it.

    ``<exp> Intensity`` is the primary quantity (``Total Intensity`` in older
    reports, else MaxLFQ). A zero intensity is FragPipe's "not quantified" and is
    left out; spectral counts are kept as counted, zeros included. A group with
    spectral counts but no intensity in the experiment keeps a row with a null
    intensity.
    """

    def value(suffix: str) -> float | None:
        return safe_float(row.get(f"{experiment} {suffix}"))

    def positive(suffix: str) -> float | None:
        quantity = value(suffix)
        return quantity if quantity is not None and quantity > 0 else None

    primary = next(filter(None, map(positive, ("Intensity", "Total Intensity", "MaxLFQ Intensity"))), None)
    quantities = [(name, positive(suffix)) for suffix, name in _ADDITIONAL_INTENSITIES]
    counts = [(name, value(suffix)) for suffix, name in _SPECTRAL_COUNTS]
    if primary is None and not any(count for _, count in counts):
        return None
    return primary, [
        {"intensity_name": name, "intensity_value": quantity} for name, quantity in quantities + counts if quantity is not None
    ]


def normalize_run_name(value: str) -> str:
    """Return a QPX run name for a FragPipe raw-file value."""
    name = run_file_stem(str(value))
    return name[: -len(".htrms")] if name.casefold().endswith(".htrms") else name


def normalize_experiment_mapping(mapping: dict[str, list[str]]) -> dict[str, list[str]]:
    """Normalize and validate a FragPipe experiment-to-runs mapping."""
    normalized: dict[str, list[str]] = {}
    for experiment, runs in mapping.items():
        experiment_name = str(experiment).strip()
        if not experiment_name:
            raise ValueError("FragPipe experiment name cannot be empty")
        if isinstance(runs, (str, bytes)):
            raise ValueError(f"Runs for FragPipe experiment {experiment_name!r} must be a list")
        run_names = sorted({normalize_run_name(str(run)) for run in runs if str(run).strip()})
        if not run_names:
            raise ValueError(f"FragPipe experiment {experiment_name!r} has no raw files")
        normalized[experiment_name] = run_names
    return normalized


def load_experiment_annotation(path: str) -> dict[str, list[str]]:
    """Read FragPipe ``experiment_annotation.tsv`` into grouped runs."""
    annotation = pd.read_csv(path, sep="\t", dtype=str)
    columns = {str(column).strip().casefold(): column for column in annotation.columns}
    missing = {"file", "sample"} - columns.keys()
    if missing:
        raise ValueError("FragPipe experiment annotation is missing column(s): " + ", ".join(sorted(missing)))

    mapping: dict[str, list[str]] = {}
    file_col = columns["file"]
    sample_col = columns["sample"]
    for row_number, row in annotation.iterrows():
        experiment = str(row[sample_col]).strip()
        file_value = str(row[file_col]).strip()
        if not experiment or not file_value or experiment == "nan" or file_value == "nan":
            raise ValueError(f"FragPipe experiment annotation has an empty file/sample value on row {row_number + 2}")
        mapping.setdefault(experiment, []).append(file_value)
    if not mapping:
        raise ValueError("FragPipe experiment annotation contains no file/sample rows")
    return normalize_experiment_mapping(mapping)


def require_label_free(sdrf_path: str | None) -> None:
    """Refuse an isobaric (TMT/iTRAQ) SDRF.

    ``combined_ion.tsv`` and ``combined_protein.tsv`` hold MS1 intensities per
    experiment; the reporter-ion quantities are in the channel columns of the
    one-experiment reports, which are not converted. Writing the MS1 intensities
    as the label-free quantity of an isobaric experiment would be wrong.
    """
    experiment_type = experiment_type_from_labels(read_sdrf_labels(sdrf_path))
    if experiment_type in ("TMT", "iTRAQ"):
        raise ValueError(
            f"FragPipe {experiment_type} reporter-ion quantities are not converted: combined_ion.tsv and "
            "combined_protein.tsv hold MS1 intensities, which would be written as label-free. Convert each plex's ion.tsv and protein.tsv with --annotation-file, or psm.tsv only."
        )


def resolve_experiment_runs(
    experiments: list[str],
    experiment_to_runs: dict[str, list[str]] | None,
    experiment_annotation_path: str | None,
    sdrf_path: str | None,
) -> dict[str, list[str]]:
    """Map each FragPipe experiment to its raw runs, raising when one has none.

    The results mapping (``experiment_annotation.tsv`` or an explicit
    ``experiment_to_runs``) comes first, then the SDRF ``source name`` of the same
    name. An experiment in neither cannot be tied to run.parquet, so conversion
    fails with the full list rather than writing the experiment as a run.
    """
    if experiment_to_runs is not None and experiment_annotation_path is not None:
        raise ValueError("Provide either experiment_to_runs or experiment_annotation_path, not both")
    if experiment_annotation_path is not None:
        mapping = load_experiment_annotation(experiment_annotation_path)
    else:
        mapping = normalize_experiment_mapping(experiment_to_runs or {})
    sdrf_mapping = experiment_runs_from_sdrf(sdrf_path) or {}
    runs = {experiment: mapping.get(experiment) or sdrf_mapping.get(experiment) for experiment in experiments}
    missing = sorted(experiment for experiment, members in runs.items() if not members)
    if missing:
        raise ValueError(
            "No raw runs for FragPipe experiment(s) "
            + ", ".join(missing)
            + ": no experiment->runs mapping. Provide experiment_annotation.tsv "
            "(experiment_annotation_path), an explicit experiment_to_runs mapping, "
            "or an SDRF (sdrf_path) whose 'source name' matches the experiment(s)."
        )
    return runs


class FragPipePgAdapter(BaseConverter):
    """Convert FragPipe ``combined_protein.tsv`` to ``pg.parquet``.

    Usage::

        with FragPipePgAdapter() as adapter:
            adapter.convert(
                protein_path="combined_protein.tsv",
                output_path="pg.parquet",
            )
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._experiment_runs: dict[str, list[str]] = {}
        self._channels: dict[str, str] = {}
        self._tmt_report: dict[str, dict[str, float]] = {}

    def load_channels(self, channels: dict[str, str], tmt_report: dict[str, dict[str, float]] | None = None) -> None:
        """Read reporter-ion quantities from these columns of a one-plex protein.tsv.

        *tmt_report* holds the TMT-Integrator log2 abundance of each protein and label.
        """
        self._channels = channels
        self._tmt_report = tmt_report or {}

    def convert(
        self,
        protein_path: str,
        output_path: str,
        chunksize: int = 100_000,
        creator: str = "fragpipe",
        experiment_to_runs: dict[str, list[str]] | None = None,
        experiment_annotation_path: str | None = None,
        sdrf_path: str | None = None,
    ) -> None:
        """Run the combined_protein.tsv -> pg.parquet conversion.

        Args:
            protein_path: Path to FragPipe ``combined_protein.tsv``.
            output_path: Destination Parquet path.
            chunksize: Records per batch; a row gives one record per
                experiment, so ``chunksize // experiments`` rows are read at a time.
            creator: Creator tag in Parquet metadata.
            experiment_to_runs: Optional mapping of each FragPipe *experiment* (the
                token that names the per-experiment intensity columns in
                ``combined_protein.tsv``) to its member ``run_file_name`` values.
                A FragPipe experiment can aggregate several raw files (fractions),
                whereas run.parquet is keyed by raw file, so ``grouped_runs`` must
                be expanded to those files to survive the sample join.
            experiment_annotation_path: Optional FragPipe
                ``experiment_annotation.tsv``. Its ``sample`` column names the
                experiment and ``file`` lists the member raw files.
            sdrf_path: Optional SDRF, used as the *fallback* source for
                ``grouped_runs`` when a FragPipe experiment has no
                ``experiment_annotation.tsv`` / ``experiment_to_runs`` mapping: the
                experiment is matched to a ``source name`` and grouped_runs becomes
                that sample's raw files. If neither the results mapping nor the
                SDRF resolves an experiment, conversion raises rather than emitting
                a dangling experiment token.
        """
        if not self._channels:
            require_label_free(sdrf_path)

        # Step 1: Load protein file into DuckDB
        self._load_protein_file(protein_path)

        # Step 2: Resolve column mappings against actual input columns
        actual_cols = {
            c[0]
            for c in self._conn.execute(
                "SELECT column_name FROM information_schema.columns WHERE table_name='fragpipe_proteins'"
            ).fetchall()
        }
        self._resolved = resolve_columns(_PG_MAP, actual_cols)

        # Step 3: Detect per-experiment intensity columns and resolve each one to
        # its raw files before writing anything.
        if self._channels:
            require_channel_columns(self._channels, actual_cols, protein_path)
        experiment_cols = list(experiment_to_runs) if self._channels else self._detect_experiment_columns()
        self._experiment_runs = resolve_experiment_runs(
            experiment_cols, experiment_to_runs, experiment_annotation_path, sdrf_path
        )

        # Step 4: Stream and transform
        self.logger.info("Transforming FragPipe protein groups ...")

        with PgWriter(output_path, creator=creator) as writer:
            for batch in self._query_batched("SELECT * FROM fragpipe_proteins", max(1, chunksize // len(experiment_cols))):
                df = batch.to_pandas()
                records = self._transform_batch(df, experiment_cols)
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)

        self.logger.info(f"FragPipe PG conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_protein_file(self, path: str) -> None:
        """Load combined_protein.tsv into DuckDB."""
        self._conn.execute(
            sql_build(
                """CREATE TABLE fragpipe_proteins AS
            SELECT * FROM read_csv_auto('$path',
                delim='\t', header=true, auto_detect=true)""",
                path=escape_path(path),
            )
        )
        count = self._conn.execute("SELECT COUNT(*) FROM fragpipe_proteins").fetchall()[0][0]
        self.logger.info(f"Loaded {count:,} FragPipe protein group rows")

    def _detect_experiment_columns(self) -> list[str]:
        """Detect the experiments from any per-experiment column.

        FragPipe ``combined_protein.tsv`` has columns like ``<experiment> Intensity``,
        ``<experiment> MaxLFQ Intensity`` and ``<experiment> Spectral Count``; a
        spectral-count-only report has just the counts. ``Combined *`` columns
        cover all experiments.
        """
        cols = self._conn.execute(
            "SELECT column_name FROM information_schema.columns WHERE table_name='fragpipe_proteins'"
        ).fetchall()
        experiments = set()
        for (col,) in cols:
            # Unprefixed columns of a one-experiment protein.tsv name no experiment.
            if col in _EXPERIMENT_SUFFIXES or col.startswith("Combined "):
                continue
            suffix = next((suffix for suffix in _EXPERIMENT_SUFFIXES if col.endswith(f" {suffix}")), None)
            if suffix:
                experiments.add(col[: -len(suffix) - 1])
        if not experiments:
            raise ValueError(
                "FragPipe combined_protein.tsv has no per-experiment columns such as '<experiment> Intensity'; "
                "a one-plex protein.tsv holds its reporter quantities in channel columns, which --annotation-file names"
            )
        return sorted(experiments)

    # ------------------------------------------------------------------
    # Transform
    # ------------------------------------------------------------------

    def _transform_batch(self, df: pd.DataFrame, experiments: list[str]) -> list[dict]:
        return self._transform_rows_with_skip(
            df.to_dict("records"), lambda row: self._transform_row(row, experiments), "FragPipe PG"
        )

    def _transform_row(self, row, experiments: list[str]) -> list[dict]:
        """Transform a single combined_protein.tsv row.

        Expands into one record per experiment (run).
        """
        records: list[dict] = []
        common_fields = self._pg_common_fields(row)
        if self._channels:
            return self._channel_records(row, common_fields, experiments[0])

        # Expand per experiment
        for experiment in experiments:
            quantities = _experiment_quantities(row, experiment)
            if quantities is None:
                continue
            primary_intensity, extra_vals = quantities

            # Intensities (new schema: {label, intensity})
            label = "LFQ"
            intensities = [{"label": label, "intensity": primary_intensity}]

            # Additional intensities: MaxLFQ, unique intensities and spectral counts
            additional_intensities = [{"label": label, "intensities": extra_vals}] if extra_vals else []

            # Per-experiment peptide counts; FragPipe reports no feature counts
            unique_peptides = safe_float(row.get(f"{experiment} Unique Peptides"))
            total_peptides = safe_float(row.get(f"{experiment} Total Peptides"))

            rec = {
                **common_fields,
                # Expand the FragPipe experiment to its member run files (results
                # mapping, else SDRF); the tool-reported per-experiment intensity
                # is kept as-is (FragPipe already aggregated across the fractions).
                "grouped_runs": self._experiment_runs[experiment],
                "intensities": intensities,
                "additional_intensities": additional_intensities or None,
                "peptide_counts": (
                    None
                    if unique_peptides is None or total_peptides is None
                    else {"unique_sequences": int(unique_peptides), "total_sequences": int(total_peptides)}
                ),
                "feature_counts": None,
            }
            records.append(rec)

        return records

    def _channel_records(self, row, common_fields: dict, plex: str) -> list[dict]:
        """One record per quantified channel of a one-plex protein.tsv."""
        unique_peptides, total_peptides = safe_float(row.get("Unique Peptides")), safe_float(row.get("Total Peptides"))
        peptide_counts = (
            None
            if unique_peptides is None or total_peptides is None
            else {"unique_sequences": int(unique_peptides), "total_sequences": int(total_peptides)}
        )
        report = self._tmt_report.get(row.get("Protein ID"), {})
        records = []
        for entry in channel_intensities(row, self._channels) or []:
            if entry["intensity"] is None:
                continue
            abundance = report.get(entry["label"])
            records.append(
                {
                    **common_fields,
                    "grouped_runs": self._experiment_runs[plex],
                    "intensities": [entry],
                    "additional_intensities": (
                        None
                        if abundance is None
                        else [
                            {
                                "label": entry["label"],
                                "intensities": [{"intensity_name": TMT_REPORT_INTENSITY, "intensity_value": abundance}],
                            }
                        ]
                    ),
                    "peptide_counts": peptide_counts,
                    "feature_counts": None,
                }
            )
        return records

    def _pg_qvalue(self, row: dict) -> float | None:
        """Resolve FDR or convert the source protein probability."""
        r = self._resolved
        # Protein group q-value. Real FDR/Q-Value columns are preferred (see the
        # ordered mapping in column_mappings.yaml). "Protein Probability" is a
        # confidence score (~1.0 = best), so when it is the only column present it
        # is converted to an FDR-like value as 1 - probability, never used raw.
        pg_qvalue_col = r.get("pg_qvalue")
        pg_qvalue_raw = safe_float(row.get(pg_qvalue_col)) if pg_qvalue_col else None
        if pg_qvalue_raw is None:
            pg_qvalue = None
        elif pg_qvalue_col == "Protein Probability":
            pg_qvalue = 1.0 - pg_qvalue_raw
        else:
            pg_qvalue = pg_qvalue_raw
        return pg_qvalue

    def _pg_common_fields(self, row: dict) -> dict:
        """Build protein-group fields shared across experiments."""
        r = self._resolved
        # Protein group identity
        protein_raw = row.get(r.get("pg_accessions", "Protein"))
        pg_accessions = protein_accessions(protein_raw, row.get("Indistinguishable Proteins"))

        gene_raw = row.get(r.get("gg_accessions", "Gene"), "")
        gg_accessions = str(gene_raw).split(",") if pd.notna(gene_raw) and gene_raw else None

        # Protein description/name
        description = row.get(r.get("pg_names", "Description"), "")
        pg_names = None
        if pd.notna(description) and description:
            pg_names = [str(description)]

        anchor_protein = pg_accessions[0] if pg_accessions else ""

        # Is decoy: detected from a decoy prefix on any RAW accession
        is_decoy = is_decoy_accession(protein_raw)

        # Combined count, the same for every member of an indistinguishable group
        total_peptides = int(row.get(r.get("peptide_count_total", "Combined Total Peptides"), 0) or 0)

        # Sequence coverage
        seq_coverage = safe_float(row.get(r.get("sequence_coverage", "Percent Coverage")))

        # Molecular weight
        mol_weight = safe_float(row.get(r.get("molecular_weight", "Protein Molecular Weight (Da)")))
        if mol_weight is not None:
            mol_weight = mol_weight / 1000.0  # Convert Da to kDa

        pg_qvalue = self._pg_qvalue(row)

        # Peptides per protein
        peptides = [{"protein_name": acc, "peptide_count": total_peptides} for acc in pg_accessions]
        return {
            "pg_accessions": pg_accessions,
            "pg_names": pg_names,
            "gg_accessions": gg_accessions,
            "gg_names": gg_accessions,
            "gg_qvalue": None,
            "anchor_protein": anchor_protein,
            "global_qvalue": None,
            "pg_qvalue": pg_qvalue,
            "is_decoy": is_decoy,
            "contaminant": is_contaminant_accession(protein_raw),
            "peptides": peptides,
            "sequence_coverage": seq_coverage,
            "molecular_weight": mol_weight,
            **self._pg_annotations(row),
        }

    @staticmethod
    def _pg_annotations(row: dict) -> dict:
        """Protein length and existence as cv_params, and the ProteinProphet probabilities as scores."""
        cv_params = []
        length = safe_float(row.get("Protein Length", row.get("Length")))
        if length is not None:
            cv_params.append({"cv_name": "sequence_length", "cv_value": str(int(length))})
        existence = row.get("Protein Existence")
        if pd.notna(existence) and str(existence).strip():
            # "1:Experimental evidence at protein level" -> "1"
            cv_params.append({"cv_name": "protein_existence", "cv_value": str(existence).split(":", 1)[0].strip()})
        scores = [
            {"score_name": name, "score_value": score, "higher_better": True}
            for column, name in (
                ("Protein Probability", "protein_probability"),
                ("Top Peptide Probability", "top_peptide_probability"),
            )
            if (score := safe_float(row.get(column))) is not None
        ]
        return {"cv_params": cv_params or None, "additional_scores": scores or None}
