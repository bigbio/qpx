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
from qpx.converters.channel_labels import experiment_runs_from_sdrf
from qpx.converters.fragpipe.constants import is_decoy_accession
from qpx.converters.mappings import get_field_mappings
from qpx.converters.utils import safe_float
from qpx.core.sql import escape_path, sql_build
from qpx.writers.pg import PgWriter

logger = logging.getLogger(__name__)

# Derive field map from central YAML mappings
_PG_MAP = get_field_mappings("fragpipe", "pg")
_RUN_EXTENSIONS = (".mzml", ".raw", ".d", ".wiff", ".htrms")


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
        self._experiment_to_runs: dict[str, list[str]] = {}
        self._sdrf_experiment_to_runs: dict[str, list[str]] | None = None

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
            chunksize: Rows per batch.
            creator: Creator tag in Parquet metadata.
            experiment_to_runs: Optional mapping of each FragPipe *experiment* (the
                token that names the per-experiment intensity columns in
                ``combined_protein.tsv``) to its member ``run_file_name`` values.
                A FragPipe experiment can aggregate several raw files (fractions),
                whereas run.parquet is keyed by raw file, so ``grouped_runs`` must
                be expanded to those files to survive the sample join. When the
                mapping is absent, ``grouped_runs`` falls back to the single
                experiment token, which is correct only when experiment equals the
                raw-file name (no fractions).
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
        if experiment_to_runs is not None and experiment_annotation_path is not None:
            raise ValueError("Provide either experiment_to_runs or experiment_annotation_path, not both")
        if experiment_annotation_path is not None:
            self._experiment_to_runs = self._load_experiment_annotation(experiment_annotation_path)
        else:
            self._experiment_to_runs = self._normalize_experiment_mapping(experiment_to_runs or {})
        self._sdrf_experiment_to_runs = experiment_runs_from_sdrf(sdrf_path)

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

        # Step 3: Detect per-experiment intensity columns and require every one to
        # resolve to raw files (results mapping first, then SDRF) — fail early with
        # the full list rather than emitting a dangling experiment token.
        experiment_cols = self._detect_experiment_columns()
        missing = sorted(exp for exp in experiment_cols if self._lookup_grouped_runs(exp) is None)
        if missing:
            raise ValueError(
                "Cannot build grouped_runs for FragPipe experiment(s) "
                + ", ".join(missing)
                + ": no experiment->runs mapping. Provide experiment_annotation.tsv "
                "(experiment_annotation_path), an explicit experiment_to_runs mapping, "
                "or an SDRF (sdrf_path) whose 'source name' matches the experiment(s)."
            )

        # Step 4: Stream and transform
        self.logger.info("Transforming FragPipe protein groups ...")

        with PgWriter(output_path, creator=creator, compression=self._compression) as writer:
            for batch in self._query_batched("SELECT * FROM fragpipe_proteins", chunksize):
                df = batch.to_pandas()
                records = self._transform_batch(df, experiment_cols)
                if records:
                    self._track_scores(records)
                    writer.write_batch(records)

        self.logger.info(f"FragPipe PG conversion complete -> {output_path}")

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _lookup_grouped_runs(self, experiment: str) -> list[str] | None:
        """Resolve an experiment to its raw files, or None if unmapped.

        Results mapping (experiment_annotation.tsv / experiment_to_runs) first,
        then the SDRF ``source name`` fallback. No dangling experiment token.
        """
        runs = self._experiment_to_runs.get(experiment)
        if runs:
            return runs
        if self._sdrf_experiment_to_runs:
            runs = self._sdrf_experiment_to_runs.get(experiment)
            if runs:
                return runs
        return None

    def _resolve_grouped_runs(self, experiment: str) -> list[str]:
        """Like :meth:`_lookup_grouped_runs` but raises when unmapped."""
        runs = self._lookup_grouped_runs(experiment)
        if runs is None:
            raise ValueError(
                f"Cannot build grouped_runs for FragPipe experiment {experiment!r}: no "
                "experiment->runs mapping. Provide experiment_annotation.tsv, an explicit "
                "experiment_to_runs mapping, or an SDRF (sdrf_path) whose 'source name' "
                "matches the experiment."
            )
        return runs

    @staticmethod
    def _normalize_run_name(value: str) -> str:
        """Return a QPX run_file_name from a FragPipe annotation file value."""
        name = str(value).strip().replace("\\", "/").rsplit("/", 1)[-1]
        lower_name = name.casefold()
        for extension in _RUN_EXTENSIONS:
            if lower_name.endswith(extension):
                return name[: -len(extension)]
        return name

    @classmethod
    def _normalize_experiment_mapping(
        cls,
        mapping: dict[str, list[str]],
    ) -> dict[str, list[str]]:
        """Normalize and validate an explicit experiment-to-runs mapping."""
        normalized: dict[str, list[str]] = {}
        for experiment, runs in mapping.items():
            experiment_name = str(experiment).strip()
            if not experiment_name:
                raise ValueError("FragPipe experiment name cannot be empty")
            if isinstance(runs, (str, bytes)):
                raise ValueError(f"Runs for FragPipe experiment {experiment_name!r} must be a list")
            run_names = sorted({cls._normalize_run_name(run) for run in runs if str(run).strip()})
            if not run_names:
                raise ValueError(f"FragPipe experiment {experiment_name!r} has no raw files")
            normalized[experiment_name] = run_names
        return normalized

    @classmethod
    def _load_experiment_annotation(
        cls,
        path: str,
    ) -> dict[str, list[str]]:
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
        return cls._normalize_experiment_mapping(mapping)

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
        count = self._conn.execute("SELECT COUNT(*) FROM fragpipe_proteins").fetchone()[0]
        self.logger.info(f"Loaded {count:,} FragPipe protein group rows")

    def _detect_experiment_columns(self) -> list[str]:
        """Detect per-experiment intensity columns.

        FragPipe ``combined_protein.tsv`` has columns like:
            ``<experiment> Total Intensity``
            ``<experiment> Spectral Count``
            ``<experiment> Unique Spectral Count``
        """
        cols = self._conn.execute(
            "SELECT column_name FROM information_schema.columns WHERE table_name='fragpipe_proteins'"
        ).fetchall()
        col_names = [c[0] for c in cols]

        # Find columns ending with " Total Intensity" or " MaxLFQ Intensity"
        intensity_cols = [c for c in col_names if c.endswith(" Total Intensity") or c.endswith(" MaxLFQ Intensity")]

        # Extract experiment names
        experiments = set()
        for col in intensity_cols:
            for suffix in [" Total Intensity", " MaxLFQ Intensity"]:
                if col.endswith(suffix):
                    experiments.add(col[: -len(suffix)])

        return sorted(experiments)

    # ------------------------------------------------------------------
    # Transform
    # ------------------------------------------------------------------

    def _transform_batch(self, df: pd.DataFrame, experiments: list[str]) -> list[dict]:
        records: list[dict] = []
        skipped = 0
        for row in df.to_dict("records"):
            try:
                recs = self._transform_row(row, experiments)
                records.extend(recs)
            except Exception as e:
                skipped += 1
                self.logger.debug(f"Skipping FragPipe PG row: {e}")
        if skipped:
            total = skipped + len(records)
            self.logger.warning(
                "Skipped %d / %d rows (%.1f%%) in batch",
                skipped,
                total,
                100 * skipped / total if total else 0,
            )
        return records

    def _transform_row(self, row, experiments: list[str]) -> list[dict]:
        """Transform a single combined_protein.tsv row.

        Expands into one record per experiment (run).
        """
        records: list[dict] = []
        r = self._resolved  # shorthand for resolved column mappings

        # Protein group identity
        protein_raw = str(row.get(r.get("pg_accessions", "Protein"), ""))
        pg_accessions = protein_raw.split(",") if protein_raw else []

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

        # Combined counts
        total_peptides = int(row.get(r.get("peptide_count_total", "Combined Total Peptides"), 0) or 0)
        unique_peptides = int(row.get(r.get("peptide_count_unique", "Combined Unique Peptides"), 0) or 0)

        # Sequence coverage
        seq_coverage = safe_float(row.get(r.get("sequence_coverage", "Percent Coverage")))

        # Molecular weight
        mol_weight = safe_float(row.get(r.get("molecular_weight", "Protein Molecular Weight (Da)")))
        if mol_weight is not None:
            mol_weight = mol_weight / 1000.0  # Convert Da to kDa

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

        # Peptides per protein
        peptides = [{"protein_name": acc, "peptide_count": total_peptides} for acc in pg_accessions]

        # Expand per experiment
        for experiment in experiments:
            total_int_col = f"{experiment} Total Intensity"
            maxlfq_col = f"{experiment} MaxLFQ Intensity"
            spec_count_col = f"{experiment} Spectral Count"
            unique_spec_col = f"{experiment} Unique Spectral Count"

            total_intensity = safe_float(row.get(total_int_col)) or 0.0
            maxlfq_val = safe_float(row.get(maxlfq_col))
            # MaxLFQ-only experiments (detected from a "<exp> MaxLFQ Intensity"
            # column) have no Total Intensity; fall back to MaxLFQ as the primary
            # intensity rather than dropping the protein group entirely.
            if total_intensity > 0:
                primary_intensity = total_intensity
            elif maxlfq_val is not None and maxlfq_val > 0:
                primary_intensity = maxlfq_val
            else:
                continue

            # Intensities (new schema: {label, intensity})
            label = "LFQ"
            intensities = [{"label": label, "intensity": float(primary_intensity)}]

            # Additional intensities pre-computed by FragPipe (MaxLFQ)
            additional_intensities = []
            extra_vals = []
            if maxlfq_val is not None:
                extra_vals.append({"intensity_name": "maxlfq", "intensity_value": float(maxlfq_val)})
            if extra_vals:
                additional_intensities.append({"label": label, "intensities": extra_vals})

            # Additional scores
            additional_scores = []
            spec_val = safe_float(row.get(spec_count_col))
            if spec_val is not None:
                additional_scores.append(
                    {
                        "score_name": "spectral_count",
                        "score_value": spec_val,
                        "higher_better": True,
                    }
                )

            # Per-experiment peptide/feature counts
            exp_unique = safe_float(row.get(unique_spec_col)) or 0
            exp_total = safe_float(row.get(spec_count_col)) or 0

            rec = {
                "pg_accessions": pg_accessions,
                "pg_names": pg_names,
                "gg_accessions": gg_accessions,
                "gg_names": gg_accessions,  # Gene symbols serve as both accession and name
                "gg_qvalue": None,
                "anchor_protein": anchor_protein,
                # Expand the FragPipe experiment to its member run files (results
                # mapping, else SDRF); the tool-reported per-experiment intensity
                # is kept as-is (FragPipe already aggregated across the fractions).
                "grouped_runs": self._resolve_grouped_runs(experiment),
                "global_qvalue": None,
                "pg_qvalue": pg_qvalue,
                "intensities": intensities,
                "additional_intensities": additional_intensities or None,
                "is_decoy": is_decoy,
                "contaminant": None,
                "peptides": peptides,
                "peptide_counts": {
                    "unique_sequences": unique_peptides,
                    "total_sequences": total_peptides,
                },
                "feature_counts": {
                    "unique_features": int(exp_unique),
                    "total_features": int(exp_total),
                },
                "sequence_coverage": seq_coverage,
                "molecular_weight": mol_weight,
                "additional_scores": additional_scores or None,
                "cv_params": None,
            }
            records.append(rec)

        return records
