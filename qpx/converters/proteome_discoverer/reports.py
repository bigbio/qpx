"""What the Proteome Discoverer exports share: modifications, accessions, files and the PSM table."""

from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path, PureWindowsPath

import pandas as pd

from qpx.converters.base import BaseConverter, resolve_columns
from qpx.converters.mappings import get_field_mappings
from qpx.converters.ptm import build_proforma, from_proforma, modification_tag, position_site
from qpx.core.files import run_file_stem
from qpx.core.sql import escape_path, sql_build

PSM_MAP = get_field_mappings("proteome_discoverer", "psm")

_MOD_TOKEN = re.compile(r"^(?:(N-Term)|(C-Term)|[A-Z](\d+))(?:\(Prot\))?\((.+)\)$")
_LIST_SEPARATOR = re.compile(r"\s*;\s*")


def split_accessions(text: object) -> list[str]:
    """The accessions of a ``;``-separated cell, empty for a missing one."""
    if not isinstance(text, str):
        return []
    return [accession for accession in _LIST_SEPARATOR.split(text.strip()) if accession]


@lru_cache(maxsize=100_000)
def peptide_identity(sequence: str, modifications: str) -> tuple[str, list[dict] | None]:
    """Peptidoform and modifications from PD's ``N-Term(TMTpro); T4(Phospho)``."""
    mods = []
    for token in split_accessions(modifications):
        match = _MOD_TOKEN.match(token)
        if match is None:
            raise ValueError(f"Unreadable Proteome Discoverer modification {token!r} on {sequence}")
        n_term, c_term, residue_position, name = match.groups()
        position = 0 if n_term else len(sequence) + 1 if c_term else int(residue_position)
        residue, terminus = position_site(sequence, position)
        mods.append((position, modification_tag(name, residue, terminus)))
    return from_proforma(build_proforma(sequence, mods), sequence, meta=None)


def runs_by_file_id(input_files_path: str | Path) -> dict[str, str]:
    """File ID (``F16.1``) -> run from an ``InputFiles`` export, whose names are Windows paths."""
    files = pd.read_csv(input_files_path, sep="\t", dtype=str)
    id_column = "Study File ID" if "Study File ID" in files.columns else "File ID"
    return {file_id: run_file_stem(PureWindowsPath(name).name) for file_id, name in zip(files[id_column], files["File Name"])}


class PdPsmTable(BaseConverter):
    """The PSM export in DuckDB, with the columns the views need under one spelling."""

    _resolved: dict[str, str] = {}

    def load_psms(self, path: str, runs: dict[str, str] | None) -> None:
        """Load ``PSMs.txt``; a run comes from ``Spectrum File``, else from *runs* by ``File ID``."""
        self._conn.execute(
            sql_build(
                """CREATE TABLE _pd_raw AS SELECT * FROM read_csv('$path', delim='\t', header=true,
                    quote='"', all_varchar=true, max_line_size=20000000, parallel=false)""",
                path=escape_path(path),
            )
        )
        columns = self.get_table_columns("_pd_raw")
        self._resolved = resolve_columns(PSM_MAP, columns)
        self._conn.execute(sql_build("CREATE TABLE pd_psms AS SELECT *, $derived FROM _pd_raw", derived=self._derived(runs)))
        self._conn.execute("DROP TABLE _pd_raw")

    def _column(self, field: str) -> str:
        column = self._resolved.get(field)
        if column is None:
            raise ValueError(f"The Proteome Discoverer PSM export has no column for {field}: {PSM_MAP[field]}")
        return '"' + column.replace('"', '""') + '"'

    def _optional_double(self, field: str) -> str:
        column = self._resolved.get(field)
        return f'TRY_CAST("{column}" AS DOUBLE)' if column else "CAST(NULL AS DOUBLE)"

    def _derived(self, runs: dict[str, str] | None) -> str:
        return ", ".join(
            (
                f"{self._run_expression(runs)} AS _run",
                # "[R].hSDtDSDRR.[E]" and "hSDtDSDRR" spell the same peptide.
                f"upper(regexp_replace(regexp_replace({self._column('sequence')}, '^\\[.*?\\]\\.', ''),"
                " '\\.\\[.*?\\]$', '')) AS _sequence",
                f"COALESCE({self._column('modifications')}, '') AS _modifications",
                f"CAST({self._column('charge')} AS INTEGER) AS _charge",
                f"CAST({self._column('scan')} AS INTEGER) AS _scan",
                f"COALESCE({self._optional_double('rank')}, 1) AS _rank",
                f"{self._optional_double('qvalue')} AS _qvalue",
            )
        )

    def _run_expression(self, runs: dict[str, str] | None) -> str:
        if "run_file_name" in self._resolved:
            return (
                f"regexp_replace(regexp_replace({self._column('run_file_name')}, '^.*[\\\\/]', ''), '\\.(raw|mzML|mzml|d)$', '')"
            )
        if not runs:
            raise ValueError("This PSM export has no Spectrum File column; pass its InputFiles export for the raw runs")
        cases = " ".join(f"WHEN {_literal(file_id)} THEN {_literal(run)}" for file_id, run in runs.items())
        return f"CASE {self._column('file_id')} {cases} END"

    def runs_by_file_id(self) -> dict[str, str]:
        """File ID -> run as the PSM export itself names them, for the file groups of the pg view."""
        rows = self._conn.execute(f"SELECT DISTINCT {self._column('file_id')}, _run FROM pd_psms").fetchall()
        return dict(rows)

    def spectrum_fields(self, row: dict) -> dict:
        """The PSM's precursor, retention time and search values."""
        r = self._resolved
        return {
            "run_file_name": row["_run"],
            "scan": [int(row["_scan"])],
            "observed_mz": _number(row, r, "observed_mz"),
            "mass_error_ppm": _number(row, r, "mass_error_ppm"),
            "rt": _seconds(_number(row, r, "rt")),
            "posterior_error_probability": _number(row, r, "posterior_error_probability"),
            "missed_cleavages": _integer(_number(row, r, "missed_cleavages")),
            "ion_mobility": _number(row, r, "ion_mobility"),
        }

    def search_scores(self, row: dict) -> list[dict] | None:
        """The search engine scores, the rank, and the precursor's share of the isolation window."""
        scores = [
            {"score_name": name, "score_value": value, "higher_better": higher_better}
            for column, name, higher_better in _SCORES
            if (value := _float(row.get(column))) is not None
        ]
        interference = _number(row, self._resolved, "isolation_interference")
        if interference is not None:
            scores.append({"score_name": "parent_ion_fraction", "score_value": 1 - interference / 100, "higher_better": True})
        return scores or None


_SCORES = (
    ("Rank", "rank", False),
    ("XCorr", "sequest_xcorr", True),
    ("Comet XCorr", "comet_xcorr", True),
    ("Ions Score", "mascot_ions_score", True),
    ("Delta Score", "delta_score", True),
    ("DeltaScore", "delta_score", True),
    ("SVM Score", "percolator_svm_score", True),
    ("Percolator SVMScore", "percolator_svm_score", True),
)


def _literal(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def _float(value: object) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _number(row: dict, resolved: dict, field: str) -> float | None:
    return _float(row.get(resolved.get(field)))


def _integer(value: float | None) -> int | None:
    return None if value is None else int(value)


def _seconds(minutes: float | None) -> float | None:
    return None if minutes is None else minutes * 60
