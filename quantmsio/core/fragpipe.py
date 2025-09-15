import logging
import uuid
from datetime import date

from typing import Iterator, List, Optional, Any, Tuple
from dataclasses import dataclass, field, asdict
from pathlib import Path
from array import array

import pandas as pd

import pyarrow as pa
import pyarrow.parquet as pq

from pyteomics import proforma

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


@dataclass
class Parameter:
    name: str
    controlled_vocabulary: Optional[str] = None
    curie: Optional[str] = None
    value: Optional[Any] = None

    @classmethod
    def parse(cls, line: str) -> "Parameter":
        cv, curie, name, value = line.split(",", 3)
        return cls(
            cv.strip() if cv else cv,
            curie.strip() if curie else curie,
            name.strip(),
            value.strip(),
        )

    def format(self) -> str:
        return f"[{self.controlled_vocabulary or ''}, {self.curie or ''}, {self.name}, {self.value or ''}]"


@dataclass
class ModificationPositions:
    positions: List[int] = field(default_factory=list)

    def append(self, i):
        self.positions.append(i)

    def extend(self, it):
        self.positions.extend(it)

    def __nonzero__(self):
        return bool(self.positions)

    def __getitem__(self, i):
        return self.positions[i]

    def __len__(self) -> int:
        return len(self.positions)

    def __iter__(self):
        yield from self.positions

    def format(self, parameters=None) -> str:
        if not parameters:
            return "|".join(map(str, self.positions))
        else:
            if len(parameters) != len(self):
                raise ValueError(
                    f"{self.__class__.__name__}.format requires the parameter list"
                    f" length ({len(parameters)}) match its length ({len(self)})"
                )
            return "|".join(
                f"{i}{p.format()}" for i, p in zip(self.positions, parameters)
            )


@dataclass
class ModificationIdentifier:
    curie: str
    mass_value: Optional[float] = field(init=False)

    def __post_init__(self):
        if isinstance(self.curie, float):
            self.mass_value = self.curie
            self.curie = f"CHEMMOD:{self.mass_value:0.4f}"
        if ":" not in self.curie:
            raise ValueError(
                f"Modification identifier must either be a mass or a CURIE-formatted string, got {self.curie}"
            )

    def is_chemmod(self) -> bool:
        return self.curie.startswith("CHEMMOD")

    def __hash__(self):
        return hash(self.curie)

    def __str__(self):
        return self.curie


@dataclass
class MzTabModification:
    identifier: ModificationIdentifier
    position: Optional[ModificationPositions] = None
    parameter: Optional[List[Parameter]] = None
    neutral_loss: Optional[Parameter] = None

    def is_chemmod(self) -> bool:
        return self.identifier.is_chemmod()

    def format(self) -> str:
        if self.position is None:
            position = "null"
        else:
            position = self.position.format(self.parameter)

        if not self.neutral_loss:
            return f"{position}-{self.identifier}"
        else:
            raise NotImplementedError("Neutral losses are not yet implemented")

    def __str__(self) -> str:
        return self.format()

    @classmethod
    def arrow_type(cls):
        localization_struct = pa.struct(
            [
                pa.field("position", pa.int32(), nullable=False),
                pa.field("localization_probability", pa.float32(), nullable=False),
            ]
        )
        return pa.struct(
            [
                pa.field("name", pa.string(), nullable=False),
                pa.field("fields", pa.list_(localization_struct), nullable=False),
            ]
        )

    def to_arrow_single(self):
        fields = []
        if self.position:
            for pos in self.position:
                fields.append({"position": pos, "localization_probability": 1.0})
        return {"name": str(self.identifier), "fields": fields}

    @classmethod
    def to_arrow(cls, batch: Iterator[List["MzTabModification"]]):
        modifications = []
        for block in batch:
            modifications.append([p.to_arrow_single() for p in block])
        return pa.array(modifications, type=cls.arrow_type())


@dataclass
class AssignedModification:
    position: int
    amino_acid: str
    mass: float

    @classmethod
    def parse(cls, text: str) -> List["AssignedModification"]:
        items = []
        if not text or not isinstance(text, str):
            return items
        for token in text.split(","):
            pos_aa, mass_ = token.split("(")
            pos_aa = pos_aa.strip()
            mass = float(mass_[:-1])
            amino_acid = pos_aa[-1]
            if pos_aa == "N-term":
                position = 0
            else:
                position = int(pos_aa[:-1])
            items.append(cls(position, amino_acid, mass))
        return items

    def as_mztab(self) -> MzTabModification:
        position = ModificationPositions([self.position])
        mod_ident = ModificationIdentifier(self.mass)
        return MzTabModification(mod_ident, position=position)


@dataclass
class Peptide:
    sequence: str
    modifications: List[MzTabModification] = field(default_factory=list)

    protein_positions: List[str] = field(default_factory=list)
    protein_accessions: List[str] = field(default_factory=list)

    gene_names: List[str] = field(default_factory=list)
    gene_accessions: List[str] = field(default_factory=list)

    @property
    def peptidoform(self) -> str:
        peptide = proforma.ProForma.parse(self.sequence)

        mod: MzTabModification
        for mod in self.modifications:
            if len(mod.position) > 1:
                raise NotImplementedError("Variably localized modifications")
            pos = mod.position[0]
            mod_tag = None
            if mod.is_chemmod():
                mod_tag = proforma.MassModification(mod.identifier.mass_value)
            else:
                raise NotImplementedError("CV-based modifications")
            if pos == 0:
                # n-term
                if not peptide.n_term:
                    peptide.n_term = []
                peptide.n_term.append(mod_tag)
            else:
                index = pos - 1
                (aa, mods_at_index) = peptide.sequence[index]
                if not mods_at_index:
                    mods_at_index = []
                mods_at_index.append(mod_tag)
                peptide.sequence[index] = (aa, mods_at_index)
        return peptide

    @classmethod
    def to_arrow(cls, batch: Iterator["Peptide"]):
        sequences = []
        peptidoforms = []
        modifications = []

        for peptide in batch:
            sequences.append(peptide.sequence)
            peptidoforms.append(str(peptide.peptidoform))
            modifications.append(
                [mod.to_arrow_single() for mod in peptide.modifications]
            )

        return {
            "sequence": pa.array(sequences, type=pa.string()),
            "peptidoform": pa.array(peptidoforms, type=pa.string()),
            "modifications": pa.array(
                modifications, type=pa.list_(MzTabModification.arrow_type())
            ),
        }

    @classmethod
    def from_series(cls, row: pd.Series):
        return peptide_from_row(row)


@dataclass
class Spectrum:
    observed_mz: float
    precursor_charge: int
    retention_time: float
    scan_number: str
    source_file: str
    ion_mobility: Optional[float] = None

    @classmethod
    def to_arrow(cls, batch: Iterator["Spectrum"]):
        observed_mz = array("f")
        precursor_charge = array("i")
        retention_time = array("f")
        scan_number = []
        source_file = []
        ion_mobility = []

        for s in batch:
            observed_mz.append(s.observed_mz)
            precursor_charge.append(s.precursor_charge)
            retention_time.append(s.retention_time)
            scan_number.append(s.scan_number)
            source_file.append(s.source_file)
            ion_mobility.append(s.ion_mobility)
        mz_arrays = pa.array(
            [None] * len(precursor_charge), type=pa.list_(pa.float32())
        )
        intensity_arrays = pa.array(
            [None] * len(precursor_charge), type=pa.list_(pa.float32())
        )
        num_peaks_array = pa.array([None] * len(precursor_charge), type=pa.int32())
        charge_arrays = pa.array(
            [None] * len(precursor_charge), type=pa.list_(pa.int32())
        )
        ion_type_arrays = pa.array(
            [None] * len(precursor_charge), type=pa.list_(pa.string())
        )
        ion_mobility_arrays = pa.array(
            [None] * len(precursor_charge), type=pa.list_(pa.float32())
        )
        return {
            "precursor_charge": pa.array(precursor_charge, type=pa.int32()),
            "observed_mz": pa.array(observed_mz, type=pa.float32()),
            "rt": pa.array(retention_time, type=pa.float32()),
            "scan": pa.array(scan_number, type=pa.string()),
            "source_file": pa.array(source_file, type=pa.string(), from_pandas=True),
            "ion_mobility": pa.array(ion_mobility, type=pa.float32(), from_pandas=True),
            "num_peaks": num_peaks_array,
            "mz_array": mz_arrays,
            "intensity_array": intensity_arrays,
            "charge_array": charge_arrays,
            "ion_type_array": ion_type_arrays,
            "ion_mobility_array": ion_mobility_arrays,
        }

    @classmethod
    def from_series(cls, row: pd.Series):
        return spectrum_from_row(row)


@dataclass
class Score:
    name: str
    score: float

    @classmethod
    def arrow_type(cls):
        return pa.struct(
            [
                pa.field("name", pa.string(), nullable=False),
                pa.field("score", pa.float32(), nullable=False),
            ]
        )


@dataclass
class Identification:
    spectrum: Spectrum
    peptide: Peptide
    calculated_mz: float

    cv_params: List[Parameter] = field(default_factory=list)
    additional_scores: List[Score] = field(default_factory=list)
    global_qvalue: Optional[float] = None
    posterior_error_probability: Optional[float] = None
    rank: Optional[int] = None
    predicted_retention_time: Optional[float] = None

    @classmethod
    def to_arrow(cls, batch: List["Identification"]) -> pa.Table:
        data_arrays = {
            "global_qvalue": pa.array(
                [i.global_qvalue for i in batch], type=pa.float64()
            ),
            "posterior_error_probability": pa.array(
                [i.posterior_error_probability for i in batch], type=pa.float32()
            ),
            "calculated_mz": pa.array(
                [i.calculated_mz for i in batch], type=pa.float32()
            ),
            "additional_scores": pa.array(
                [[asdict(s) for s in i.additional_scores] for i in batch],
                type=pa.list_(Score.arrow_type()),
            ),
            "rank": pa.array([i.rank for i in batch], type=pa.int32()),
            "predicted_rt": pa.array(
                [i.predicted_retention_time for i in batch],
                type=pa.float32(),
                from_pandas=True,
            ),
        }

        data_arrays.update(Peptide.to_arrow(i.peptide for i in batch))

        data_arrays.update(Spectrum.to_arrow(i.spectrum for i in batch))
        return pa.table(data_arrays)

    @classmethod
    def from_series(cls, row: pd.Series):
        return identification_from_row(row)

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame):
        return df.apply(cls.from_series, 1).tolist()


@dataclass
class FragPipe:
    output_directory: Path

    def write_psms_to_parquet(
        self,
        file_path: Path,
        batch_size: int = 10000,
        output_prefix_file: Optional[str] = None,
        **metadata,
    ):
        if not file_path.exists():
            raise FileNotFoundError(file_path)
        if not self.output_directory.exists():
            self.output_directory.mkdir(parents=True)
        if not output_prefix_file:
            output_prefix_file = "psm"

        file_uuid = uuid.uuid4()
        output_path = (
            self.output_directory / f"{output_prefix_file}-{file_uuid}.psm.parquet"
        )

        metadata["file_type"] = "psm"
        metadata["uuid"] = str(file_uuid)
        metadata["creation_date"] = date.today().isoformat()

        logger.debug(
            "Writing FragPipe PSMs to %s with a batch size of %d",
            output_path,
            batch_size,
        )
        writer = None

        file_metadata = []
        try:
            for i, batch in enumerate(
                self.convert_psms(file_path, batch_size=batch_size)
            ):
                logger.debug("Converting batch %d with %d entries", i, batch.num_rows)
                if writer is None:
                    logger.debug(
                        "Initializing ParquetWriter with schema %r", batch.schema
                    )
                    writer = pq.ParquetWriter(
                        output_path,
                        schema=batch.schema,
                        metadata_collector=file_metadata,
                    )
                    writer.add_key_value_metadata(metadata)

                writer.write_batch(batch)
        finally:
            if writer is not None:
                writer.close()
            else:
                logger.warning("No PSMs found. Not writing PSM parquet file")
        return file_metadata

    @staticmethod
    def convert_psms(
        file_path: Path,
        batch_size: int = 10000,
    ) -> Iterator[pa.RecordBatch]:
        iterator = pd.read_csv(file_path, iterator=True, chunksize=batch_size, sep="\t")
        for batch in iterator:
            idents = Identification.from_dataframe(batch)
            table = Identification.to_arrow(idents)
            yield from table.to_batches(batch_size)


def parse_spectrum_id(identifier: str) -> Tuple[str, str]:
    tokens = identifier.split(".")
    source_file_name = tokens[0]
    scan_number = tokens[1].lstrip("0")
    return source_file_name, scan_number


def spectrum_from_row(row: pd.Series) -> Spectrum:
    observed_mz = row["Calibrated Observed M/Z"]
    precursor_charge = row["Charge"]
    retention_time = row["Retention"]
    source_file, scan_number = parse_spectrum_id(row["Spectrum"])
    return Spectrum(
        observed_mz,
        precursor_charge,
        retention_time,
        scan_number,
        source_file,
    )


def peptide_from_row(row: pd.Series) -> Peptide:
    protein_accession = row["Protein ID"]
    gene_name = row["Gene"]
    start = row["Protein Start"]
    end = row["Protein End"]
    return Peptide(
        row["Peptide"],
        [
            mod.as_mztab()
            for mod in AssignedModification.parse(row["Assigned Modifications"])
        ],
        [f"{start}:{end}"],
        [protein_accession],
        [gene_name],
        [],
    )


def identification_from_row(row: pd.Series) -> Identification:
    spectrum = spectrum_from_row(row)
    peptide = peptide_from_row(row)

    calculated_mz = row["Calculated M/Z"]

    posterior = row["PeptideProphet Probability"]

    cv_params = [
        Parameter(
            "number of missed cleavages",
            "PSI-MS",
            str(row["Number of Missed Cleavages"]),
            "MS:1003044",
        ),
        Parameter(
            "number of enzymatic termini",
            "PSI-MS",
            str(row["Number of Enzymatic Termini"]),
            "MS:1003048",
        ),
    ]

    additional_scores = [
        Score("MSFragger:Hyperscore", row["Hyperscore"]),
        Score("MSFragger:Expectation", row["Expectation"]),
    ]

    return Identification(
        spectrum,
        peptide,
        calculated_mz,
        cv_params=cv_params,
        additional_scores=additional_scores,
        posterior_error_probability=posterior,
    )
