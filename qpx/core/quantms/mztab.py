# Standard library imports
import gzip
import os
import shutil
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, Generator, List, Optional, Union
import logging

# Third-party imports
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import re

# Local application imports
from qpx.core.common import OPENMS_PROTEIN_QVALUE_WORDS
from qpx.core.duckdb import DuckDB
from qpx.core.sdrf import SDRFHandler
from qpx.utils.constants import (
    ACCESSION,
    AMBIGUITY_MEMBERS,
    ANCHOR_PROTEIN,
    DECOY_PREFIXES,
    INDISTINGUISHABLE_GROUP,
    ITRAQ_CHANNEL,
    MSSTATS_INTENSITY,
    MSSTATS_PEPTIDE_SEQUENCE,
    MSSTATS_PROTEIN_NAME,
    MSSTATS_REFERENCE,
    MSSTATS_REFERENCE_NAME,
    MSSTATS_NAME_MAP,
    MSSTATS_USECOLS,
    MZTAB_PROTEIN_BEST_SEARCH_ENGINE_SCORE,
    OPT_GLOBAL_RESULT_TYPE,
    MZTAB_PROTEIN_COLUMNS,
    PROTEIN_DETAILS_MZTAB,
    SINGLE_PROTEIN_MZTAB,
    SPECTRA_REF,
    SPECTRA_REF_FILE,
    SPECTRA_REF_SCAN,
    TMT_CHANNELS,
)
from qpx.utils.file_utils import (
    validate_extension,
    validate_file,
)
from qpx.utils.pride_utils import (
    generate_scan_number,
    standardize_protein_string_accession,
)
from qpx.utils.mztab_utils import (
    parse_metadata_line,
    parse_header_line,
    parse_protein_line,
    parse_psm_line,
    is_mztab_line_type,
)
from qpx.utils.system import check_disk_space


@contextmanager
def log_time(logger, step: str) -> Generator[None, None, None]:
    """Context manager to log time taken for a step.

    Args:
        logger: Logger instance to use for logging
        step: Description of the step being timed

    Yields:
        None: Context manager yields nothing

    Example:
        with log_time(logger, "processing proteins"):
            process_proteins()
    """
    start = time.time()
    try:
        yield
    finally:
        duration = time.time() - start
        logger.debug(f"Time taken for {step}: {duration:.2f} seconds")


@contextmanager
def parquet_writer_context(
    writer: Optional[pq.ParquetWriter],
) -> Generator[Optional[pq.ParquetWriter], None, None]:
    """Context manager for parquet writers to ensure proper cleanup.

    Args:
        writer: ParquetWriter instance or None

    Yields:
        ParquetWriter instance

    Example:
        with parquet_writer_context(writer) as w:
            w.write_table(table)
    """
    try:
        yield writer
    finally:
        if writer:
            try:
                writer.close()
            except Exception as e:
                logging.getLogger(__name__).warning(
                    f"Error closing parquet writer: {e}"
                )


class MzTabIndexer(DuckDB):
    """
    DuckDB implementation specific to mzTab data with advanced MSstats analysis capabilities.

    This class provides high-performance analysis of mzTab and MSstats data using either
    DuckDB or Parquet backends. It includes sophisticated analysis methods for proteomics
    data that leverage SQL for optimal performance.

    Key Features:
    - Fast indexing and querying of mzTab files
    - Advanced MSstats analysis methods
    - Support for both DuckDB and Parquet backends
    - SDRF integration for sample mapping
    - Comprehensive protein and peptide statistics
    - FDR calculations and quality metrics
    """

    _MZTAB_INDEXER_FOLDER = "mztab_indexer_"

    _MZTAB_INDEXER_TABLE_METADATA = "metadata"
    _MZTAB_INDEXER_TABLE_PROTEINS = "proteins"
    _MZTAB_INDEXER_TABLE_PROTEIN_DETAILS = "protein_details"
    _MZTAB_INDEXER_TABLE_PSMS = "psms"
    _MZTAB_INDEXER_TABLE_MSSTATS = "msstats"

    @classmethod
    def create(
        cls,
        mztab_path: Union[Path, str],
        msstats_path: Optional[Union[Path, str]] = None,
        sdrf_path: Optional[Union[Path, str]] = None,
        database_path: Optional[Union[Path, str]] = None,
        **kwargs,
    ) -> "MzTabIndexer":
        """
        Create a new MzTabIndexer database or store.

        This factory method creates a new MzTabIndexer instance and processes the mzTab file
        to create the database or parquet store. Use this when you have a mzTab file to process.

        Args:
            mztab_path: Path to mzTab file (required for creation)
            msstats_path: Optional path to MSstats file for additional analysis
            sdrf_path: Optional path to SDRF file for sample mapping
            database_path: Path for database file (.duckdb) or directory
            **kwargs: Additional arguments passed to __init__

        Returns:
            New MzTabIndexer instance with created database/store

        Raises:
            ValueError: If required parameters are missing or invalid
            FileNotFoundError: If mzTab file doesn't exist

        Example:
            >>> indexer = MzTabIndexer.create(
            ...     mztab_path="data.mzTab",
            ...     database_path="output.duckdb",
            ... )
        """
        database_path = str(database_path)
        return cls(
            mztab_path=mztab_path,
            msstats_path=msstats_path,
            database_path=database_path,
            sdrf_path=sdrf_path,
            **kwargs,
        )

    @classmethod
    def open(
        cls,
        database_path: Union[Path, str],
        sdrf_path: Optional[Union[Path, str]] = None,
        **kwargs,
    ) -> "MzTabIndexer":
        """
        Open an existing MzTabIndexer database or store.

        This factory method opens an existing MzTabIndexer database or parquet store.
        Use this when you have already created a database and want to access it.

        Args:
            database_path: Path to existing database file (.duckdb) or directory (for parquet backend)
            sdrf_path: Optional path to SDRF file for sample mapping
            **kwargs: Additional arguments passed to __init__

        Returns:
            MzTabIndexer instance connected to existing database/store

        Raises:
            ValueError: If database doesn't exist or is invalid
            FileNotFoundError: If database file/directory doesn't exist

        Example:
            >>> indexer = MzTabIndexer.open(
            ...     database_path="existing.duckdb",
            ... )
        """
        database_path = str(database_path)
        return cls(
            mztab_path=None,
            msstats_path=None,
            database_path=database_path,
            sdrf_path=sdrf_path,
            **kwargs,
        )

    def __init__(
        self,
        mztab_path: Optional[Union[Path, str]] = None,
        msstats_path: Optional[Union[Path, str]] = None,
        database_path: Optional[Union[Path, str]] = None,
        max_memory: str = "16GB",
        worker_threads: int = 8,
        batch_size: int = 100000,
        protein_columns: Optional[List[str]] = MZTAB_PROTEIN_COLUMNS,
        psm_columns: Optional[List[str]] = None,
        sdrf_path: Optional[Union[Path, str]] = None,
    ):
        """Initialize MzTabIndexer.

        Args:
            mztab_path: Path to mzTab file. Required if database_path is not provided.
            msstats_path: Optional path to MSstats file.
            backend: The storage and query backend to use ('duckdb' or 'parquet').
            database_path: Path to the DuckDB database file (required for 'duckdb' backend).
            max_memory: Maximum memory to use.
            worker_threads: Number of worker threads.
            batch_size: Number of rows to process in each batch.
            protein_columns: Optional list of protein columns to include (if None, include all).
            psm_columns: Optional list of PSM columns to include (if None, include all).
            sdrf_path: Optional path to SDRF file for sample mapping.
        """
        if not database_path:
            raise ValueError("`database_path` is required for backend.")

        self._database_path = str(database_path)
        self._mztab_path = mztab_path
        self._msstats_path = msstats_path

        self._batch_size = batch_size
        self._worker_threads = worker_threads
        self._max_memory = max_memory
        self._protein_columns = protein_columns
        self._psm_columns = psm_columns

        if sdrf_path:
            if not validate_file(sdrf_path):
                raise ValueError("A valid sdrf_path must be provided.")

            self._sdrf_handler = SDRFHandler(sdrf_path)
            self.experiment_type = self._sdrf_handler.get_experiment_type_from_sdrf()
            self._sample_map = self._sdrf_handler.get_sample_map_run()

        # Initialize temporary directory as None
        self._temp_parquet_dir = None

        if self._mztab_path and not validate_file(self._mztab_path):
            raise ValueError("A valid mztab_path must be provided.")

        if self._msstats_path and not validate_file(self._msstats_path):
            raise ValueError("A valid msstats_path must be provided.")

        # Case 1: Opening existing database/store (no mzTab file provided)
        if not self._mztab_path:
            # Initialize parent class first to get logger
            super().__init__(
                database_path=self._database_path,
                max_memory=self._max_memory,
                worker_threads=self._worker_threads,
            )

            self.logger.info(
                f"Opening existing database in DuckDB mode: {self._database_path}"
            )
            if not validate_extension(self._database_path, ".duckdb"):
                raise ValueError(
                    f"Database file {self._database_path} is not a DuckDB database file."
                )
            if not validate_file(self._database_path):
                raise ValueError(
                    f"Database file {self._database_path} does not exist or is empty."
                )
            if not self.check_tables_in_db(
                [
                    self._MZTAB_INDEXER_TABLE_METADATA,
                    self._MZTAB_INDEXER_TABLE_PROTEINS,
                    self._MZTAB_INDEXER_TABLE_PROTEIN_DETAILS,
                    self._MZTAB_INDEXER_TABLE_PSMS,
                ],
            ):
                raise ValueError(
                    f"Database file {self._database_path} does not contain the expected mztab sections tables: metadata, proteins, protein_details, psms"
                )

        # Case 2: Creating new database/store (mzTab file provided)
        else:
            # For DuckDB backend, we need to check if the file should have .duckdb extension
            if not str(self._database_path).endswith(".duckdb"):
                raise ValueError(
                    f"Database file {self._database_path} is not a DuckDB database file, use .duckdb extension."
                )
            # Check if file already exists
            if Path(self._database_path).exists():
                raise ValueError(
                    f"Database file {self._database_path} already exists, use a different file name."
                )

            # Initialize DuckDB parent class first to get logger
            super().__init__(
                database_path=self._database_path,
                max_memory=self._max_memory,
                worker_threads=self._worker_threads,
            )

            self.logger.info(
                f"Creating new database in DuckDB mode: {self._database_path}"
            )
            # Create tables from mzTab file
            self.create_database_tables()

        # Initialize metadata caching attributes
        self._metadata_cache = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        """Release resources, if any."""
        if hasattr(self, "_sdrf_handler") and self._sdrf_handler:
            self._sdrf_handler = None
        # Close DuckDB connection
        super().destroy_database()

    def _get_output_dir(self) -> Path:
        """
        Get the output directory for parquet files based on backend type.

        This method determines the appropriate directory for storing parquet files
        based on the current backend configuration.

        Returns:
            Path to the output directory:
            - For DuckDB backend: Temporary directory in parent of database file
            - For Parquet backend: The database_path directory itself

        Note:
            For DuckDB backend, parquet files are temporary and are cleaned up
            after loading into the database. For Parquet backend, the files
            are the final storage format.
        """
        # For DuckDB backend, use for parquet temp files the parent directory of the database file
        if self._temp_parquet_dir is None:
            self._temp_parquet_dir = (
                Path(self._database_path).parent / self._MZTAB_INDEXER_FOLDER
            )
            self._temp_parquet_dir.mkdir(parents=True, exist_ok=True)
        return self._temp_parquet_dir

    def _load_parquet_tables_to_duckdb(self, parquet_dir: Path):
        """
        Map parquet tables into the in-memory DuckDB instance.

        This method loads existing parquet files into the DuckDB database as tables.
        It's used when opening an existing parquet store in DuckDB backend mode.

        Args:
            parquet_dir: Path to directory containing parquet files

        Raises:
            ValueError: If required parquet files are missing from the directory

        Note:
            This method expects the following parquet files to exist:
            - metadata.parquet
            - proteins.parquet
            - protein_details.parquet
            - psms.parquet
            - msstats.parquet (optional)
        """
        parquet_dir = Path(parquet_dir)  # Ensure it's a Path object
        self.logger.info(
            f"Mapping parquet tables into the in memory DuckDB from {self._database_path}"
        )
        # Create tables from parquet files
        if (
            parquet_dir.exists()
            and Path(parquet_dir / "metadata.parquet").exists()
            and Path(parquet_dir / "proteins.parquet").exists()
            and Path(parquet_dir / "protein_details.parquet").exists()
            and Path(parquet_dir / "psms.parquet").exists()
        ):
            self._duckdb.execute(
                f"CREATE TABLE {self._MZTAB_INDEXER_TABLE_METADATA} AS SELECT * FROM '{parquet_dir / 'metadata.parquet'}'"
            )
            self._duckdb.execute(
                f"CREATE TABLE {self._MZTAB_INDEXER_TABLE_PROTEINS} AS SELECT * FROM '{parquet_dir / 'proteins.parquet'}'"
            )
            self._duckdb.execute(
                f"CREATE TABLE {self._MZTAB_INDEXER_TABLE_PROTEIN_DETAILS} AS SELECT * FROM '{parquet_dir / 'protein_details.parquet'}'"
            )
            self._duckdb.execute(
                f"CREATE TABLE {self._MZTAB_INDEXER_TABLE_PSMS} AS SELECT * FROM '{parquet_dir / 'psms.parquet'}'"
            )
        else:
            raise ValueError(
                f"Parquet directory {parquet_dir} does not contain the expected mztab sections files: metadata.parquet, proteins.parquet, protein_details.parquet, psms.parquet"
            )

        # Add MSstats table if available
        if Path(parquet_dir / "msstats.parquet").exists():
            self._duckdb.execute(
                f"CREATE TABLE {self._MZTAB_INDEXER_TABLE_MSSTATS} AS SELECT * FROM '{parquet_dir / 'msstats.parquet'}'"
            )

    # def _build_parquet_database(self):
    #     """
    #     Build parquet database from mzTab file.

    #     This method processes the mzTab file to create parquet files and then
    #     loads them into the DuckDB database. It handles both the initial processing
    #     and the database table creation.

    #     The method processes the following sections:
    #     - Metadata (MTD lines)
    #     - Proteins (PRT lines)
    #     - Protein details (PRT lines with specific result types)
    #     - PSMs (PSM lines)
    #     - MSstats (if provided)

    #     Raises:
    #         ValueError: If mzTab file path is not provided
    #         RuntimeError: If processing fails

    #     Note:
    #         For DuckDB backend, temporary parquet files are created and then
    #         cleaned up after loading into the database.
    #     """
    #     if not self._mztab_path:
    #         raise ValueError(
    #             "A valid mztab_path must be provided for creating a new database/store."
    #         )

    #     self.logger.info(
    #         f"Building parquet database from mzTab file: {self._mztab_path}"
    #     )

    #     # Process mzTab to parquet files
    #     (
    #         metadata_parquet,
    #         proteins_parquet,
    #         protein_details_parquet,
    #         psms_parquet,
    #     ) = self._process_mztab_to_parquet()

    #     if self._msstats_path:
    #         msstats_parquet = self._process_msstats_to_parquet()

    #     # Create DuckDB tables from parquet files
    #     self.logger.debug("Creating DuckDB tables from parquet files")

    #     if self._MZTAB_INDEXER_TABLE_METADATA not in [
    #         t[0] for t in self._duckdb.execute("SHOW TABLES").fetchall()
    #     ]:
    #         self.logger.info(f"Creating table: {self._MZTAB_INDEXER_TABLE_METADATA}")
    #         self._duckdb.execute(
    #             f"CREATE TABLE {self._MZTAB_INDEXER_TABLE_METADATA} AS SELECT * FROM read_parquet('{metadata_parquet}')"
    #         )

    #     if self._MZTAB_INDEXER_TABLE_PROTEINS not in [
    #         t[0] for t in self._duckdb.execute("SHOW TABLES").fetchall()
    #     ]:
    #         if (
    #             os.path.exists(proteins_parquet)
    #             and os.path.getsize(proteins_parquet) > 0
    #         ):
    #             self.logger.info(
    #                 f"Creating table: {self._MZTAB_INDEXER_TABLE_PROTEINS}"
    #             )
    #             self._duckdb.execute(
    #                 f"CREATE TABLE {self._MZTAB_INDEXER_TABLE_PROTEINS} AS SELECT * FROM read_parquet('{proteins_parquet}')"
    #             )

    #     if self._MZTAB_INDEXER_TABLE_PROTEIN_DETAILS not in [
    #         t[0] for t in self._duckdb.execute("SHOW TABLES").fetchall()
    #     ]:
    #         if (
    #             os.path.exists(protein_details_parquet)
    #             and os.path.getsize(protein_details_parquet) > 0
    #         ):
    #             self.logger.info(
    #                 f"Creating table: {self._MZTAB_INDEXER_TABLE_PROTEIN_DETAILS}"
    #             )
    #             self._duckdb.execute(
    #                 f"CREATE TABLE {self._MZTAB_INDEXER_TABLE_PROTEIN_DETAILS} AS SELECT * FROM read_parquet('{protein_details_parquet}')"
    #             )

    #     if self._MZTAB_INDEXER_TABLE_PSMS not in [
    #         t[0] for t in self._duckdb.execute("SHOW TABLES").fetchall()
    #     ]:
    #         if os.path.exists(psms_parquet) and os.path.getsize(psms_parquet) > 0:
    #             self.logger.info(f"Creating table: {self._MZTAB_INDEXER_TABLE_PSMS}")
    #             self._duckdb.execute(
    #                 f"CREATE TABLE {self._MZTAB_INDEXER_TABLE_PSMS} AS SELECT * FROM read_parquet('{psms_parquet}')"
    #             )

    #     # Handle MSstats if available
    #     if self._msstats_path:
    #         if self._MZTAB_INDEXER_TABLE_MSSTATS not in [
    #             t[0] for t in self._duckdb.execute("SHOW TABLES").fetchall()
    #         ]:
    #             self.logger.info(f"Creating table: {self._MZTAB_INDEXER_TABLE_MSSTATS}")
    #             self._duckdb.execute(
    #                 f"CREATE TABLE {self._MZTAB_INDEXER_TABLE_MSSTATS} AS SELECT * FROM read_parquet('{msstats_parquet}')"
    #             )

    #     # Clean up temporary parquet files for DuckDB backend
    #     if self._backend == "duckdb" and self._temp_parquet_dir is not None:
    #         self._cleanup_temp_parquet_dir()

    # def _setup_paths_and_create_data(self):
    #     """
    #     Set up paths and create parquet files for parquet backend.

    #     This method is called only when creating a new parquet store (parquet backend).
    #     It processes the mzTab file to create parquet files in the specified output directory.

    #     The method handles:
    #     - Processing mzTab file to parquet format
    #     - Creating parquet files for all mzTab sections
    #     - Adding MSstats data if provided

    #     Note:
    #         This method is not called when opening existing stores.
    #         For existing stores, the parquet files should already exist.
    #     """
    #     # This method is only called for parquet backend when creating new store
    #     if self._mztab_path is None:
    #         self.logger.info(
    #             "Opening existing database/store - skipping table creation"
    #         )
    #         return

    #     # Create parquet files from mzTab
    #     self.logger.info(
    #         f"Processing mzTab file to parquet format for parquet backend..."
    #     )
    #     (
    #         metadata_parquet,
    #         proteins_parquet,
    #         protein_details_parquet,
    #         psms_parquet,
    #     ) = self._process_mztab_to_parquet()
    #     self.logger.info("Finished processing mzTab file to parquet format.")

    #     # Handle MSstats if available
    #     if self._msstats_path:
    #         self.logger.debug(f"Adding MSstats parquet file from {self._msstats_path}")
    #         self.add_msstats_table(self._msstats_path)

    def _cleanup_temp_parquet_dir(self):
        """
        Clean up temporary directory used in DuckDB backend mode.

        This method removes the temporary parquet directory that was created
        during the processing of mzTab files for DuckDB backend. The parquet
        files are only temporary since the data is loaded into the DuckDB database.

        Note:
            This method is only called for DuckDB backend. For parquet backend,
            the parquet files are the final storage format and are not cleaned up.
        """
        if self._temp_parquet_dir and os.path.exists(self._temp_parquet_dir):
            self.logger.debug(
                f"Cleaning up temporary directory: {self._temp_parquet_dir}"
            )
            shutil.rmtree(self._temp_parquet_dir)
            self._temp_parquet_dir = None

    def _annotate_decoy_hit(
        self, protein_dict: Dict[str, str]
    ) -> Optional[Dict[str, str]]:
        """
        Annotates the decoy hit column based on protein accession if the value is null.

        A protein is considered a decoy if its accession starts with 'DECOY', 'REV', or 'RANDOM'.
        This method modifies the protein_dict in-place and returns it.

        Args:
            protein_dict: Dictionary containing protein data with accession and decoy hit fields

        Returns:
            Modified protein_dict with decoy hit annotation

        Note:
            The decoy hit column may not be present in all mzTab files or may not be selected
            for processing. In such cases, the original dictionary is returned unchanged.
        """
        decoy_col = "opt_global_cv_PRIDE:0000303_decoy_hit"

        # The column may not be present in all mzTab files, or may not be selected.
        if decoy_col not in protein_dict or (
            self._protein_columns and decoy_col not in self._protein_columns
        ):
            return protein_dict

        if protein_dict.get(decoy_col) == "null":
            accession = protein_dict.get(ACCESSION)
            if accession and accession != "null":
                accession_upper = accession.upper()
                is_decoy = any(
                    accession_upper.startswith(prefix) for prefix in DECOY_PREFIXES
                )
                protein_dict[decoy_col] = "1" if is_decoy else "0"
            else:
                # If accession is also null, we cannot determine, so we set it to 0 (not a decoy).
                protein_dict[decoy_col] = "0"
        return protein_dict

    def _process_mztab_to_parquet(self):
        """
        Process mzTab file to parquet files in batches.

        This method orchestrates the processing of mzTab files by coordinating
        setup, file processing, batch writing, and cleanup operations.

        Returns:
            tuple: A tuple containing paths to the created parquet files:
                - metadata_parquet (Path): Path to metadata.parquet
                - proteins_parquet (Path): Path to proteins.parquet
                - protein_details_parquet (Path): Path to protein_details.parquet
                - psms_parquet (Path): Path to psms.parquet

        Raises:
            FileNotFoundError: If the mzTab file doesn't exist
            PermissionError: If there are permission issues accessing the file
            MemoryError: If insufficient memory to process the file
            RuntimeError: If there are insufficient disk space or other processing errors
        """
        # Setup processing environment
        context = self._setup_processing_environment()

        try:
            # Process the mzTab file
            self._process_mztab_file(context)

            # Write final batches
            self._write_final_batches(context)

        except FileNotFoundError:
            raise FileNotFoundError(f"mzTab file not found: {self._mztab_path}")
        except PermissionError:
            raise PermissionError(
                f"Permission denied accessing file: {self._mztab_path}"
            )
        except MemoryError:
            raise MemoryError(
                "Insufficient memory to process mzTab file. Try reducing batch_size."
            )
        except Exception as e:
            self.logger.error(f"Unexpected error processing mzTab file: {e}")
            raise RuntimeError(f"Failed to process mzTab file: {e}")
        finally:
            # Clean up resources
            self._cleanup_processing_resources(context)

        # Log completion statistics
        self._log_processing_completion(context)

        return (
            context["parquet_paths"]["metadata"],
            context["parquet_paths"]["proteins"],
            context["parquet_paths"]["protein_details"],
            context["parquet_paths"]["psms"],
        )

    def _setup_processing_environment(self):
        """Setup processing environment and return context dictionary."""
        # Determine if file is gzipped
        is_gzipped = str(self._mztab_path).endswith(".gz")
        output_dir = self._get_output_dir()

        self.logger.debug(f"Processing mzTab file: {self._mztab_path}")
        self.logger.debug(f"File is gzipped: {is_gzipped}")
        self.logger.debug(f"Output directory: {output_dir}")
        self.logger.debug(f"Batch size: {self._batch_size}")

        # Check available disk space
        if not check_disk_space(output_dir):
            raise RuntimeError("Insufficient disk space for processing mzTab file")

        return {
            "is_gzipped": is_gzipped,
            "output_dir": output_dir,
            "parquet_paths": {
                "metadata": output_dir / "metadata.parquet",
                "proteins": output_dir / "proteins.parquet",
                "protein_details": output_dir / "protein_details.parquet",
                "psms": output_dir / "psms.parquet",
            },
            "writers": {"proteins": None, "protein_details": None, "psms": None},
            "counters": {
                "metadata": 0,
                "protein": 0,
                "protein_details": 0,
                "psm": 0,
                "batch": 0,
            },
            "data_buffers": {
                "metadata": [],
                "proteins": [],
                "protein_details": [],
                "psms": [],
            },
        }

    def _process_mztab_file(self, context):
        """Process the mzTab file line by line."""
        open_func = gzip.open if context["is_gzipped"] else open

        with open_func(self._mztab_path, "rt") as f:
            for line_num, line in enumerate(f, 1):
                try:
                    line = line.strip()
                    if not line:
                        continue

                    self._process_single_line(line, context, line_num)
                    self._write_batches_if_needed(context)

                except Exception as e:
                    self.logger.error(f"Error processing line {line_num}: {e}")
                    self.logger.error(f"Line content: {line[:100]}...")
                    continue

    def _process_single_line(self, line, context, line_num):
        """Process a single mzTab line."""
        if is_mztab_line_type(line, "MTD"):
            self._process_metadata_line(line, context)
        elif is_mztab_line_type(line, "PRH"):
            self._process_protein_header_line(line)
        elif is_mztab_line_type(line, "PRT"):
            self._process_protein_line(line, context)
        elif is_mztab_line_type(line, "PSH"):
            self._process_psm_header_line(line)
        elif is_mztab_line_type(line, "PSM"):
            self._process_psm_line(line, context)

    def _process_metadata_line(self, line, context):
        """Process metadata line."""
        context["data_buffers"]["metadata"].append(parse_metadata_line(line))
        context["counters"]["metadata"] += 1

    def _process_protein_header_line(self, line):
        """Process protein header line."""
        self._protein_header = parse_header_line(line)
        self.logger.debug(
            f"Found protein header with {len(self._protein_header)} columns"
        )

    def _process_protein_line(self, line, context):
        """Process protein line."""
        protein_dict = self._parse_protein_line_with_processing(line)
        if protein_dict is None:
            return

        if OPT_GLOBAL_RESULT_TYPE in protein_dict:
            result_type = protein_dict[OPT_GLOBAL_RESULT_TYPE]

            if result_type in [INDISTINGUISHABLE_GROUP, SINGLE_PROTEIN_MZTAB]:
                context["data_buffers"]["proteins"].append(protein_dict)
                context["counters"]["protein"] += 1
            elif result_type == PROTEIN_DETAILS_MZTAB:
                context["data_buffers"]["protein_details"].append(protein_dict)
                context["counters"]["protein_details"] += 1

    def _process_psm_header_line(self, line):
        """Process PSM header line."""
        self._psm_header = parse_header_line(line)
        self.logger.debug(f"Found PSM header with {len(self._psm_header)} columns")

    def _process_psm_line(self, line, context):
        """Process PSM line."""
        context["data_buffers"]["psms"].append(
            self._parse_psm_line_with_processing(line)
        )
        context["counters"]["psm"] += 1

    def _write_batches_if_needed(self, context):
        """Write batches to parquet files if they've reached the batch size limit."""
        # Write protein batch if needed
        if len(context["data_buffers"]["proteins"]) >= self._batch_size:
            context["writers"]["proteins"] = self._write_protein_batch(
                context["data_buffers"]["proteins"],
                context["parquet_paths"]["proteins"],
                context["writers"]["proteins"],
            )
            context["counters"]["batch"] += 1
            self.logger.debug(
                f"Wrote protein batch {context['counters']['batch']} "
                f"({len(context['data_buffers']['proteins'])} proteins)"
            )
            context["data_buffers"]["proteins"] = []

        # Write protein details batch if needed
        if len(context["data_buffers"]["protein_details"]) >= self._batch_size:
            context["writers"]["protein_details"] = self._write_protein_details_batch(
                context["data_buffers"]["protein_details"],
                context["parquet_paths"]["protein_details"],
                context["writers"]["protein_details"],
            )
            self.logger.debug(
                f"Wrote protein_details batch ({len(context['data_buffers']['protein_details'])} entries)"
            )
            context["data_buffers"]["protein_details"] = []

        # Write PSM batch if needed
        if len(context["data_buffers"]["psms"]) >= self._batch_size:
            context["writers"]["psms"] = self._write_psm_batch(
                context["data_buffers"]["psms"],
                context["parquet_paths"]["psms"],
                context["writers"]["psms"],
            )
            context["counters"]["batch"] += 1
            self.logger.debug(
                f"Wrote PSM batch {context['counters']['batch']} ({len(context['data_buffers']['psms'])} PSMs)"
            )
            context["data_buffers"]["psms"] = []

    def _write_final_batches(self, context):
        """Write any remaining data batches to parquet files."""
        if context["data_buffers"]["metadata"]:
            self._write_metadata_batch(
                context["data_buffers"]["metadata"],
                context["parquet_paths"]["metadata"],
            )
            self.logger.debug(
                f"Wrote metadata batch ({len(context['data_buffers']['metadata'])} entries)"
            )

        if context["data_buffers"]["proteins"]:
            context["writers"]["proteins"] = self._write_protein_batch(
                context["data_buffers"]["proteins"],
                context["parquet_paths"]["proteins"],
                context["writers"]["proteins"],
            )
            self.logger.debug(
                f"Wrote final protein batch ({len(context['data_buffers']['proteins'])} proteins)"
            )

        if context["data_buffers"]["protein_details"]:
            context["writers"]["protein_details"] = self._write_protein_details_batch(
                context["data_buffers"]["protein_details"],
                context["parquet_paths"]["protein_details"],
                context["writers"]["protein_details"],
            )
            self.logger.debug(
                f"Wrote final protein_details batch ({len(context['data_buffers']['protein_details'])} entries)"
            )

        if context["data_buffers"]["psms"]:
            context["writers"]["psms"] = self._write_psm_batch(
                context["data_buffers"]["psms"],
                context["parquet_paths"]["psms"],
                context["writers"]["psms"],
            )
            self.logger.debug(
                f"Wrote final PSM batch ({len(context['data_buffers']['psms'])} PSMs)"
            )

    def _cleanup_processing_resources(self, context):
        """Close all parquet writers safely."""
        try:
            for writer_name, writer in context["writers"].items():
                if writer:
                    try:
                        writer.close()
                        self.logger.debug(f"Successfully closed {writer_name} writer")
                    except Exception as writer_error:
                        self.logger.warning(
                            f"Error closing {writer_name} writer: {writer_error}"
                        )
        except Exception as e:
            self.logger.error(f"Error closing parquet writers: {e}")

    def _log_processing_completion(self, context):
        """Log processing completion statistics."""
        counters = context["counters"]
        self.logger.debug(
            f"Processing complete. Total counts: metadata={counters['metadata']}, "
            f"proteins={counters['protein']}, protein_details={counters['protein_details']}, "
            f"psms={counters['psm']}"
        )

    def _parse_protein_line_with_processing(
        self, line: str
    ) -> Optional[Dict[str, str]]:
        """
        Parse protein line into dictionary with additional processing.

        Uses the generic parse_protein_line utility function and adds
        protein group processing and decoy annotation specific to this class.

        Args:
            line: Tab-separated protein line from mzTab file

        Returns:
            Dictionary with protein data, including processed protein groups and decoy annotation

        Note:
            This method assumes self._protein_header has been set by a previous PRH line.
        """
        try:
            # Validate input
            if not line or not line.strip():
                raise ValueError("Empty or invalid protein line")

            if not hasattr(self, "_protein_header") or not self._protein_header:
                raise ValueError(
                    "Protein header not found. Ensure PRH line is processed first."
                )

            protein_dict = parse_protein_line(line, self._protein_header)

            # Validate required fields
            if not protein_dict:
                raise ValueError("Failed to parse protein line into dictionary")

            # Validate accession field
            if ACCESSION not in protein_dict:
                self.logger.warning(
                    f"Protein line missing accession field: {line[:100]}..."
                )
                return None

            if OPT_GLOBAL_RESULT_TYPE in protein_dict:
                result_type = protein_dict[OPT_GLOBAL_RESULT_TYPE]
                if result_type == INDISTINGUISHABLE_GROUP:
                    anchor = protein_dict.get(ACCESSION)
                    ambiguity_members = protein_dict.get(AMBIGUITY_MEMBERS)
                    if ambiguity_members and ambiguity_members != "null":
                        try:
                            standardized_accession = (
                                standardize_protein_string_accession(
                                    ambiguity_members, is_sorted=True
                                )
                            )
                            protein_dict[ACCESSION] = standardized_accession
                        except Exception as e:
                            self.logger.warning(
                                f"Failed to standardize protein accession: {e}"
                            )
                            # Keep original accession if standardization fails
                    else:
                        self.logger.warning(
                            f"Indistinguishable group {anchor} has no ambiguity members."
                        )
                    protein_dict[ANCHOR_PROTEIN] = anchor
                elif result_type == SINGLE_PROTEIN_MZTAB:
                    protein_dict[ANCHOR_PROTEIN] = protein_dict.get(ACCESSION)

            protein_dict = self._annotate_decoy_hit(protein_dict)
            return protein_dict

        except Exception as e:
            self.logger.error(f"Error parsing protein line: {e}")
            self.logger.error(f"Line content: {line[:100]}...")
            # Return a minimal valid dictionary to prevent processing failure
            return {ACCESSION: "error", "error": str(e)}

    def _parse_psm_line_with_processing(self, line: str) -> Dict[str, str]:
        """
        Parse PSM line into dictionary with additional processing.

        Uses the generic parse_psm_line utility function and adds
        protein accession standardization and spectra reference parsing.

        Args:
            line: Tab-separated PSM line from mzTab file

        Returns:
            Dictionary with PSM data, including standardized protein accessions
            and parsed spectra reference information

        Note:
            This method assumes self._psm_header has been set by a previous PSH line.
        """
        try:
            # Validate input
            if not line or not line.strip():
                raise ValueError("Empty or invalid PSM line")

            if not hasattr(self, "_psm_header") or not self._psm_header:
                raise ValueError(
                    "PSM header not found. Ensure PSH line is processed first."
                )

            psm_dict = parse_psm_line(line, self._psm_header)

            # Validate required fields
            if not psm_dict:
                raise ValueError("Failed to parse PSM line into dictionary")

            # Validate accession field
            if ACCESSION not in psm_dict:
                self.logger.warning(
                    f"PSM line missing accession field: {line[:100]}..."
                )
                psm_dict[ACCESSION] = "unknown"

            if (
                ACCESSION in psm_dict
                and psm_dict[ACCESSION]
                and psm_dict[ACCESSION] != "null"
            ):
                try:
                    psm_dict[ACCESSION] = standardize_protein_string_accession(
                        psm_dict[ACCESSION], is_sorted=True
                    )
                except Exception as e:
                    self.logger.warning(
                        f"Failed to standardize PSM protein accession: {e}"
                    )
                    # Keep original accession if standardization fails

            # Split spectra_ref into file and scan info for plex experiments
            if SPECTRA_REF in psm_dict and psm_dict[SPECTRA_REF]:
                try:
                    spectra_ref_parts = psm_dict[SPECTRA_REF].split(":")
                    psm_dict[SPECTRA_REF_FILE] = spectra_ref_parts[0]
                    if len(spectra_ref_parts) > 1:
                        psm_dict[SPECTRA_REF_SCAN] = generate_scan_number(
                            spectra_ref_parts[1]
                        )
                    else:
                        psm_dict[SPECTRA_REF_SCAN] = None
                except Exception as e:
                    self.logger.warning(f"Failed to parse spectra reference: {e}")
                    psm_dict[SPECTRA_REF_FILE] = psm_dict.get(SPECTRA_REF, "unknown")
                    psm_dict[SPECTRA_REF_SCAN] = None

            return psm_dict

        except Exception as e:
            self.logger.error(f"Error parsing PSM line: {e}")
            self.logger.error(f"Line content: {line[:100]}...")
            # Return a minimal valid dictionary to prevent processing failure
            return {ACCESSION: "error", "error": str(e)}

    def _write_metadata_batch(
        self, metadata: List[Dict[str, Optional[str]]], output_path: Path
    ) -> None:
        """
        Write metadata batch to parquet file.

        Args:
            metadata: List of metadata dictionaries to write
            output_path: Path to the output parquet file
        """
        df = pd.DataFrame(metadata)
        df.to_parquet(output_path, index=False)

    def _write_protein_batch(
        self, proteins: List[Dict[str, str]], output_path: Path, writer=None
    ) -> Optional[pq.ParquetWriter]:
        """
        Write protein batch to parquet file.

        Args:
            proteins: List of protein dictionaries to write
            output_path: Path to the output parquet file
            writer: Existing ParquetWriter instance or None to create new one

        Returns:
            ParquetWriter instance for continued writing
        """
        if not proteins:
            return writer
        df = pd.DataFrame(proteins)
        if self._protein_columns:
            # Only select columns that exist in the DataFrame
            cols = [col for col in self._protein_columns if col in df.columns]
            df = df[cols]
        table = pa.Table.from_pandas(df)
        if writer is None:
            writer = pq.ParquetWriter(output_path, table.schema)
        writer.write_table(table)
        return writer

    def _write_protein_details_batch(
        self, protein_details: List[Dict[str, str]], output_path: Path, writer=None
    ) -> Optional[pq.ParquetWriter]:
        """
        Write protein details batch to parquet file.

        Args:
            protein_details: List of protein details dictionaries to write
            output_path: Path to the output parquet file
            writer: Existing ParquetWriter instance or None to create new one

        Returns:
            ParquetWriter instance for continued writing
        """
        df = pd.DataFrame(protein_details)
        table = pa.Table.from_pandas(df)

        if writer is None:
            # First batch or new writer needed
            writer = pq.ParquetWriter(output_path, table.schema)

        writer.write_table(table)
        return writer

    def _write_psm_batch(
        self, psms: List[Dict[str, str]], output_path: Path, writer=None
    ) -> Optional[pq.ParquetWriter]:
        """
        Write PSM batch to parquet file.

        Args:
            psms: List of PSM dictionaries to write
            output_path: Path to the output parquet file
            writer: Existing ParquetWriter instance or None to create new one

        Returns:
            ParquetWriter instance for continued writing
        """
        df = pd.DataFrame(psms)

        # Filter columns if specified
        if self._psm_columns:
            df = df.reindex(columns=self._psm_columns)

        table = pa.Table.from_pandas(df)

        if writer is None:
            # First batch or new writer needed
            writer = pq.ParquetWriter(output_path, table.schema)

        writer.write_table(table)
        return writer

    def create_database_tables(self):
        """
        Create database tables from mzTab file processing.

        This method processes the mzTab file and creates the appropriate storage
        based on the backend type:

        - For DuckDB backend: Creates parquet files temporarily, then loads them
          into DuckDB tables and cleans up the temporary files
        - For Parquet backend: Creates parquet files as the final storage format

        The method handles all mzTab sections:
        - Metadata (MTD lines)
        - Proteins (PRT lines)
        - Protein details (PRT lines with specific result types)
        - PSMs (PSM lines)
        - MSstats (if provided)

        Raises:
            RuntimeError: If processing fails or required files are missing

        Note:
            This method is called during initialization when creating a new database/store.
            For existing databases/stores, this method is not called.
        """
        self.logger.info(f"Processing mzTab file to parquet format...")
        (
            metadata_parquet,
            proteins_parquet,
            protein_details_parquet,
            psms_parquet,
        ) = self._process_mztab_to_parquet()
        self.logger.info("Finished processing mzTab file to parquet format.")

        # If using duckdb backend, ingest parquet files into the database
        self.logger.debug("Using DuckDB backend - creating tables from parquet files")
        tables = self._duckdb.execute("SHOW TABLES").fetchall()
        table_names = [table[0] for table in tables]
        self.logger.debug(f"Existing tables: {table_names}")

        if self._MZTAB_INDEXER_TABLE_METADATA not in table_names:
            self.logger.info(f"[Creating table] {self._MZTAB_INDEXER_TABLE_METADATA}")
            self._duckdb.execute(
                f"CREATE TABLE {self._MZTAB_INDEXER_TABLE_METADATA} AS SELECT * FROM read_parquet('{metadata_parquet}')"
            )
            self.logger.debug(f"Created metadata table from {metadata_parquet}")

        if (
            self._MZTAB_INDEXER_TABLE_PROTEINS not in table_names
            and os.path.exists(proteins_parquet)
            and os.path.getsize(proteins_parquet) > 0
        ):
            self.logger.info(f"[Creating table] {self._MZTAB_INDEXER_TABLE_PROTEINS}")
            self._duckdb.execute(
                f"CREATE TABLE {self._MZTAB_INDEXER_TABLE_PROTEINS} AS SELECT * FROM read_parquet('{proteins_parquet}')"
            )
            self.logger.debug(f"Created proteins table from {proteins_parquet}")

        # if (
        #     self._MZTAB_INDEXER_TABLE_PROTEIN_DETAILS not in table_names
        #     and os.path.exists(protein_details_parquet)
        #     and os.path.getsize(protein_details_parquet) > 0
        # ):
        #     self.logger.info(
        #         f"[Creating table] {self._MZTAB_INDEXER_TABLE_PROTEIN_DETAILS}"
        #     )
        #     self._duckdb.execute(
        #         f"CREATE TABLE {self._MZTAB_INDEXER_TABLE_PROTEIN_DETAILS} AS SELECT * FROM read_parquet('{protein_details_parquet}')"
        #     )
        #     self.logger.debug(
        #         f"Created protein_details table from {protein_details_parquet}"
        #     )

        if (
            self._MZTAB_INDEXER_TABLE_PSMS not in table_names
            and os.path.exists(psms_parquet)
            and os.path.getsize(psms_parquet) > 0
        ):
            self.logger.info(f"[Creating table] {self._MZTAB_INDEXER_TABLE_PSMS}")
            self._duckdb.execute(
                f"CREATE TABLE {self._MZTAB_INDEXER_TABLE_PSMS} AS SELECT * FROM read_parquet('{psms_parquet}')"
            )
            self.logger.debug(f"Created psms table from {psms_parquet}")

        if self._msstats_path:
            self.logger.debug(f"Adding MSstats table from {self._msstats_path}")
            self.add_msstats_table(self._msstats_path)

        # Clean up temporary parquet files for DuckDB backend
        if self._temp_parquet_dir is not None:
            self._cleanup_temp_parquet_dir()

    def _get_table_source(self, table_name: str) -> str:
        """
        Get the source for a table in a query based on backend type.

        This method returns the appropriate table source string for SQL queries
        based on the current backend configuration.

        Args:
            table_name: Name of the table to get source for

        Returns:
            Table source string:
            - For DuckDB backend: Returns the table name directly
            - For Parquet backend: Returns "read_parquet('path/to/file.parquet')"

        Raises:
            FileNotFoundError: If required parquet file doesn't exist (parquet backend)
            ValueError: If backend type is unsupported

        Note:
            This method abstracts the backend differences so that query methods
            can work with both backends without modification.
        """
        return table_name

    def get_metadata(self) -> pd.DataFrame:
        """
        Get metadata as DataFrame.

        Returns:
            DataFrame containing all metadata entries from the mzTab file.
            The DataFrame has columns for 'key' and 'value' pairs.

        Note:
            Metadata includes information about the experiment, instruments,
            software versions, and other configuration details.
        """
        if self._metadata_cache is None:
            source = self._get_table_source(self._MZTAB_INDEXER_TABLE_METADATA)
            self._metadata_cache = self.query_to_df(f"SELECT * FROM {source}")
        return self._metadata_cache

    def _is_protein_score_qvalue(self) -> bool:
        """
        Check if the protein score is a q-value.

        Highly optimized version using vectorized operations and regex for maximum performance.
        """
        metadata_df = self.get_metadata()

        # Filter rows that contain 'protein_search_engine_score' in the key
        protein_score_rows = metadata_df[
            metadata_df["key"].str.contains(
                "protein_search_engine_score", case=False, na=False
            )
        ]

        if protein_score_rows.empty:
            return False

        # Convert values to lowercase once and check all words at once
        values_lower = protein_score_rows["value"].astype(str).str.lower()

        # Check if any row contains all required words using vectorized operations
        # This is much faster than iterating through each row
        for word in OPENMS_PROTEIN_QVALUE_WORDS:
            values_lower = values_lower[values_lower.str.contains(word, na=False)]
            if values_lower.empty:
                return False

        return True

    def get_proteins(self) -> pd.DataFrame:
        """
        Get proteins as DataFrame.

        Returns:
            DataFrame containing all protein entries from the mzTab file.
            Includes protein accessions, descriptions, scores, and other
            protein-level information.

        Note:
            This includes both single proteins and protein groups (indistinguishable
            groups). Protein details are stored separately.
        """
        source = self._MZTAB_INDEXER_TABLE_PROTEINS
        return self.query_to_df(f"SELECT * FROM {source}")

    def get_protein_details(self) -> pd.DataFrame:
        """
        Get protein_details as DataFrame.

        Returns:
            DataFrame containing all protein entries from the mzTab file.
            Includes protein accessions, descriptions, scores, and other
            protein-level information.

        Note:
            This includes both single proteins and protein groups (indistinguishable
            groups). Protein details are stored separately.
        """
        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_PROTEIN_DETAILS)
        return self.query_to_df(f"SELECT * FROM {source}")

    def get_psms(self) -> pd.DataFrame:
        """
        Get psms as DataFrame.
        """
        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_PSMS)
        return self.query_to_df(f"SELECT * FROM {source}")

    def get_protein_best_searchengine_score(self) -> pd.DataFrame:
        """
        Get protein q-value as DataFrame.
        """
        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_PROTEINS)
        return self.query_to_df(
            f'SELECT {ACCESSION}, "{MZTAB_PROTEIN_BEST_SEARCH_ENGINE_SCORE}" AS opt_global_qvalue FROM {source}'
        )

    def get_protein_details(self) -> pd.DataFrame:
        """
        Get protein details as DataFrame.

        Returns:
            DataFrame containing protein detail entries from the mzTab file.
            These are typically individual proteins within protein groups
            that have additional detailed information.

        Note:
            Protein details are separate from main protein entries and contain
            more granular information about individual proteins within groups.
        """
        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_PROTEIN_DETAILS)
        return self.query_to_df(f"SELECT * FROM {source}")

    def _get_num_psms(self) -> int:
        """
        Get number of PSMs in the database.
        """
        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_PSMS)
        return self.query_to_df(f"SELECT COUNT(*) as count FROM {source}")[
            "count"
        ].iloc[0]

    def get_msstats(self) -> Optional[pd.DataFrame]:
        """
        Get MSstats data as DataFrame if available.

        Returns:
            DataFrame containing MSstats data if available, None otherwise.
            Includes protein names, peptide sequences, intensities, and
            other quantitative information.

        Note:
            MSstats data is optional and may not be present in all databases.
            This data is used for quantitative analysis and statistical testing.
        """
        if self._has_msstats_table():
            source = self._get_table_source(self._MZTAB_INDEXER_TABLE_MSSTATS)
            return self.query_to_df(f"SELECT * FROM {source}")
        return None

    def get_protein_by_accession(self, accession: str) -> pd.DataFrame:
        """Get protein data by accession.

        Args:
            accession: Protein accession

        Returns:
            DataFrame with protein data
        """
        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_PROTEINS)
        query = f"""
        SELECT *
        FROM {source}
        WHERE {ACCESSION} = ?
        """
        return self.secure_query_to_df(query, {"accession": accession})

    def get_psms_by_protein(self, accession: str) -> pd.DataFrame:
        """Get PSMs for a specific protein.

        Args:
            accession: Protein accession

        Returns:
            DataFrame with PSM data
        """
        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_PSMS)
        query = f"""
        SELECT *
        FROM {source}
        WHERE {ACCESSION} LIKE ?
        """
        return self.secure_query_to_df(query, {"accession": f"%{accession}%"})

    def get_metadata_value(self, key: str) -> Optional[str]:
        """Get metadata value by key.

        Args:
            key: Metadata key

        Returns:
            Metadata value if found, None otherwise
        """
        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_METADATA)
        query = f"""
        SELECT value
        FROM {source}
        WHERE key = ?
        LIMIT 1
        """
        result = self.secure_query_to_df(query, {"key": key})
        return result["value"].iloc[0] if not result.empty else None

    def get_protein_count(self) -> int:
        """
        Get total number of proteins in the database.

        Returns:
            Total count of protein entries including both single proteins
            and protein groups (indistinguishable groups).

        Note:
            This count includes all protein entries regardless of their
            result type or identification status.
        """
        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_PROTEINS)
        query = f"SELECT COUNT(*) as count FROM {source}"
        result = self.query_to_df(query)
        return result["count"].iloc[0]

    def get_protein_details_count(self) -> int:
        """
        Get total number of protein details in the database.

        Returns:
            Total count of protein detail entries. These are typically
            individual proteins within protein groups that have additional
            detailed information.

        Note:
            Protein details are separate from main protein entries and
            contain more granular information about individual proteins.
        """
        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_PROTEIN_DETAILS)
        query = f"SELECT COUNT(*) as count FROM {source}"
        result = self.query_to_df(query)
        return result["count"].iloc[0]

    def get_psm_count(self) -> int:
        """
        Get total number of PSMs (Peptide-Spectrum Matches) in the database.

        Returns:
            Total count of PSM entries representing individual peptide
            identifications with their associated spectral evidence.

        Note:
            PSMs are the most granular level of identification and include
            all peptide-spectrum matches regardless of their quality scores.
        """
        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_PSMS)
        query = f"SELECT COUNT(*) as count FROM {source}"
        result = self.query_to_df(query)
        return result["count"].iloc[0]

    def stream_section(self, section: str, chunk_size: int = 1000000):
        """
        Stream data from a specific section in chunks for memory-efficient processing.

        This method allows processing large datasets without loading everything
        into memory at once. It uses SQL LIMIT and OFFSET to paginate through
        the data in configurable chunk sizes.

        Args:
            section: Section to stream. Valid values:
                - "PSM": Peptide-Spectrum Matches
                - "PRT": Proteins (main protein entries)
                - "PEP": Peptides (uses PSM data as proxy)
            chunk_size: Number of rows to process in each chunk (default: 1,000,000)

        Yields:
            DataFrame containing data from the specified section for each chunk

        Raises:
            ValueError: If section name is invalid

        Note:
            This method is useful for processing very large datasets that don't
            fit in memory. Each yielded DataFrame contains up to chunk_size rows.
            Empty DataFrames indicate the end of the data.
        """
        section = section.upper()

        if section == "PSM":
            source = self._get_table_source(self._MZTAB_INDEXER_TABLE_PSMS)
        elif section == "PRT":
            source = self._get_table_source(self._MZTAB_INDEXER_TABLE_PROTEINS)
        elif section == "PEP":
            # For PEP, we'll use PSM data since PEP is typically derived from PSM
            source = self._get_table_source(self._MZTAB_INDEXER_TABLE_PSMS)
        else:
            raise ValueError(f"Invalid section: {section}")

        offset = 0
        while True:
            query = f"SELECT * FROM {source} LIMIT {chunk_size} OFFSET {offset}"
            df = self.query_to_df(query)
            if df.empty:
                break
            yield df
            offset += chunk_size

    def _has_msstats_table(self) -> bool:
        """
        Check if MSstats table exists in the database.

        Returns:
            True if MSstats table exists, False otherwise
        """
        try:
            # For DuckDB backend, check if table exists in database
            tables = self._duckdb.execute("SHOW TABLES").fetchall()
            table_names = [table[0] for table in tables]
            return self._MZTAB_INDEXER_TABLE_MSSTATS in table_names
        except Exception:
            return False

    def get_msstats_by_peptide(self, peptide: str) -> Optional[pd.DataFrame]:
        """Get MSstats data for a specific peptide.

        Args:
            peptide: Peptide sequence

        Returns:
            DataFrame with MSstats data if available
        """
        if not self._has_msstats_table():
            return None
        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_MSSTATS)
        query = f"""
        SELECT *
        FROM {source}
        WHERE {MSSTATS_PEPTIDE_SEQUENCE} = ?
        """
        return self.secure_query_to_df(query, {"peptide": peptide})

    def get_msstats_runs(self) -> Optional[List[str]]:
        """Get list of runs from MSstats data.

        Returns:
            List of run names if MSstats data is available
        """
        if not self._has_msstats_table():
            return None

        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_MSSTATS)
        query = f"""
        SELECT DISTINCT {MSSTATS_REFERENCE}
        FROM {source}
        ORDER BY {MSSTATS_REFERENCE}
        """
        result = self.query_to_df(query)
        return result[MSSTATS_REFERENCE].tolist() if not result.empty else []

    def get_msstats_protein_intensities(self, protein: str) -> Optional[pd.DataFrame]:
        """Get protein intensities across runs from MSstats data.

        Args:
            protein: Protein accession or name

        Returns:
            DataFrame with protein intensities if available
        """
        if not self._has_msstats_table():
            return None

        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_MSSTATS)
        query = f"""
        SELECT {MSSTATS_REFERENCE}, {MSSTATS_INTENSITY}
        FROM {source}
        WHERE {MSSTATS_PROTEIN_NAME} LIKE ?
        ORDER BY {MSSTATS_REFERENCE}
        """
        return self.secure_query_to_df(query, {"protein": f"%{protein}%"})

    def add_msstats_table(self, msstats_path: str):
        """
        Add or replace MSstats table in the database.

        This method processes an MSstats CSV file and adds it to the database
        based on the backend type:

        - For DuckDB backend: Creates a table in the database and adds indices
        - For Parquet backend: Creates a parquet file in the database directory

        The method handles:
        - Processing the MSstats file to parquet format with transformations
        - Dropping existing MSstats table if it exists (DuckDB backend)
        - Creating optimized indices for common query patterns
        - Data transformations (protein name standardization, reference name cleaning)

        Args:
            msstats_path: Path to the MSstats CSV file to process

        Raises:
            RuntimeError: If database is not initialized
            FileNotFoundError: If MSstats file doesn't exist
            RuntimeError: If processing fails

        Note:
            This method can be called multiple times to update MSstats data.
            Each call will replace the existing MSstats data completely.
        """
        if not self._duckdb:
            raise RuntimeError("Database not initialized.")

        self.logger.info("Attempting to add/replace msstats data.")
        self._msstats_path = msstats_path
        msstats_parquet = self._process_msstats_to_parquet()

        # Check if the msstats table already exists and drop it
        tables = self._duckdb.execute("SHOW TABLES").fetchall()
        table_names = [table[0] for table in tables]
        if self._MZTAB_INDEXER_TABLE_MSSTATS in table_names:
            self.logger.debug("Dropping existing msstats table")
            self._duckdb.execute(f"DROP TABLE {self._MZTAB_INDEXER_TABLE_MSSTATS}")

        # Create table from parquet
        self._duckdb.execute(
            f"CREATE TABLE {self._MZTAB_INDEXER_TABLE_MSSTATS} AS SELECT * FROM read_parquet('{msstats_parquet}')"
        )
        self.logger.info("[Creating table] msstats")

        # Create indices
        table_columns = [
            col[0]
            for col in self._duckdb.execute(
                f"SELECT name FROM pragma_table_info('{self._MZTAB_INDEXER_TABLE_MSSTATS}')"
            ).fetchall()
        ]
        columns_to_index = [MSSTATS_PROTEIN_NAME]
        for column in columns_to_index:
            if column in table_columns:
                safe_column_name = column.replace(".", "_")
                self._duckdb.execute(
                    f'CREATE INDEX IF NOT EXISTS idx_msstats_{safe_column_name} ON {self._MZTAB_INDEXER_TABLE_MSSTATS} ("{column}")'
                )

    def _process_msstats_to_parquet(self) -> Path:
        """
        Process MSstats file to parquet format with transformations.

        This method reads an MSstats CSV file in chunks and applies data transformations
        before writing to a parquet file. The processing is done in chunks to handle
        large files efficiently and manage memory usage.

        The method applies the following transformations to the MSstats data:
        - Protein name standardization: Sorts and joins protein names separated by semicolons
        - Reference name creation: Removes file extensions and controller information from
          reference names to create cleaner identifiers

        The method uses optimized parquet writing settings:
        - ZSTD compression for better compression ratio
        - Dictionary encoding for string columns to reduce file size
        - Configurable compression level for balance between size and speed

        The method includes comprehensive error handling for:
        - File corruption and malformed data
        - Memory issues during processing
        - Resource cleanup (parquet writer)

        Returns:
            Path: Path to the created msstats.parquet file

        Raises:
            FileNotFoundError: If the MSstats file doesn't exist
            PermissionError: If there are permission issues accessing the file
            MemoryError: If insufficient memory to process the file
            RuntimeError: If there are processing errors during transformation or writing

        Note:
            - The method uses chunked processing with configurable batch size (self._batch_size)
            - Progress is logged at DEBUG level for monitoring large file processing
            - The parquet writer is properly closed even if errors occur
            - Memory usage is optimized by processing data in chunks rather than loading
              the entire file into memory at once
            - The method assumes the MSstats file is in CSV format with standard columns
        """
        self.logger.debug(
            f"Starting MSstats processing to parquet: {self._msstats_path}"
        )
        output_dir = self._get_output_dir()
        msstats_parquet = output_dir / "msstats.parquet"
        writer = None
        total_rows = 0

        try:
            # Process file in chunks to avoid memory issues
            self.logger.debug(
                f"Reading MSstats file in chunks of size {self._batch_size}"
            )
            for chunk_idx, df in enumerate(
                pd.read_csv(self._msstats_path, chunksize=self._batch_size)
            ):
                chunk_size = len(df)
                total_rows += chunk_size

                # Sort and join protein names
                if MSSTATS_PROTEIN_NAME in df.columns:
                    df[MSSTATS_PROTEIN_NAME] = df[MSSTATS_PROTEIN_NAME].apply(
                        lambda x: (
                            ";".join(sorted(str(x).split(";"))) if pd.notna(x) else x
                        )
                    )

                # Create Reference_Name by removing file extensions and controller info
                if MSSTATS_REFERENCE in df.columns:
                    df[MSSTATS_REFERENCE_NAME] = df[MSSTATS_REFERENCE].apply(
                        lambda x: (
                            x.split("_controllerType")[0].rsplit(".", 1)[0]
                            if isinstance(x, str)
                            else x
                        )
                    )

                # Match sample_accession and channel from the SDRF.
                df.rename(columns=MSSTATS_NAME_MAP, inplace=True)

                if self.experiment_type == "LFQ":
                    df["file_channel"] = df[MSSTATS_REFERENCE_NAME]
                    df["channel"] = "LFQ"
                else:
                    channel_mapping = self._create_channel_mapping_table()
                    df["channel"] = df["Channel"].map(channel_mapping)
                    df.drop(columns=["Channel"], inplace=True)
                    df.loc[:, "file_channel"] = (
                        df[MSSTATS_REFERENCE_NAME] + "-" + df["channel"]
                    )

                df["sample_accession"] = df["file_channel"].map(self._sample_map)

                if "RetentionTime" in df.columns:
                    df.rename(columns={"RetentionTime": "rt"}, inplace=True)
                    need_cols = MSSTATS_USECOLS + ["rt"]
                else:
                    need_cols = MSSTATS_USECOLS

                # Convert to pyarrow table
                table = pa.Table.from_pandas(df[need_cols])

                if writer is None:
                    # First chunk - create writer with schema
                    writer = pq.ParquetWriter(
                        msstats_parquet,
                        table.schema,
                        compression="ZSTD",  # Use ZSTD compression for better compression ratio
                        compression_level=3,  # Moderate compression level for balance
                        use_dictionary=True,  # Enable dictionary encoding for strings
                        dictionary_pagesize_limit=1048576,  # 1MB dictionary page size
                    )

                # Write chunk to parquet
                writer.write_table(table)

                # Log progress
                self.logger.debug(f"Processed chunk {chunk_idx + 1} ({len(df)} rows)")
        except Exception as e:
            self.logger.error(f"Error processing MSstats file: {str(e)}")
            if writer:
                writer.close()
            raise
        finally:
            if writer:
                self.logger.debug("Closing parquet writer")
                writer.close()

        self.logger.info(
            f"Completed MSstats processing. Total rows processed: {total_rows}"
        )

        return msstats_parquet

    def _create_channel_mapping_table(self):
        """Create a channel mapping table in DuckDB for efficient lookups."""
        if "TMT" in self.experiment_type:
            channels = TMT_CHANNELS[self.experiment_type]
        else:
            channels = ITRAQ_CHANNEL[self.experiment_type]

        # Create mapping data
        mapping_data = {i + 1: channel for i, channel in enumerate(channels)}

        return mapping_data

    # ============================================================================
    # ADVANCED MSSTATS ANALYSIS METHODS
    # ============================================================================

    def get_msstats_experiment_type(self) -> Optional[str]:
        """
        Determine the experiment type from MSstats data.

        Returns:
            Experiment type ('LFQ', 'TMT10plex', 'TMT11plex', etc.) or None
        """
        if not self._has_msstats_table():
            self.logger.warning("No MSstats data available")
            return None

        msstats_df = self.get_msstats()
        if msstats_df is None or msstats_df.empty:
            return None

        # Check if Channel column exists and what values it contains
        if "channel" in msstats_df.columns:
            unique_channels = set(msstats_df["channel"].dropna().unique())

            # Check LFQ
            if len(unique_channels) == 1 and "LFQ" in unique_channels:
                return "LFQ"

            # Check for TMT channels
            for tmt_type, channels in TMT_CHANNELS.items():
                if set(channels).issubset(unique_channels) or unique_channels.issubset(
                    set(channels)
                ):
                    return tmt_type

            # Check for iTRAQ channels
            for itraq_type, channels in ITRAQ_CHANNEL.items():
                if set(channels).issubset(unique_channels) or unique_channels.issubset(
                    set(channels)
                ):
                    return itraq_type

        # Default to LFQ if no multiplexed channels found
        return "LFQ"

    def get_msstats_file_statistics(self) -> Optional[pd.DataFrame]:
        """
        Get comprehensive file-level statistics from MSstats data using SQL.

        Returns:
            DataFrame with file-level statistics or None if no MSstats data
        """
        if not self._has_msstats_table():
            self.logger.warning("No MSstats data available")
            return None

        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_MSSTATS)

        experiment_type = self.get_msstats_experiment_type()

        if experiment_type == "LFQ":
            sql = f"""
            SELECT 
                reference_file_name,
                channel,
                COUNT(*) as feature_count,
                COUNT(DISTINCT pg_accessions) as protein_count,
                COUNT(DISTINCT peptidoform) as peptide_count,
                AVG(intensity) as avg_intensity,
                MIN(intensity) as min_intensity,
                MAX(intensity) as max_intensity,
                COUNT(DISTINCT reference_file_name) as sample_count
            FROM {source}
            GROUP BY reference_file_name, channel
            ORDER BY reference_file_name, channel
            """
        else:
            sql = f"""
            SELECT 
                reference_file_name,
                channel,
                COUNT(*) as feature_count,
                COUNT(DISTINCT pg_accessions) as protein_count,
                COUNT(DISTINCT peptidoform) as peptide_count,
                AVG(intensity) as avg_intensity,
                MIN(intensity) as min_intensity,
                MAX(intensity) as max_intensity,
                COUNT(DISTINCT reference_file_name) as sample_count
            FROM {source}
            GROUP BY reference_file_name, channel
            ORDER BY reference_file_name, channel
            """

        try:
            return self.query_to_df(sql)
        except Exception as e:
            self.logger.error(f"Error getting file statistics: {e}")
            return None

    def get_msstats_protein_file_matrix(
        self, protein_filter: Optional[str] = None
    ) -> Optional[pd.DataFrame]:
        """
        Get protein presence matrix across files from MSstats data using SQL.

        Args:
            protein_filter: Optional protein filter string

        Returns:
            DataFrame with protein presence matrix or None if no MSstats data
        """
        if not self._has_msstats_table():
            self.logger.warning("No MSstats data available")
            return None

        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_MSSTATS)

        where_clause = ""
        if protein_filter:
            where_clause = f"WHERE ProteinName LIKE '%{protein_filter}%'"

        sql = f"""
        SELECT 
            ProteinName as protein_accession,
            COUNT(DISTINCT Reference_Name) as file_count,
            string_agg(DISTINCT Reference_Name, ';') as files_present,
            COUNT(*) as total_features,
            AVG(Intensity) as avg_intensity
        FROM {source}
        {where_clause}
        GROUP BY ProteinName
        HAVING COUNT(DISTINCT Reference_Name) > 1
        ORDER BY file_count DESC, total_features DESC
        """

        try:
            return self.query_to_df(sql)
        except Exception as e:
            self.logger.error(f"Error getting protein file matrix: {e}")
            return None

    def iter_msstats_files(
        self, file_batch_size: int = 10
    ) -> Generator[pd.DataFrame, None, None]:
        """
        Iterate through MSstats data in file batches using SQL.

        Args:
            file_batch_size: Number of files to process in each batch

        Yields:
            DataFrame with MSstats data for each batch or None if no MSstats data
        """
        if not self._has_msstats_table():
            self.logger.warning("No MSstats data available")
            return None

        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_MSSTATS)

        # Get unique files
        files_sql = f"""SELECT {MSSTATS_REFERENCE_NAME},
            ROW_NUMBER() OVER (ORDER BY {MSSTATS_REFERENCE_NAME}) AS file_rank
        FROM {source}
        GROUP BY {MSSTATS_REFERENCE_NAME}
        """

        try:
            files_df = self.query_to_df(files_sql)
            total_files = len(files_df)

            for batch_start in range(0, total_files, file_batch_size):
                batch_end = min(batch_start + file_batch_size, total_files)

                # Get files for this batch
                batch_files = files_df.iloc[batch_start:batch_end][
                    MSSTATS_REFERENCE_NAME
                ].tolist()
                file_list_str = "', '".join(batch_files)

                # Use SQL to get batch data efficiently
                batch_sql = f"""
                SELECT *
                FROM {source} 
                WHERE {MSSTATS_REFERENCE_NAME} IN ('{file_list_str}')
                ORDER BY {MSSTATS_REFERENCE_NAME}, peptidoform
                """

                yield self.query_to_df(batch_sql)

        except Exception as e:
            self.logger.error(f"Error iterating MSstats files: {e}")
            return None

    def get_msstats_sample_mapping(self) -> Optional[Dict[str, str]]:
        """
        Get sample mapping from MSstats data if SDRF is available.

        Returns:
            Dictionary mapping file-channel combinations to sample accessions or None
        """
        if not self._sdrf_handler:
            self.logger.warning("No SDRF data available for sample mapping")
            return None

        try:
            return self._sdrf_handler.get_sample_map_run()
        except Exception as e:
            self.logger.error(f"Error getting sample mapping: {e}")
            return None

    def aggregate_msstats_channels(
        self, file_batch_size: int = 10
    ) -> Generator[pd.DataFrame, None, None]:
        """
        Aggregate channels per peptide for TMT/iTRAQ experiments using SQL.

        Args:
            file_batch_size: Number of files to process in each batch

        Yields:
            DataFrame with aggregated channel data or None if no MSstats data
        """
        if not self._has_msstats_table():
            self.logger.warning("No MSstats data available")
            return None

        experiment_type = self.get_msstats_experiment_type()
        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_MSSTATS)

        if experiment_type == "LFQ":
            # For LFQ, just return the processed data with intensities array
            sql = f"""
            SELECT 
                Reference_Name as reference_file_name,
                ProteinName as pg_accessions_raw,
                PeptideSequence as peptidoform,
                COALESCE(PrecursorCharge, Charge, 3) as precursor_charge,
                -- Extract protein accession
                CASE 
                    WHEN ProteinName LIKE '%|%|%' 
                    THEN split_part(ProteinName, '|', 2)
                    ELSE ProteinName
                END as anchor_protein,
                -- Check if unique protein
                CASE 
                    WHEN ProteinName LIKE '%;%' OR ProteinName LIKE '%,%' 
                    THEN 0 ELSE 1 
                END as unique_protein,
                -- Create intensities array (single element for LFQ)
                array_agg(STRUCT_PACK(
                    sample_accession := Reference,
                    channel := Channel,
                    intensity := Intensity
                )) as intensities
            FROM {source}
            GROUP BY Reference_Name, ProteinName, PeptideSequence, precursor_charge, 
                     anchor_protein, unique_protein
            ORDER BY Reference_Name, PeptideSequence
            """
        else:
            # For TMT/iTRAQ, aggregate channels using SQL
            sql = f"""
            WITH channel_aggregated AS (
                SELECT 
                    Reference_Name as reference_file_name,
                    ProteinName as pg_accessions_raw,
                    PeptideSequence as peptidoform,
                    COALESCE(PrecursorCharge, Charge, 3) as precursor_charge,
                    -- Extract protein accession
                    CASE 
                        WHEN ProteinName LIKE '%|%|%' 
                        THEN split_part(ProteinName, '|', 2)
                        ELSE ProteinName
                    END as anchor_protein,
                    -- Check if unique protein
                    CASE 
                        WHEN ProteinName LIKE '%;%' OR ProteinName LIKE '%,%' 
                        THEN 0 ELSE 1 
                    END as unique_protein,
                    -- Aggregate all channels for each peptide
                    array_agg(STRUCT_PACK(
                        sample_accession := Reference,
                        channel := Channel,
                        intensity := Intensity
                    )) as intensities
                FROM {source}
                GROUP BY Reference_Name, ProteinName, PeptideSequence, precursor_charge, 
                         anchor_protein, unique_protein
            )
            SELECT 
                reference_file_name,
                pg_accessions_raw,
                peptidoform, 
                precursor_charge,
                anchor_protein,
                unique_protein,
                intensities
            FROM channel_aggregated
            ORDER BY reference_file_name, peptidoform
            """

        try:
            return self.query_to_df(sql)
        except Exception as e:
            self.logger.error(f"Error aggregating MSstats channels: {e}")
            return None

    def get_msstats_peptide_summary(self) -> Optional[pd.DataFrame]:
        """
        Get peptide-level summary statistics from MSstats data.

        Returns:
            DataFrame with peptide summary or None if no MSstats data
        """
        if not self._has_msstats_table():
            self.logger.warning("No MSstats data available")
            return None

        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_MSSTATS)

        sql = f"""
        SELECT 
            peptidoform as peptide_sequence,
            COUNT(DISTINCT pg_accessions) as protein_count,
            COUNT(DISTINCT reference_file_name) as file_count,
            COUNT(*) as total_observations,
            AVG(intensity) as avg_intensity,
            MIN(intensity) as min_intensity,
            MAX(intensity) as max_intensity,
            STDDEV(intensity) as std_intensity,
            -- Check if peptide is unique to single protein
            CASE 
                WHEN COUNT(DISTINCT pg_accessions) = 1 THEN 1 
                ELSE 0 
            END as is_unique_peptide
        FROM {source}
        WHERE peptidoform IS NOT NULL
        GROUP BY peptidoform
        ORDER BY total_observations DESC
        """

        try:
            return self.query_to_df(sql)
        except Exception as e:
            self.logger.error(f"Error getting peptide summary: {e}")
            return None

    def get_msstats_protein_summary(self) -> Optional[pd.DataFrame]:
        """
        Get protein-level summary statistics from MSstats data.

        Returns:
            DataFrame with protein summary or None if no MSstats data
        """
        if not self._has_msstats_table():
            self.logger.warning("No MSstats data available")
            return None

        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_MSSTATS)

        sql = f"""
        SELECT 
            pg_accessions as protein_name,
            COUNT(DISTINCT peptidoform) as peptide_count,
            COUNT(DISTINCT reference_file_name) as file_count,
            COUNT(*) as total_observations,
            AVG(intensity) as avg_intensity,
            MIN(intensity) as min_intensity,
            MAX(intensity) as max_intensity,
            STDDEV(intensity) as std_intensity,
            -- Calculate coefficient of variation
            CASE 
                WHEN AVG(intensity) > 0 THEN STDDEV(intensity) / AVG(intensity)
                ELSE NULL
            END as cv_intensity
        FROM {source}
        WHERE pg_accessions IS NOT NULL
        GROUP BY pg_accessions
        ORDER BY total_observations DESC
        """

        try:
            return self.query_to_df(sql)
        except Exception as e:
            self.logger.error(f"Error getting protein summary: {e}")
            return None

    def get_msstats_intensity_distribution(
        self, log_transform: bool = True
    ) -> Optional[pd.DataFrame]:
        """
        Get intensity distribution statistics from MSstats data.

        Args:
            log_transform: Whether to log-transform intensities for analysis

        Returns:
            DataFrame with intensity distribution or None if no MSstats data
        """
        if not self._has_msstats_table():
            self.logger.warning("No MSstats data available")
            return None

        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_MSSTATS)

        if log_transform:
            intensity_expr = (
                "LOG10(CASE WHEN Intensity > 0 THEN Intensity ELSE NULL END)"
            )
            intensity_label = "log10_intensity"
        else:
            intensity_expr = "Intensity"
            intensity_label = "intensity"

        sql = f"""
        WITH intensity_stats AS (
            SELECT 
                {intensity_expr} as intensity_value
            FROM {source}
            WHERE Intensity IS NOT NULL AND Intensity > 0
        ),
        percentiles AS (
            SELECT 
                MIN(intensity_value) as min_intensity,
                MAX(intensity_value) as max_intensity,
                AVG(intensity_value) as mean_intensity,
                MEDIAN(intensity_value) as median_intensity,
                STDDEV(intensity_value) as std_intensity,
                percentile_cont(0.25) WITHIN GROUP (ORDER BY intensity_value) as q25,
                percentile_cont(0.75) WITHIN GROUP (ORDER BY intensity_value) as q75,
                percentile_cont(0.05) WITHIN GROUP (ORDER BY intensity_value) as q05,
                percentile_cont(0.95) WITHIN GROUP (ORDER BY intensity_value) as q95
            FROM intensity_stats
        )
        SELECT 
            '{intensity_label}' as measure_type,
            min_intensity,
            q05,
            q25,
            median_intensity,
            mean_intensity,
            q75,
            q95,
            max_intensity,
            std_intensity,
            (q75 - q25) as iqr
        FROM percentiles
        """

        try:
            return self.query_to_df(sql)
        except Exception as e:
            self.logger.error(f"Error getting intensity distribution: {e}")
            return None

    # ============================================================================
    # ENHANCED QUERY METHODS
    # ============================================================================

    def search_proteins_by_keyword(
        self, keyword: str, case_sensitive: bool = False
    ) -> pd.DataFrame:
        """
        Search proteins by keyword in accession or description.

        Args:
            keyword: Search keyword
            case_sensitive: Whether search should be case sensitive

        Returns:
            DataFrame with matching proteins
        """
        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_PROTEINS)
        if case_sensitive:
            query = f"""
            SELECT *
            FROM {source}
            WHERE accession LIKE '%{keyword}%' OR description LIKE '%{keyword}%'
            """
        else:
            query = f"""
            SELECT *
            FROM {source}
            WHERE LOWER(accession) LIKE LOWER('%{keyword}%') OR LOWER(description) LIKE LOWER('%{keyword}%')
            """
        return self.query_to_df(query)

    def get_proteins_by_score_threshold(
        self, score_column: str, threshold: float, operator: str = ">="
    ) -> pd.DataFrame:
        """
        Get proteins filtered by score threshold.

        Args:
            score_column: Column name for scoring
            threshold: Score threshold value
            operator: Comparison operator ('>', '>=', '<', '<=', '=')

        Returns:
            DataFrame with filtered proteins
        """
        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_PROTEINS)
        query = f"""
        SELECT *
        FROM {source}
        WHERE CAST({score_column} AS DOUBLE) {operator} {threshold}
        ORDER BY CAST({score_column} AS DOUBLE) DESC
        """
        return self.query_to_df(query)

    def get_psms_by_charge(self, charge: int) -> pd.DataFrame:
        """
        Get PSMs for a specific charge state.

        Args:
            charge: Charge state to filter by

        Returns:
            DataFrame with PSMs for the specified charge
        """
        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_PSMS)
        query = f"""
        SELECT *
        FROM {source}
        WHERE charge = {charge}
        """
        return self.query_to_df(query)

    def get_psms_by_mass_error(
        self, min_error: float, max_error: float
    ) -> pd.DataFrame:
        """
        Get PSMs within a mass error range.

        Args:
            min_error: Minimum mass error (ppm)
            max_error: Maximum mass error (ppm)

        Returns:
            DataFrame with PSMs in the mass error range
        """
        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_PSMS)
        query = f"""
        SELECT *
        FROM {source}
        WHERE CAST(mass_error AS DOUBLE) BETWEEN {min_error} AND {max_error}
        """
        return self.query_to_df(query)

    def get_proteins_by_coverage(
        self, min_coverage: float = 0.0, max_coverage: float = 100.0
    ) -> pd.DataFrame:
        """
        Get proteins filtered by coverage range.

        Args:
            min_coverage: Minimum protein coverage (%)
            max_coverage: Maximum protein coverage (%)

        Returns:
            DataFrame with proteins in the coverage range
        """
        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_PROTEINS)
        query = f"""
        SELECT *
        FROM {source}
        WHERE CAST(protein_coverage AS DOUBLE) BETWEEN {min_coverage} AND {max_coverage}
        ORDER BY CAST(protein_coverage AS DOUBLE) DESC
        """
        return self.query_to_df(query)

    def get_decoy_proteins(self) -> pd.DataFrame:
        """
        Get all decoy proteins.

        Returns:
            DataFrame with decoy proteins
        """
        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_PROTEINS)
        query = f"""
        SELECT *
        FROM {source}
        WHERE opt_global_cv_PRIDE:0000303_decoy_hit = 1
        """
        return self.query_to_df(query)

    def get_target_proteins(self) -> pd.DataFrame:
        """
        Get all target (non-decoy) proteins.

        Returns:
            DataFrame with target proteins
        """
        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_PROTEINS)
        query = f"""
        SELECT *
        FROM {source}
        WHERE opt_global_cv_PRIDE:0000303_decoy_hit = 0 OR opt_global_cv_PRIDE:0000303_decoy_hit IS NULL
        """
        return self.query_to_df(query)

    def export_summary_report(self, output_path: Union[Path, str]) -> Path:
        """
        Export a comprehensive summary report.

        Args:
            output_path: Path for the summary report

        Returns:
            Path to the created report
        """
        output_path = Path(output_path)

        # Gather statistics
        stats = {
            "Total Proteins": self.get_protein_count(),
            "Total Protein Details": self.get_protein_details_count(),
            "Total PSMs": self.get_psm_count(),
            "Target Proteins": len(self.get_target_proteins()),
            "Decoy Proteins": len(self.get_decoy_proteins()),
        }

        if self._has_msstats_table():
            msstats_df = self.get_msstats()
            if msstats_df is not None:
                stats["MSstats Entries"] = len(msstats_df)
                if "Reference" in msstats_df.columns:
                    stats["MSstats Runs"] = msstats_df["Reference"].nunique()

        # Get metadata summary
        metadata_df = self.get_metadata()
        if not metadata_df.empty:
            stats["Metadata Entries"] = len(metadata_df)

        # Write report
        with open(output_path, "w") as f:
            f.write("MzTabIndexer Summary Report\n")
            f.write("=" * 50 + "\n\n")

            for key, value in stats.items():
                f.write(f"{key}: {value:,}\n")

            f.write("\n" + "=" * 50 + "\n")
            f.write(f"Report generated: {pd.Timestamp.now()}\n")
            f.write(f"Backend: {self._backend}\n")
            if self._database_path:
                f.write(f"Database: {self._database_path}\n")

        self.logger.info(f"Summary report exported to {output_path}")
        return output_path

    # ============================================================================
    # PERFORMANCE OPTIMIZATION METHODS
    # ============================================================================

    def create_indices(self, indices_config: Optional[Dict[str, List[str]]] = None):
        """
        Create database indices for improved query performance.

        This method creates database indices on specified columns to accelerate
        common query patterns. Indices are only created for DuckDB backend.

        Args:
            indices_config: Dictionary mapping table names to lists of columns to index.
                If None, uses default indices for common query patterns:
                - proteins: accession, opt_global_cv_PRIDE:0000303_decoy_hit
                - psms: accession, charge, spectra_ref
                - msstats: ProteinName, Reference

        Note:
            - Only works with DuckDB backend. For parquet backend, indices are not applicable.
            - Indices are created with "IF NOT EXISTS" to avoid errors on repeated calls.
            - Column names are sanitized for index creation (special characters replaced).
            - Missing tables or columns are logged as warnings but don't stop processing.
        """
        if self._backend != "duckdb":
            self.logger.warning("Indices can only be created for DuckDB backend")
            return

        if not indices_config:
            # Default indices for common query patterns
            indices_config = {
                self._MZTAB_INDEXER_TABLE_PROTEINS: [
                    "accession",
                    "opt_global_cv_PRIDE:0000303_decoy_hit",
                ],
                self._MZTAB_INDEXER_TABLE_PSMS: ["accession", "charge", "spectra_ref"],
                self._MZTAB_INDEXER_TABLE_MSSTATS: ["ProteinName", "Reference"],
            }

        for table_name, columns in indices_config.items():
            # Check if table exists
            tables = self._duckdb.execute("SHOW TABLES").fetchall()
            table_names = [t[0] for t in tables]

            if table_name not in table_names:
                self.logger.warning(
                    f"Table {table_name} does not exist, skipping indices"
                )
                continue

            # Get actual columns in the table
            table_columns = [
                col[0]
                for col in self._duckdb.execute(
                    f"PRAGMA table_info({table_name})"
                ).fetchall()
            ]

            for column in columns:
                if column in table_columns:
                    safe_column_name = column.replace(".", "_").replace(":", "_")
                    index_name = f"idx_{table_name}_{safe_column_name}"

                    try:
                        self._duckdb.execute(
                            f'CREATE INDEX IF NOT EXISTS {index_name} ON {table_name} ("{column}")'
                        )
                        self.logger.info(
                            f"Created index {index_name} on {table_name}.{column}"
                        )
                    except Exception as e:
                        self.logger.warning(
                            f"Failed to create index on {table_name}.{column}: {e}"
                        )
                else:
                    self.logger.warning(
                        f"Column {column} not found in table {table_name}"
                    )

    def analyze_tables(self):
        """
        Run ANALYZE on all tables to update statistics for query optimization.

        This method checks the backend type and calls the appropriate parent method.
        For DuckDB backend, it updates table statistics for query optimization.
        For Parquet backend, it logs a warning since ANALYZE is not applicable.
        """
        if self._backend == "duckdb":
            super().analyze_tables()
        else:
            self.logger.warning(
                "ANALYZE can only be run for DuckDB backend, not parquet backend"
            )

    def vacuum_database(self):
        """
        Run VACUUM to reclaim storage and optimize the database.

        This method checks the backend type and calls the appropriate parent method.
        For DuckDB backend, it reclaims storage and optimizes the database structure.
        For Parquet backend, it logs a warning since VACUUM is not applicable.
        """
        if self._backend == "duckdb":
            super().vacuum_database()
        else:
            self.logger.warning(
                "VACUUM can only be run for DuckDB backend, not parquet backend"
            )

    def get_data_quality_metrics(self) -> Dict[str, Any]:
        """
        Calculate comprehensive data quality metrics.

        Returns:
            Dictionary with various quality metrics
        """
        metrics = {}

        # Protein metrics
        proteins_df = self.get_proteins()
        if not proteins_df.empty:
            metrics["proteins"] = {
                "total_count": len(proteins_df),
                "target_count": len(self.get_target_proteins()),
                "decoy_count": len(self.get_decoy_proteins()),
                "unique_accessions": proteins_df["accession"].nunique(),
                "null_accessions": proteins_df["accession"].isnull().sum(),
                "null_descriptions": (
                    proteins_df["description"].isnull().sum()
                    if "description" in proteins_df.columns
                    else 0
                ),
            }

            # Coverage statistics if available
            if "protein_coverage" in proteins_df.columns:
                coverage_values = pd.to_numeric(
                    proteins_df["protein_coverage"], errors="coerce"
                )
                metrics["proteins"]["coverage_stats"] = {
                    "mean": coverage_values.mean(),
                    "median": coverage_values.median(),
                    "min": coverage_values.min(),
                    "max": coverage_values.max(),
                    "std": coverage_values.std(),
                }

        # PSM metrics
        psms_df = self.get_psms()
        if not psms_df.empty:
            metrics["psms"] = {
                "total_count": len(psms_df),
                "unique_accessions": psms_df["accession"].nunique(),
                "null_accessions": psms_df["accession"].isnull().sum(),
            }

            # Charge distribution if available
            if "charge" in psms_df.columns:
                charge_counts = psms_df["charge"].value_counts().to_dict()
                metrics["psms"]["charge_distribution"] = charge_counts

        # MSstats metrics if available
        if self._has_msstats_table():
            msstats_df = self.get_msstats()
            if msstats_df is not None and not msstats_df.empty:
                metrics["msstats"] = {
                    "total_count": len(msstats_df),
                    "unique_proteins": (
                        msstats_df["ProteinName"].nunique()
                        if "ProteinName" in msstats_df.columns
                        else 0
                    ),
                    "unique_runs": (
                        msstats_df["Reference"].nunique()
                        if "Reference" in msstats_df.columns
                        else 0
                    ),
                }

                # Intensity statistics if available
                if "Intensity" in msstats_df.columns:
                    intensity_values = pd.to_numeric(
                        msstats_df["Intensity"], errors="coerce"
                    )
                    metrics["msstats"]["intensity_stats"] = {
                        "mean": intensity_values.mean(),
                        "median": intensity_values.median(),
                        "min": intensity_values.min(),
                        "max": intensity_values.max(),
                        "std": intensity_values.std(),
                    }

        return metrics

    # ============================================================================
    # STATISTICAL ANALYSIS METHODS
    # ============================================================================

    def calculate_fdr(
        self,
        score_column: str = "best_search_engine_score[1]",
        target_column: str = "opt_global_cv_PRIDE:0000303_decoy_hit",
    ) -> pd.DataFrame:
        """
        Calculate False Discovery Rate (FDR) for proteins.

        Args:
            score_column: Column name for scoring
            target_column: Column name indicating target/decoy status

        Returns:
            DataFrame with FDR calculations
        """
        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_PROTEINS)

        query = f"""
        WITH ranked_proteins AS (
            SELECT 
                accession,
                {score_column} as score,
                {target_column} as is_decoy,
                ROW_NUMBER() OVER (ORDER BY CAST({score_column} AS DOUBLE) DESC) as rank
            FROM {source}
            WHERE {score_column} IS NOT NULL AND {score_column} != 'null'
        ),
        cumulative_counts AS (
            SELECT 
                rank,
                score,
                is_decoy,
                SUM(CASE WHEN is_decoy = 1 THEN 1 ELSE 0 END) OVER (ORDER BY rank) as cumulative_decoys,
                SUM(CASE WHEN is_decoy = 0 THEN 1 ELSE 0 END) OVER (ORDER BY rank) as cumulative_targets
            FROM ranked_proteins
        )
        SELECT 
            rank,
            score,
            is_decoy,
            cumulative_decoys,
            cumulative_targets,
            CASE 
                WHEN cumulative_targets > 0 THEN 
                    CAST(cumulative_decoys AS DOUBLE) / CAST(cumulative_targets AS DOUBLE)
                ELSE 0 
            END as fdr
        FROM cumulative_counts
        ORDER BY rank
        """

        return self.query_to_df(query)

    def get_protein_coverage_distribution(self) -> pd.DataFrame:
        """
        Get distribution of protein coverage values.

        Returns:
            DataFrame with coverage distribution statistics
        """
        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_PROTEINS)

        query = f"""
        SELECT 
            CASE 
                WHEN CAST(protein_coverage AS DOUBLE) < 10 THEN '0-10%'
                WHEN CAST(protein_coverage AS DOUBLE) < 25 THEN '10-25%'
                WHEN CAST(protein_coverage AS DOUBLE) < 50 THEN '25-50%'
                WHEN CAST(protein_coverage AS DOUBLE) < 75 THEN '50-75%'
                WHEN CAST(protein_coverage AS DOUBLE) < 100 THEN '75-100%'
                ELSE '100%'
            END as coverage_range,
            COUNT(*) as protein_count,
            AVG(CAST(protein_coverage AS DOUBLE)) as avg_coverage
        FROM {source}
        WHERE protein_coverage IS NOT NULL AND protein_coverage != 'null'
        GROUP BY coverage_range
        ORDER BY 
            CASE coverage_range
                WHEN '0-10%' THEN 1
                WHEN '10-25%' THEN 2
                WHEN '25-50%' THEN 3
                WHEN '50-75%' THEN 4
                WHEN '75-100%' THEN 5
                WHEN '100%' THEN 6
            END
        """

        return self.query_to_df(query)

    def get_charge_state_distribution(self) -> pd.DataFrame:
        """
        Get distribution of PSM charge states.

        Returns:
            DataFrame with charge state distribution
        """
        source = self._get_table_source(self._MZTAB_INDEXER_TABLE_PSMS)

        query = f"""
        SELECT 
            charge,
            COUNT(*) as psm_count,
            COUNT(*) * 100.0 / (SELECT COUNT(*) FROM {source}) as percentage
        FROM {source}
        WHERE charge IS NOT NULL
        GROUP BY charge
        ORDER BY charge
        """

        return self.query_to_df(query)

    def get_protein_identification_statistics(self) -> Dict[str, Any]:
        """
        Get comprehensive protein identification statistics.

        Returns:
            Dictionary with various identification statistics
        """
        stats = {}

        # Basic counts
        stats["total_proteins"] = self.get_protein_count()
        stats["target_proteins"] = len(self.get_target_proteins())
        stats["decoy_proteins"] = len(self.get_decoy_proteins())

        # Calculate FDR
        fdr_df = self.calculate_fdr()
        if not fdr_df.empty:
            # Find FDR at different thresholds
            for fdr_threshold in [0.01, 0.05, 0.10]:
                threshold_data = fdr_df[fdr_df["fdr"] <= fdr_threshold]
                if not threshold_data.empty:
                    stats[f"proteins_at_fdr_{int(fdr_threshold*100)}"] = len(
                        threshold_data
                    )
                else:
                    stats[f"proteins_at_fdr_{int(fdr_threshold*100)}"] = 0

        # Coverage statistics
        coverage_df = self.get_protein_coverage_distribution()
        if not coverage_df.empty:
            stats["coverage_distribution"] = coverage_df.to_dict("records")

        # Charge state statistics
        charge_df = self.get_charge_state_distribution()
        if not charge_df.empty:
            stats["charge_distribution"] = charge_df.to_dict("records")

        return stats

    # ============================================================================
    # INTEGRATION AND UTILITY METHODS
    # ============================================================================

    def get_database_info(self) -> Dict[str, Any]:
        """
        Get comprehensive information about the database.

        Returns:
            Dictionary with database information
        """
        info = {
            "backend": self._backend,
            "database_path": self._database_path,
            "mztab_path": self._mztab_path,
            "msstats_path": self._msstats_path,
            "batch_size": self._batch_size,
            "worker_threads": self._worker_threads,
            "max_memory": self._max_memory,
        }

        # Get table information
        if self._backend == "duckdb":
            tables = self._duckdb.execute("SHOW TABLES").fetchall()
            table_info = {}
            for table in tables:
                table_name = table[0]
                try:
                    count = self._duckdb.execute(
                        f"SELECT COUNT(*) FROM {table_name}"
                    ).fetchone()[0]
                    table_info[table_name] = {"row_count": count}
                except Exception as e:
                    table_info[table_name] = {"error": str(e)}
            info["tables"] = table_info
        else:
            # For parquet backend, check files
            if self._database_path and Path(self._database_path).exists():
                parquet_files = list(Path(self._database_path).glob("*.parquet"))
                info["parquet_files"] = [f.name for f in parquet_files]

        return info

    def backup_database(self, backup_path: Union[Path, str]) -> Path:
        """
        Create a backup of the database.

        Args:
            backup_path: Path for the backup file

        Returns:
            Path to the created backup
        """
        backup_path = Path(backup_path)

        if self._backend == "duckdb":
            # For DuckDB, copy the database file
            import shutil

            shutil.copy2(self._database_path, backup_path)
            self.logger.info(f"Database backed up to {backup_path}")
        else:
            # For parquet backend, create a zip archive
            import zipfile

            with zipfile.ZipFile(backup_path, "w", zipfile.ZIP_DEFLATED) as zipf:
                for parquet_file in Path(self._database_path).glob("*.parquet"):
                    zipf.write(parquet_file, parquet_file.name)
            self.logger.info(f"Parquet files backed up to {backup_path}")

        return backup_path

    def get_memory_usage(self) -> Dict[str, Any]:
        """
        Get memory usage statistics for the database.

        Returns:
            Dictionary with memory usage information
        """
        if self._backend != "duckdb":
            return {
                "backend": self._backend,
                "note": "Memory usage only available for DuckDB backend",
            }

        try:
            # Get DuckDB memory usage
            memory_stats = self._duckdb.execute("PRAGMA memory_usage").fetchall()
            memory_dict = {row[0]: row[1] for row in memory_stats}

            # Get table sizes
            tables = self._duckdb.execute("SHOW TABLES").fetchall()
            table_sizes = {}
            for table in tables:
                table_name = table[0]
                try:
                    size_query = f"SELECT SUM(column_size) FROM pragma_table_info('{table_name}')"
                    size = self._duckdb.execute(size_query).fetchone()[0]
                    table_sizes[table_name] = size
                except Exception:
                    table_sizes[table_name] = "unknown"

            return {
                "backend": self._backend,
                "memory_usage": memory_dict,
                "table_sizes": table_sizes,
            }
        except Exception as e:
            return {"backend": self._backend, "error": str(e)}

    def optimize_database(self):
        """
        Perform comprehensive database optimization.

        This method runs a series of optimization operations to improve
        database performance and storage efficiency:

        1. Create indices on commonly queried columns
        2. Analyze tables to update statistics for query optimization
        3. Vacuum database to reclaim storage and optimize structure

        Note:
            - Only applies to DuckDB backend. For parquet backend, only
              index creation is attempted (and will be skipped).
            - This method can be called multiple times safely.
            - Optimization operations are logged for monitoring.
        """
        self.logger.info("Starting database optimization...")

        # Create indices
        self.create_indices()

        # Analyze tables
        self.analyze_tables()

        # Vacuum database
        self.vacuum_database()

        self.logger.info("Database optimization completed")

    def get_unique_from_psm_table(self):
        database_query = """
                SELECT COUNT(DISTINCT database) AS database_count
                FROM psms;
                """

        if self._duckdb.execute(database_query).fetchone()[0] > 1:
            self.logger.warning(
                "Cannot calculate unique values when multiple databases are present in the PSM table."
            )

        else:
            unique_query = """
                SELECT DISTINCT
                    accession,
                    "opt_global_cv_MS:1000889_peptidoform_sequence" AS peptidoform,
                    "unique"
                FROM psms;
            """

            unique_peptide = self._duckdb.execute(unique_query).df()

            unique_peptide["pg_accessions"] = unique_peptide["accession"].apply(
                dedup_accession
            )
            unique_peptide = (
                unique_peptide[["pg_accessions", "peptidoform", "unique"]]
                .drop_duplicates()
                .reset_index(drop=True)
            )

            return unique_peptide


def dedup_accession(accession_str):
    if pd.isna(accession_str):
        return accession_str
    seen = set()
    result = []
    for item in re.split(r"[;,]", accession_str):
        item = item.strip()
        if item not in seen:
            seen.add(item)
            result.append(item)

    return ";".join(result)
