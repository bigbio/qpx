"""QPX converter layer -- tool-specific adapters that produce QPX Parquet files.

Each converter:
    1. Creates its own DuckDB connection (NOT shared with Dataset)
    2. Loads tool output into DuckDB
    3. Transforms via SQL into QPX schema
    4. Pipes results into a Writer (FeatureWriter, PsmWriter, PgWriter, etc.)
"""

from qpx.converters.base import BaseConverter
from qpx.converters.cdap.converter import CdapConverter
from qpx.converters.cdap.feature_adapter import CdapFeatureAdapter
from qpx.converters.cdap.pg_adapter import CdapPgAdapter
from qpx.converters.cdap.psm_adapter import CdapPsmAdapter
from qpx.converters.diann.converter import DiaNNConverter
from qpx.converters.diann.feature_adapter import DiannFeatureAdapter
from qpx.converters.diann.pg_adapter import DiannPgAdapter
from qpx.converters.fragpipe.converter import FragPipeConverter
from qpx.converters.fragpipe.feature_adapter import FragPipeFeatureAdapter
from qpx.converters.fragpipe.pg_adapter import FragPipePgAdapter
from qpx.converters.fragpipe.psm_adapter import FragPipePsmAdapter
from qpx.converters.maxquant.converter import MaxQuantConverter
from qpx.converters.maxquant.feature_adapter import MaxQuantFeatureAdapter
from qpx.converters.maxquant.pg_adapter import MaxQuantPgAdapter
from qpx.converters.maxquant.psm_adapter import MaxQuantPsmAdapter
from qpx.converters.mzidentml.psm_adapter import MzIdentMLPsmAdapter
from qpx.converters.openms.converter import OpenMSConverter
from qpx.converters.orchestrator import BaseOrchestrator, build_dataset_record
from qpx.converters.quantms_msstats.converter import QuantmsMsstatsConverter
from qpx.converters.quantms_msstats.feature_adapter import QuantmsMsstatsFeatureAdapter
from qpx.converters.sdrf import SdrfConverter
from qpx.converters.spectronaut.converter import SpectronautConverter
from qpx.converters.spectronaut.feature_adapter import SpectronautFeatureAdapter
from qpx.converters.spectronaut.pg_adapter import SpectronautPgAdapter

__all__ = [
    # Base
    "BaseConverter",
    "BaseOrchestrator",
    "build_dataset_record",
    # SDRF
    "SdrfConverter",
    # CDAP (CPTAC .psm)
    "CdapPsmAdapter",
    "CdapFeatureAdapter",
    "CdapPgAdapter",
    "CdapConverter",
    # DIA-NN
    "DiannFeatureAdapter",
    "DiannPgAdapter",
    "DiaNNConverter",
    # MaxQuant
    "MaxQuantPsmAdapter",
    "MaxQuantFeatureAdapter",
    "MaxQuantPgAdapter",
    "MaxQuantConverter",
    # FragPipe
    "FragPipePsmAdapter",
    "FragPipeFeatureAdapter",
    "FragPipePgAdapter",
    "FragPipeConverter",
    # mzIdentML
    "MzIdentMLPsmAdapter",
    # OpenMS
    "OpenMSConverter",
    # QuantMS MSstats
    "QuantmsMsstatsFeatureAdapter",
    "QuantmsMsstatsConverter",
    # Spectronaut
    "SpectronautFeatureAdapter",
    "SpectronautPgAdapter",
    "SpectronautConverter",
]
