"""
CRISPR Designer Platform
~~~~~~~~~~~~~~~~~~~~~~~~
An extensible, modular platform for designing guide RNA candidates across 
multiple organisms and CRISPR endonuclease/editing configurations.
"""

__version__ = "1.0.0"

from .base import CRISPRSystem
from .systems import (
    CRISPR_REGISTRY,
    register_system,
    SpCas9System,
    SaCas9System,
    Cas12aSystem,
    PrimeEditorSystem,
    Cas14aSystem
)
from .organisms import (
    ORGANISM_REGISTRY,
    fetch_genomic_sequence,
    get_organism_metadata,
    get_reference_genome_accession,
    get_local_genome_dir,
    get_local_genome_path,
    REFERENCE_GENOME_ACCESSIONS
)
from .designer import CRISPRDesigner
from .predictors import PredictorFactory, BasePredictor
from .offtarget import (
    search_offtargets,
    offtarget_risk_score,
    calc_cfd_score,
    calc_sacas9_score,
    is_cas_offinder_available,
    is_reference_genome_ready,
    CasOffinderNotFoundError,
    ReferenceGenomeMissingError,
)

__all__ = [
    "CRISPRDesigner",
    "CRISPRSystem",
    "CRISPR_REGISTRY",
    "ORGANISM_REGISTRY",
    "register_system",
    "fetch_genomic_sequence",
    "get_organism_metadata",
    "get_reference_genome_accession",
    "get_local_genome_dir",
    "get_local_genome_path",
    "REFERENCE_GENOME_ACCESSIONS",
    "search_offtargets",
    "offtarget_risk_score",
    "calc_cfd_score",
    "calc_sacas9_score",
    "is_cas_offinder_available",
    "is_reference_genome_ready",
    "CasOffinderNotFoundError",
    "ReferenceGenomeMissingError",
    "SpCas9System",
    "SaCas9System",
    "Cas12aSystem",
    "PrimeEditorSystem",
    "Cas14aSystem",
    "PredictorFactory",
    "BasePredictor"
]
