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
    Cas13dSystem,
    Cas14aSystem
)
from .organisms import (
    ORGANISM_REGISTRY,
    fetch_genomic_sequence,
    get_organism_metadata
)
from .designer import CRISPRDesigner
from .predictors import PredictorFactory, BasePredictor

__all__ = [
    "CRISPRDesigner",
    "CRISPRSystem",
    "CRISPR_REGISTRY",
    "ORGANISM_REGISTRY",
    "register_system",
    "fetch_genomic_sequence",
    "get_organism_metadata",
    "SpCas9System",
    "SaCas9System",
    "Cas12aSystem",
    "PrimeEditorSystem",
    "Cas13dSystem",
    "Cas14aSystem",
    "PredictorFactory",
    "BasePredictor"
]
