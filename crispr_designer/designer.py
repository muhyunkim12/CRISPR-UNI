"""
crispr_designer/designer.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Implements the main orchestrator class CRISPRDesigner.
Supports dynamic switching of genomes and CRISPR endonuclease configurations,
and exposes a generic interface for guide RNA candidate scanning.
"""

from typing import List, Dict, Any, Optional
from .systems import CRISPR_REGISTRY, CRISPRSystem
from .organisms import ORGANISM_REGISTRY, fetch_genomic_sequence

class CRISPRDesigner:
    """
    Main user interface and orchestration controller for the Universal CRISPR Designer.
    Allows runtime configuration of organisms and CRISPR systems, and delegates
    scanning tasks to specific systems.
    """

    def __init__(
        self, 
        organism_name: Optional[str] = None, 
        crispr_system_name: Optional[str] = None
    ):
        """
        Initializes the CRISPRDesigner.

        Args:
            organism_name (str, optional): Target organism name (e.g. '벼', 'Human').
            crispr_system_name (str, optional): CRISPR endonuclease name (e.g. 'SpCas9').
        """
        self.organism_name: Optional[str] = None
        self.organism_metadata: Optional[Dict[str, Any]] = None
        self.active_system: Optional[CRISPRSystem] = None

        if organism_name:
            self.set_organism(organism_name)
        if crispr_system_name:
            self.set_crispr_system(crispr_system_name)

    def set_organism(self, name: str) -> 'CRISPRDesigner':
        """
        Dynamically changes the target organism.

        Args:
            name (str): Organism name in Korean or English (e.g. '벼', 'Rice', 'Human').

        Returns:
            CRISPRDesigner: Self, allowing method chaining.
        """
        if name not in ORGANISM_REGISTRY:
            raise KeyError(
                f"Organism '{name}' is not registered. "
                f"Available organisms: {list(ORGANISM_REGISTRY.keys())}"
            )
        self.organism_name = name
        self.organism_metadata = ORGANISM_REGISTRY[name]
        return self

    def set_crispr_system(self, system_name: str) -> 'CRISPRDesigner':
        """
        Dynamically switches the active CRISPR system configuration.

        Args:
            system_name (str): System identifier (e.g., 'SpCas9', 'SaCas9', 'Cas12a', 'Prime_Editor').

        Returns:
            CRISPRDesigner: Self, allowing method chaining.
        """
        if system_name not in CRISPR_REGISTRY:
            raise KeyError(
                f"CRISPR system '{system_name}' is not registered. "
                f"Available CRISPR systems: {list(CRISPR_REGISTRY.keys())}"
            )
        self.active_system = CRISPR_REGISTRY[system_name]
        return self

    def register_custom_system(self, system_name: str, system_instance: CRISPRSystem) -> 'CRISPRDesigner':
        """
        Allows dynamic registration of a custom CRISPR system at runtime
        without modifying core source code files.

        Args:
            system_name (str): Unique name of the new system.
            system_instance (CRISPRSystem): An instantiated class inheriting from CRISPRSystem.
        """
        if not isinstance(system_instance, CRISPRSystem):
            raise TypeError("The custom system must be an instance of CRISPRSystem.")
        
        CRISPR_REGISTRY[system_name] = system_instance
        # Optionally make it active immediately
        self.active_system = system_instance
        return self

    def fetch_reference_sequence(self, chromosome: str, start: Optional[int] = None, end: Optional[int] = None, force_mock: bool = True) -> str:
        """
        Helper method to fetch chromosome sequences for the currently set organism.
        """
        if not self.organism_name:
            raise ValueError("No organism configured. Please call set_organism() first.")
        
        return fetch_genomic_sequence(
            organism_name=self.organism_name,
            chromosome=chromosome,
            start=start,
            end=end,
            force_mock=force_mock
        )

    def find_candidates(self, sequence: str) -> List[Dict[str, Any]]:
        """
        Finds all target guide RNA candidates in the given sequence
        using the active CRISPR endonuclease configuration.

        Args:
            sequence (str): Target nucleotide sequence.

        Returns:
            List[Dict[str, Any]]: Guide RNA candidates metadata dictionary list.
        """
        if not self.active_system:
            raise ValueError(
                "No CRISPR system configured. "
                "Please call set_crispr_system() before scanning."
            )
        
        # Delegate sequence scanning to the modular active CRISPR system class
        return self.active_system.find_candidates(sequence)
