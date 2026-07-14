"""
crispr_designer/systems.py
~~~~~~~~~~~~~~~~~~~~~~~~~~
Implements concrete CRISPR endonuclease systems:
SpCas9, SaCas9, Cas12a, and Prime Editor.
Utilizes class inheritance and defines custom features like Azimuth scoring
and sticky-end cleavage metadata.
"""

from typing import Dict, Any, Optional
from .base import CRISPRSystem
from .utils import gc_content

# Global CRISPR system registry
CRISPR_REGISTRY: Dict[str, CRISPRSystem] = {}

def register_system(name: str):
    """
    Decorator to automatically register a CRISPR system class in the global registry.
    """
    def decorator(cls):
        # We instantiate a default singleton of the system for the registry
        CRISPR_REGISTRY[name] = cls()
        return cls
    return decorator


@register_system('SpCas9')
class SpCas9System(CRISPRSystem):
    """
    Streptococcus pyogenes Cas9 (SpCas9) system.
    PAM: NGG (3' prime)
    Spacer Length: 20 nt
    Scoring: Azimuth (Bioinformatics predictive algorithm for guide efficiency)
    """

    def __init__(self):
        super().__init__(
            name='SpCas9',
            pam='NGG',
            spacer_length=20,
            pam_position='3_prime',
            scoring_method='Azimuth'
        )

    def get_context_window(self) -> Optional[tuple]:
        # Azimuth scoring uses a 30nt window: 4nt upstream + 20nt spacer + 3nt PAM + 3nt downstream
        return (4, 3)

    def score_candidate(self, spacer: str, target_site: str) -> float:
        """
        Calculates a mock 'Azimuth' on-target activity score based on biological properties:
        - GC content (optimal 40% - 60%).
        - U6 promoter compatibility (starts with a 'G' or 'A' at position 1).
        - Absence of poly-T terminators ('TTTT') which aborts transcription in standard U6/U3 promoters.
        - Nucleotide preference at critical spacer positions (e.g. high preference for G/C near PAM).
        """
        if not spacer or len(spacer) < 20:
            return 0.0

        spacer = spacer.upper()
        score = 60.0  # Base line score

        # 1. GC Content Penalty (ideal GC is 40% - 60% or 0.4 - 0.6)
        gc_ratio = gc_content(spacer)
        if 0.4 <= gc_ratio <= 0.6:
            score += 15.0
        elif gc_ratio < 0.3 or gc_ratio > 0.7:
            score -= 20.0

        # 2. Transcription Terminator Penalty
        if 'TTTT' in spacer:
            score -= 30.0  # Severe penalty: standard U6 RNA Polymerase III will terminate early

        # 3. U6/U3 Promoter compatibility (Position 1 at 5' end should ideally be 'G')
        if spacer[0] == 'G':
            score += 10.0
        elif spacer[0] == 'A':
            score += 5.0

        # 4. Position-specific preferences (nucleotides close to PAM / seed region: positions 10-20)
        # Seed region mutations are highly sensitive. Pyrimidines (C/T) at the position right next to PAM can lower activity.
        seed_region = spacer[-5:]  # last 5 bases before PAM
        gc_seed = sum(1 for b in seed_region if b in ('G', 'C'))
        score += (gc_seed * 3.0)  # Reward GC stability in the seed region

        # Clamp between 0.0 and 100.0
        return max(0.0, min(100.0, score))


@register_system('SaCas9')
class SaCas9System(CRISPRSystem):
    """
    Staphylococcus aureus Cas9 (SaCas9) system.
    PAM: NNGRRT (3' prime)
    Spacer Length: 21 nt
    Scoring: Standard
    """

    def __init__(self):
        super().__init__(
            name='SaCas9',
            pam='NNGRRT',
            spacer_length=21,
            pam_position='3_prime',
            scoring_method='Standard'
        )


@register_system('Cas12a')
class Cas12aSystem(CRISPRSystem):
    """
    Acidaminococcus / Lachnospiraceae Cas12a (Cpf1) system.
    PAM: TTTN (5' prime)
    Spacer Length: 23 nt
    Cleavage: Sticky end cut leaving a 5-nucleotide 5' overhang
    """

    def __init__(self):
        super().__init__(
            name='Cas12a',
            pam='TTTN',
            spacer_length=23,
            pam_position='5_prime',
            scoring_method='Standard',
            cut_type='Sticky end cut',
            overhang_len=5,
            overhang_direction="5_prime"
        )

    def score_candidate(self, spacer: str, target_site: str) -> float:
        """
        Cas12a has different base preferences compared to Cas9:
        - Prefers lower GC content (30% - 50% is optimal for Cas12a).
        - Avoids poly-T sequence (TTTT).
        """
        if not spacer:
            return 0.0
        
        spacer = spacer.upper()
        score = 55.0

        # GC ratio (Cas12a prefers slightly AT-richer targets)
        gc_ratio = gc_content(spacer)
        if 0.3 <= gc_ratio <= 0.5:
            score += 20.0
        elif gc_ratio > 0.6:
            score -= 15.0

        if 'TTTT' in spacer:
            score -= 25.0

        # Seed region for Cas12a is at the 5' end of the spacer (positions 1-8 downstream of 5' PAM)
        seed_region = spacer[:8]
        at_seed = sum(1 for b in seed_region if b in ('A', 'T'))
        score += (at_seed * 2.0)  # Reward AT rich seeds slightly for open strand kinetics

        return max(0.0, min(100.0, score))


@register_system('Prime_Editor')
class PrimeEditorSystem(CRISPRSystem):
    """
    Prime Editor system.
    Based on SpCas9 H840A nickase fused to reverse transcriptase.
    PAM: NGG (3' prime)
    Spacer Length: 20 nt
    Additional parameters: Primer Binding Site (PBS) length, Reverse Transcription (RT) template length.
    """

    def __init__(self, pbs_length: int = 13, rt_length: int = 15):
        # We can dynamically change these lengths or pass them on instantiation
        self.pbs_length = pbs_length
        self.rt_length = rt_length
        
        super().__init__(
            name='Prime_Editor',
            pam='NGG',
            spacer_length=20,
            pam_position='3_prime',
            scoring_method='Standard',
            pbs_len=pbs_length,
            rt_len=rt_length,
            enzyme_variant='PE2/PE3 (nCas9-M-MLV RT)'
        )

    def get_context_window(self) -> Optional[tuple]:
        # Same 30nt SpCas9-derived scaffold context window (PE2/PE3 nickase retains the NGG PAM geometry)
        return (4, 3)

    def update_prime_parameters(self, pbs_length: int, rt_length: int):
        """
        Allows dynamic updating of the Prime Editing pegRNA parameters.
        """
        self.pbs_length = pbs_length
        self.rt_length = rt_length
        self.extra_metadata['pbs_len'] = pbs_length
        self.extra_metadata['rt_len'] = rt_length


@register_system('Cas14a(Cas12f)')
class Cas14aSystem(CRISPRSystem):
    """
    Cas14a (also known as miniature Cas12f) system.
    PAM: TTTA (5' prime)
    Spacer Length: 20 nt
    Target: ssDNA (single-stranded DNA)
    """

    def __init__(self):
        super().__init__(
            name='Cas14a(Cas12f)',
            pam='TTTA',
            spacer_length=20,
            pam_position='5_prime',
            scoring_method='Standard',
            target_type='ssDNA'
        )

