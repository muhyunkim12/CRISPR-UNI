"""
crispr_designer/base.py
~~~~~~~~~~~~~~~~~~~~~~~
Defines the base abstract architecture class for CRISPR systems.
Handles dual-strand scanning, overlapping PAM search via regex lookaheads,
and genomic coordinate mapping.
"""

import re
from abc import ABC, abstractmethod
from typing import Dict, List, Any, Optional

from .utils import reverse_complement, iupac_to_regex

class CRISPRSystem(ABC):
    """
    Abstract Base Class representing a generic CRISPR-Cas system.
    Any new CRISPR endonuclease system (e.g. Cas9, Cas12a) inherits from this.
    """

    def __init__(
        self,
        name: str,
        pam: str,
        spacer_length: int,
        pam_position: str = '3_prime',  # '3_prime' (Cas9) or '5_prime' (Cas12a)
        scoring_method: str = 'Standard',
        target_type: str = 'DNA',
        **kwargs
    ):
        """
        Initializes the CRISPR System metadata.

        Args:
            name (str): Unique name of the CRISPR system (e.g., 'SpCas9').
            pam (str): IUPAC PAM pattern (e.g., 'NGG', 'NNGRRT', 'TTTN').
            spacer_length (int): Target spacer sequence length in nucleotides.
            pam_position (str): Position of PAM relative to the spacer ('3_prime' or '5_prime').
            scoring_method (str): Identifier for the scoring algorithm.
            target_type (str): Substrate target molecule ('DNA', 'ssDNA', 'RNA').
            **kwargs: Extra parameters for system-specific metadata.
        """
        self.name = name
        self.pam = pam
        self.spacer_length = spacer_length
        self.target_type = target_type
        
        if pam_position not in ('3_prime', '5_prime'):
            raise ValueError("pam_position must be either '3_prime' or '5_prime'.")
        self.pam_position = pam_position
        self.scoring_method = scoring_method
        self.extra_metadata = kwargs

        # Pre-compile the IUPAC PAM pattern to a regular expression for fast scanning
        self.pam_regex_str = iupac_to_regex(self.pam)
        # Lookahead assertion regex to allow matching of overlapping PAM sequences
        self.scanner_regex = re.compile(f"(?=({self.pam_regex_str}))", re.IGNORECASE)

    def find_candidates(self, sequence: str) -> List[Dict[str, Any]]:
        """
        Scans a given sequence for target guide RNA candidates.
        Supports double-strand scanning for double-stranded DNA targets, and
        single-strand bypass scanning for single-stranded substrates (RNA, ssDNA).

        Args:
            sequence (str): Target sequence to search (DNA or RNA).

        Returns:
            List[Dict[str, Any]]: A list of guide RNA candidates.
        """
        candidates = []
        seq_len = len(sequence)

        # Normalize RNA sequence characters (U -> T) to enable standard DNA IUPAC regex scanning
        normalized_seq = sequence.replace('U', 'T').replace('u', 't')

        # 1. Scan Forward Strand (+)
        forward_candidates = self._scan_strand(normalized_seq, strand='+', seq_len=seq_len, original_seq=sequence)
        candidates.extend(forward_candidates)

        # 2. Scan Reverse Complement Strand (-) - Bypassed for single-stranded targets (RNA, ssDNA)
        if self.target_type == 'DNA':
            rev_sequence = reverse_complement(normalized_seq)
            reverse_candidates = self._scan_strand(rev_sequence, strand='-', seq_len=seq_len, original_seq=sequence)
            candidates.extend(reverse_candidates)

        return candidates

    def _scan_strand(self, target_seq: str, strand: str, seq_len: int, original_seq: str) -> List[Dict[str, Any]]:
        """
        Internal scanning helper for a single sequence strand.
        """
        strand_candidates = []

        # Find all overlapping PAM matches using lookahead regex
        for match in self.scanner_regex.finditer(target_seq):
            pam_match_normalized = match.group(1)
            pam_len = len(pam_match_normalized)
            match_start = match.start()
            match_end = match_start + pam_len

            if self.pam_position == '3_prime':
                # 3' PAM: 5'-Spacer(20nt)-PAM(3nt)-3'
                spacer_start = match_start - self.spacer_length
                if spacer_start < 0:
                    continue  # Not enough sequence upstream for a spacer

                # On minus strand, slice from the reversed sequence.
                # On plus strand, slice from the original sequence to preserve native U bases.
                if strand == '+':
                    spacer_seq = original_seq[spacer_start:match_start]
                    target_site = original_seq[spacer_start:match_end]
                    pam_match = original_seq[match_start:match_end]
                else:
                    spacer_seq = target_seq[spacer_start:match_start]
                    target_site = target_seq[spacer_start:match_end]
                    pam_match = pam_match_normalized
                
                strand_start = spacer_start
                strand_end = match_end

            else:
                # 5' PAM: 5'-PAM(4nt)-Spacer(23nt)-3'
                spacer_end = match_end + self.spacer_length
                if spacer_end > len(target_seq):
                    continue  # Not enough sequence downstream for a spacer

                if strand == '+':
                    spacer_seq = original_seq[match_end:spacer_end]
                    target_site = original_seq[match_start:spacer_end]
                    pam_match = original_seq[match_start:match_end]
                else:
                    spacer_seq = target_seq[match_end:spacer_end]
                    target_site = target_seq[match_start:spacer_end]
                    pam_match = pam_match_normalized

                strand_start = match_start
                strand_end = spacer_end

            # Map coordinates back to the original forward strand if on the minus strand
            if strand == '+':
                final_start = strand_start
                final_end = strand_end
            else:
                # On the reverse complement strand, index x corresponds to seq_len - 1 - x
                final_start = seq_len - strand_end
                final_end = seq_len - strand_start

            # Calculate 30nt context sequence (4nt upstream + spacer + PAM + 3nt downstream)
            context_seq_30 = None
            if self.pam_position == '3_prime' and self.spacer_length == 20 and len(pam_match_normalized) == 3:
                context_start = spacer_start - 4
                context_end = match_end + 3
                if context_start >= 0 and context_end <= len(target_seq):
                    # Slice from target_seq which represents the searched strand (already 5' to 3')
                    context_seq_30 = target_seq[context_start:context_end].upper()

            # Calculate score
            score = self.score_candidate(spacer_seq, target_site)

            candidate_info = {
                'spacer': spacer_seq.upper(),
                'pam': pam_match.upper(),
                'target_site': target_site.upper(),
                'start': final_start,
                'end': final_end,
                'strand': strand,
                'score': round(score, 3),
                'system_name': self.name,
                'context_seq_30': context_seq_30,
                'metadata': self.get_system_metadata()
            }
            strand_candidates.append(candidate_info)

        return strand_candidates

    def score_candidate(self, spacer: str, target_site: str) -> float:
        """
        Base scoring algorithm. Calculates GC content as a baseline structural score.
        Concrete classes can override this to implement advanced scoring (e.g. Azimuth).

        Args:
            spacer (str): Spacer sequence.
            target_site (str): Complete target site sequence.

        Returns:
            float: A score representing target efficiency (usually between 0.0 and 1.0 or 0 to 100).
        """
        if not spacer:
            return 0.0
        # GC content calculation: standard marker for secondary structure / stability
        gc_count = sum(1 for base in spacer.upper() if base in ('G', 'C'))
        gc_content = gc_count / len(spacer)
        
        # Penalyze extreme GC contents (ideal is 40% - 60%)
        # Returns a normalized score based on GC suitability
        deviation = abs(gc_content - 0.5)
        score = max(0.0, 100.0 - (deviation * 150.0))  # Scales ideal GC to 100, extreme GC to lower
        return score

    def get_system_metadata(self) -> Dict[str, Any]:
        """
        Returns a dictionary representing system-specific characteristics.
        """
        base_meta = {
            'pam': self.pam,
            'spacer_length': self.spacer_length,
            'pam_position': self.pam_position,
            'scoring_method': self.scoring_method
        }
        # Merge with any dynamic kwargs
        base_meta.update(self.extra_metadata)
        return base_meta
