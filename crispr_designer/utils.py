"""
crispr_designer/utils.py
~~~~~~~~~~~~~~~~~~~~~~~~
Bioinformatic helper utilities for processing DNA/RNA sequences and PAM patterns.
Contains reverse complement calculations and IUPAC degenerate nucleotide pattern converters.
"""

import re

# IUPAC Degenerate Base mapping to standard regex groupings
IUPAC_MAP = {
    'A': 'A',
    'T': 'T',
    'C': 'C',
    'G': 'G',
    'U': 'T',         # Convert RNA U to DNA T for alignment/matching purposes
    'R': '[AG]',      # PuRine (A or G)
    'Y': '[CT]',      # PYrimidine (C or T)
    'S': '[GC]',      # Strong interaction (G or C)
    'W': '[AT]',      # Weak interaction (A or T)
    'K': '[GT]',      # Keto (G or T)
    'M': '[AC]',      # aMino (A or C)
    'B': '[CGT]',     # Not A (C, G, or T)
    'D': '[AGT]',     # Not C (A, G, or T)
    'H': '[ACT]',     # Not G (A, C, or T)
    'V': '[ACG]',     # Not T/U (A, C, or G)
    'N': '[ACGT]'     # aNy base (A, C, G, or T)
}

def reverse_complement(sequence: str) -> str:
    """
    Computes the reverse complement of a given DNA sequence.
    Preserves case of original nucleotides and handles standard degenerate bases.

    Args:
        sequence (str): The input DNA sequence (e.g., '5-ATCG-3').

    Returns:
        str: The reverse complement sequence in the same case.
    """
    complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
        'N': 'N', 'n': 'n', 'U': 'A', 'u': 'a',
        # Degenerates complement mapping
        'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
        'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
        'D': 'H', 'H': 'D', 'r': 'y', 'y': 'r',
        's': 's', 'w': 'w', 'k': 'm', 'm': 'k',
        'b': 'v', 'v': 'b', 'd': 'h', 'h': 'd'
    }
    return "".join(complement.get(base, base) for base in reversed(sequence))

def gc_content(sequence: str) -> float:
    """
    Calculates the GC ratio (0.0 - 1.0) of a nucleotide sequence.
    Shared helper so scoring functions across CRISPRSystem subclasses don't
    each reimplement the same GC-counting loop.

    Args:
        sequence (str): Nucleotide sequence (case-insensitive).

    Returns:
        float: GC ratio, or 0.0 for an empty sequence.
    """
    if not sequence:
        return 0.0
    seq = sequence.upper()
    gc_count = sum(1 for base in seq if base in ('G', 'C'))
    return gc_count / len(seq)

def iupac_to_regex(pattern: str) -> str:
    """
    Converts a degenerate IUPAC DNA sequence (like a PAM pattern) into a compiled regex string.
    Supports case-insensitive mapping.

    Example:
        'NGG' -> '[ACGT]GG'
        'NNGRRT' -> '[ACGT][ACGT][AG][AG][AG]T'

    Args:
        pattern (str): The IUPAC code string (e.g. 'NGG').

    Returns:
        str: A regular expression string representation.
    """
    regex_parts = []
    for char in pattern.upper():
        if char in IUPAC_MAP:
            regex_parts.append(IUPAC_MAP[char])
        else:
            # Fallback to literal if character is not standard IUPAC
            regex_parts.append(re.escape(char))
    return "".join(regex_parts)

def seq_to_onehot(seq: str) -> "np.ndarray":
    """
    Converts a DNA/RNA sequence into a multi-dimensional One-hot encoded Numpy array.
    Supports standard bases: A, C, G, T, U and degenerate base N.
    """
    import numpy as np
    mapping = {
        'A': [1.0, 0.0, 0.0, 0.0],
        'C': [0.0, 1.0, 0.0, 0.0],
        'G': [0.0, 0.0, 1.0, 0.0],
        'T': [0.0, 0.0, 0.0, 1.0],
        'U': [0.0, 0.0, 0.0, 1.0],
        'N': [0.25, 0.25, 0.25, 0.25]
    }
    encoded = [mapping.get(char.upper(), [0.0, 0.0, 0.0, 0.0]) for char in seq]
    return np.array(encoded, dtype=np.float32)
