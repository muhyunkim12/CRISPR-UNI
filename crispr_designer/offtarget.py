"""
crispr_designer/offtarget.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Off-target search for guide RNA candidates against a full local reference genome.

Unlike on-target scanning (base.py), which searches a short user-supplied sequence held
entirely in memory, off-target search needs a whole chromosome/genome-scale reference to
search against - see prepare_genome.py for how that reference gets downloaded locally
(one organism per subdirectory under ./genomes/).

Rather than reimplementing genome-scale approximate string matching (a suffix array / FM-index
style algorithm) from scratch, this wraps Cas-OFFinder (https://github.com/snugel/cas-offinder),
a widely used, actively maintained, OpenCL-accelerated open-source off-target search tool. It
must be installed separately and available on PATH as `cas-offinder` - see CasOffinderNotFoundError
below for install instructions.
"""

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

from .organisms import get_local_genome_dir, get_local_genome_path
from .cfd_scores import CFD_MISMATCH_SCORES, CFD_PAM_SCORES
from .sacas9_scores import SACAS9_MISMATCH_SCORES

CAS_OFFINDER_INSTALL_URL = "https://github.com/snugel/cas-offinder/releases"

# The CFD table only covers SpCas9's 20nt-spacer + 3' PAM geometry (see cfd_scores.py). Any
# system sharing that exact geometry (SpCas9 itself, and SpCas9-derived Prime Editor) can use
# the real published score; anything else falls back to the approximate heuristic below.
_CFD_APPLICABLE_SPACER_LENGTH = 20
_CFD_APPLICABLE_PAM_POSITION = '3_prime'

# The SaCas9 specificity model (see sacas9_scores.py) only covers SaCas9's fixed 21nt-spacer
# + 3' PAM geometry (crispr_designer.systems.SaCas9System).
_SACAS9_APPLICABLE_SPACER_LENGTH = 21
_SACAS9_APPLICABLE_PAM_POSITION = '3_prime'

_CFD_COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A'}


def _cfd_revcom(seq: str) -> str:
    return ''.join(_CFD_COMPLEMENT.get(base, base) for base in reversed(seq))


def calc_cfd_score(guide_spacer: str, hit_spacer: str, hit_pam: str) -> Optional[float]:
    """
    Computes the published CFD off-target score (Doench et al. 2016, table in cfd_scores.py)
    for one off-target hit: how likely SpCas9 is to still cut at a site that differs from the
    intended 20nt guide spacer. 1.0 means identical to a perfect on-target site; lower means
    less likely to be cut (safer).

    Args:
        guide_spacer: the intended (designed) 20nt spacer sequence.
        hit_spacer: the 20nt sequence actually found in the genome at this off-target site.
        hit_pam: the PAM found alongside hit_spacer (only the 2 PAM-proximal bases are used,
            matching the original CFD implementation).

    Returns:
        float score in [0.0, 1.0], or None if the inputs don't fit the table (wrong length,
        non-ACGU characters, or a PAM combination the table doesn't cover) - callers should
        treat None as "CFD isn't computable for this hit" and fall back to the heuristic.
    """
    guide = guide_spacer.upper().replace('T', 'U')
    hit = hit_spacer.upper().replace('T', 'U')
    if len(guide) != 20 or len(hit) != 20:
        return None

    # Note: unlike the spacer, the PAM key stays in DNA letters (T, not U) - CFD_PAM_SCORES
    # is keyed that way, matching the original reference implementation.
    pam_key = hit_pam.upper()[-2:]
    if pam_key not in CFD_PAM_SCORES:
        return None

    score = 1.0
    for i, (guide_base, hit_base) in enumerate(zip(guide, hit), start=1):
        if guide_base == hit_base:
            continue
        key = f'r{guide_base}:d{_cfd_revcom(hit_base)},{i}'
        if key not in CFD_MISMATCH_SCORES:
            return None
        score *= CFD_MISMATCH_SCORES[key]

    return score * CFD_PAM_SCORES[pam_key]


def calc_sacas9_score(guide_spacer: str, hit_spacer: str) -> Optional[float]:
    """
    Computes the published SaCas9 off-target specificity score (Tycko et al. 2018; table in
    sacas9_scores.py) for one off-target hit: predicted relative activity of SaCas9 at a site
    that differs from the intended 21nt guide spacer. 1.0 means identical to a perfect
    on-target site; lower means less likely to be cut (safer).

    Faithfully replicates the reference implementation (predict_activity_single.py from the
    paper's official repo): mismatches at spacer position 1 (5'-most, PAM-distal) are never
    scored (treated as if matched), and only single-nucleotide mismatches are supported (no
    bulges/gaps - this simply returns None if lengths don't match, same as calc_cfd_score).

    Args:
        guide_spacer: the intended (designed) 21nt SaCas9 spacer sequence.
        hit_spacer: the 21nt sequence actually found in the genome at this off-target site.

    Returns:
        float score in [0.0, 1.0], or None if the inputs don't fit the table (wrong length) -
        callers should treat None as "not computable" and fall back to the heuristic.
    """
    guide = guide_spacer.upper()
    hit = hit_spacer.upper()
    if len(guide) != _SACAS9_APPLICABLE_SPACER_LENGTH or len(hit) != _SACAS9_APPLICABLE_SPACER_LENGTH:
        return None

    score = 1.0
    for i, (hit_base, guide_base) in enumerate(zip(hit, guide), start=1):
        if i == 1:
            # Position 1 (5'-most, PAM-distal) mismatches aren't scored by this model - per
            # the original script: "Mismatches at the 5' position are not truly determined
            # by the model."
            continue
        if hit_base == guide_base:
            continue
        key = f'{hit_base}{guide_base},{i}'
        if key not in SACAS9_MISMATCH_SCORES:
            return None
        score *= SACAS9_MISMATCH_SCORES[key]

    return score


class CasOffinderNotFoundError(RuntimeError):
    """Raised when the `cas-offinder` binary isn't installed or isn't on PATH."""

    def __init__(self):
        super().__init__(
            "Cas-OFFinder is not installed or not on PATH. Off-target search relies on it "
            "instead of reimplementing genome-scale approximate string matching from scratch. "
            f"Install a prebuilt binary from: {CAS_OFFINDER_INSTALL_URL}"
        )


class ReferenceGenomeMissingError(FileNotFoundError):
    """Raised when the local reference genome for an organism hasn't been prepared yet."""

    def __init__(self, organism_name: str, genome_path: str):
        self.organism_name = organism_name
        self.genome_path = genome_path
        super().__init__(
            f"No local reference genome found for '{organism_name}' at {genome_path}. "
            f"Run `python prepare_genome.py {organism_name}` first."
        )


@dataclass
class OffTargetHit:
    """A single off-target match reported by Cas-OFFinder."""
    chromosome: str
    position: int
    strand: str
    mismatches: int
    matched_sequence: str


def is_cas_offinder_available() -> bool:
    return shutil.which("cas-offinder") is not None


def is_reference_genome_ready(organism_name: str) -> bool:
    return os.path.exists(get_local_genome_path(organism_name))


def _build_query_pattern(spacer_length: int, pam: str, pam_position: str) -> str:
    """Builds the Cas-OFFinder pattern line: spacer positions as 'N', PAM as literal IUPAC codes."""
    pam_pattern = pam.upper()
    if pam_position == '3_prime':
        return ('N' * spacer_length) + pam_pattern
    return pam_pattern + ('N' * spacer_length)


def search_offtargets(
    candidates: List[Dict[str, Any]],
    organism_name: str,
    spacer_length: int,
    pam: str,
    pam_position: str = '3_prime',
    max_mismatches: int = 3,
    device: str = 'C',
    timeout_sec: int = 1800,
) -> Dict[str, List[OffTargetHit]]:
    """
    Runs Cas-OFFinder against the organism's local reference genome for each candidate's
    full target_site (spacer + PAM), returning a dict keyed by target_site -> list of
    OffTargetHit found anywhere in the genome. This includes the guide's own on-target site
    (a 0-mismatch hit) - callers that want off-target hits only should filter that one out
    themselves (e.g. by comparing chromosome/position against where the guide was found).

    Args:
        candidates: guide candidate dicts from CRISPRSystem.find_candidates() (each must have
            a 'target_site' key - spacer + PAM concatenated, already in genomic 5'->3' order).
        organism_name: any registered organism alias/common_name_en (see ORGANISM_REGISTRY).
        spacer_length: spacer length for the active CRISPR system.
        pam: IUPAC PAM pattern for the active CRISPR system.
        pam_position: '3_prime' or '5_prime'.
        max_mismatches: maximum mismatches Cas-OFFinder should tolerate per site.
        device: 'C' for CPU, 'G'/'G0'/etc. for a GPU device if Cas-OFFinder was built with
            OpenCL GPU support. CPU is the safe default since it needs no extra setup.
        timeout_sec: safety timeout for the Cas-OFFinder subprocess.

    Returns:
        Dict[str, List[OffTargetHit]], keyed by the (uppercased) target_site sequence.

    Raises:
        CasOffinderNotFoundError: if the `cas-offinder` binary isn't on PATH.
        ReferenceGenomeMissingError: if prepare_genome.py hasn't been run for this organism yet.
    """
    if not is_cas_offinder_available():
        raise CasOffinderNotFoundError()

    genome_dir = get_local_genome_dir(organism_name)
    if not is_reference_genome_ready(organism_name):
        raise ReferenceGenomeMissingError(organism_name, get_local_genome_path(organism_name))

    pattern = _build_query_pattern(spacer_length, pam, pam_position)
    expected_len = len(pattern)

    # Cas-OFFinder requires every query line to match the pattern's length exactly, and de-dupes
    # naturally if we only submit each distinct target_site once.
    query_sequences: List[str] = []
    seen = set()
    for cand in candidates:
        target_site = cand['target_site'].upper()
        if len(target_site) != expected_len or target_site in seen:
            continue
        seen.add(target_site)
        query_sequences.append(target_site)

    results: Dict[str, List[OffTargetHit]] = {seq: [] for seq in query_sequences}
    if not query_sequences:
        return results

    with tempfile.TemporaryDirectory() as tmp_dir:
        input_path = os.path.join(tmp_dir, "offtarget_input.txt")
        output_path = os.path.join(tmp_dir, "offtarget_output.txt")

        with open(input_path, "w") as f:
            f.write(os.path.abspath(genome_dir) + "\n")
            f.write(pattern + "\n")
            for seq in query_sequences:
                f.write(f"{seq} {max_mismatches}\n")

        try:
            subprocess.run(
                ["cas-offinder", input_path, device, output_path],
                check=True,
                capture_output=True,
                text=True,
                timeout=timeout_sec,
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"Cas-OFFinder failed: {e.stderr.strip() if e.stderr else e}"
            ) from e
        except subprocess.TimeoutExpired as e:
            raise RuntimeError(
                f"Cas-OFFinder timed out after {timeout_sec}s (genome too large, or "
                f"max_mismatches={max_mismatches} too permissive for this many candidates)."
            ) from e

        if os.path.exists(output_path):
            with open(output_path) as f:
                for line in f:
                    parts = line.rstrip("\n").split("\t")
                    # Cas-OFFinder output columns: query_seq, chrom, 0-based pos, matched_seq,
                    # strand, mismatch_count.
                    if len(parts) < 6:
                        continue
                    query_seq, chrom, pos, matched_seq, strand, mismatches = parts[:6]
                    query_seq = query_seq.upper()
                    if query_seq not in results:
                        continue
                    results[query_seq].append(
                        OffTargetHit(
                            chromosome=chrom,
                            position=int(pos),
                            strand=strand,
                            mismatches=int(mismatches),
                            matched_sequence=matched_seq.upper(),
                        )
                    )

    return results


def _heuristic_risk_score(hits: List[OffTargetHit]) -> float:
    """
    Fallback used only when the real CFD table doesn't apply (see offtarget_risk_score).
    Simplified, clearly-approximate off-target risk score from 0 (many close off-target
    matches - risky) to 100 (no close off-target matches found - clean). This is a rough
    mismatch-count-weighted heuristic - it does not account for *where* in the spacer a
    mismatch falls (real off-target risk depends heavily on whether mismatches land in the
    PAM-proximal "seed" region vs. the PAM-distal end), only on how many mismatches there are.
    """
    if not hits:
        return 100.0

    penalty = 0.0
    for hit in hits:
        if hit.mismatches == 0:
            penalty += 40.0
        else:
            penalty += max(0.0, 20.0 - (hit.mismatches * 5.0))

    return max(0.0, 100.0 - penalty)


def offtarget_risk_score(
    hits: List[OffTargetHit],
    guide_spacer: Optional[str] = None,
    spacer_length: Optional[int] = None,
    pam_position: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Aggregates a guide's off-target hits into a single specificity score (0-100, higher is
    safer/cleaner).

    If the CRISPR system's geometry matches what the published CFD (Cutting Frequency
    Determination) score was empirically fit to - a 20nt spacer with a 3' PAM, i.e. SpCas9
    and SpCas9-derived systems like Prime Editor - this computes the real per-hit CFD score
    (Doench et al. 2016; see cfd_scores.py) and aggregates them with the same formula used by
    established tools like FlashFry:

        specificity = 100 / (1 + sum(CFD score of every off-target hit))

    If instead the geometry matches SaCas9's fixed 21nt spacer + 3' PAM, this uses the real
    published SaCas9 specificity model (Tycko et al. 2018; see sacas9_scores.py) with the same
    aggregation formula.

    For any other system (Cas12a, Cas14a, or a custom-registered one), no published,
    portable scoring table is known to be available, so this falls back to the approximate
    _heuristic_risk_score instead of misapplying a table fit to a different system's geometry.

    Args:
        hits: off-target hits for one guide candidate, from search_offtargets().
        guide_spacer: the candidate's own spacer sequence (cand['spacer']). Required for CFD
            and for the SaCas9 model.
        spacer_length: the active system's spacer_length. Required for CFD and SaCas9.
        pam_position: the active system's pam_position. Required for CFD and SaCas9.

    Returns:
        {'score': float, 'method': 'CFD' | 'SaCas9-specificity' | 'heuristic'}

    Caveat: if the guide's own on-target sequence also exists verbatim elsewhere in the
    reference genome (or, as is normal, at its own true genomic location), that shows up as a
    0-mismatch hit indistinguishable from a real duplicate site - this module has no way to
    map the guide back to "where it was actually meant to cut" in genome coordinates, so the
    resulting score is a conservative (slightly pessimistic) estimate in that case.
    """
    cfd_applicable = (
        guide_spacer is not None
        and spacer_length == _CFD_APPLICABLE_SPACER_LENGTH
        and pam_position == _CFD_APPLICABLE_PAM_POSITION
    )

    if cfd_applicable:
        cfd_scores = []
        for hit in hits:
            matched = hit.matched_sequence.upper()
            if len(matched) < spacer_length + 2:
                cfd_applicable = False
                break
            cfd = calc_cfd_score(guide_spacer, matched[:spacer_length], matched[spacer_length:])
            if cfd is None:
                cfd_applicable = False
                break
            cfd_scores.append(cfd)

        if cfd_applicable:
            specificity = 100.0 / (1.0 + sum(cfd_scores))
            return {'score': specificity, 'method': 'CFD'}

    sacas9_applicable = (
        guide_spacer is not None
        and spacer_length == _SACAS9_APPLICABLE_SPACER_LENGTH
        and pam_position == _SACAS9_APPLICABLE_PAM_POSITION
    )

    if sacas9_applicable:
        sacas9_scores = []
        for hit in hits:
            matched = hit.matched_sequence.upper()
            if len(matched) < spacer_length:
                sacas9_applicable = False
                break
            sa_score = calc_sacas9_score(guide_spacer, matched[:spacer_length])
            if sa_score is None:
                sacas9_applicable = False
                break
            sacas9_scores.append(sa_score)

        if sacas9_applicable:
            specificity = 100.0 / (1.0 + sum(sacas9_scores))
            return {'score': specificity, 'method': 'SaCas9-specificity'}

    return {'score': _heuristic_risk_score(hits), 'method': 'heuristic'}
