#!/usr/bin/env python
"""
prepare_genome.py
~~~~~~~~~~~~~~~~~~
Downloads and unpacks a full local reference genome FASTA for one or more of
CRISPR-UNI's registered organisms, using NCBI's `datasets` command-line tool.

Why this exists: on-target guide scanning only ever needs a short user-supplied
sequence, so crispr_designer/organisms.py fetches small slices on demand via the
Ensembl REST API (or falls back to mock sequences). Off-target analysis is a
different problem - it needs the *entire* reference genome sitting locally so a
candidate guide can be searched against it. Whole-genome FASTA files are large
(tens of MB for E. coli up to several GB for human), so this is deliberately a
separate, explicit, one-time setup step rather than something that happens
silently as a side effect of a guide search.

Each organism already registered in crispr_designer/organisms.py has a curated
representative NCBI RefSeq assembly accession in REFERENCE_GENOME_ACCESSIONS
(e.g. Human -> GCF_000001405.40 / GRCh38.p14). This script resolves an organism
name to that accession, shells out to the NCBI `datasets` CLI to fetch the
genome data package (a zip file), and extracts just the genomic FASTA to
./genomes/<organism>/<organism>.fna - one subdirectory per organism (not a flat
./genomes/ folder), because off-target search tools like Cas-OFFinder scan every
FASTA file in a given directory; keeping organisms in separate folders is what
lets crispr_designer/offtarget.py point Cas-OFFinder at exactly one organism's
genome instead of accidentally searching everyone's genome at once.

Prerequisite: the NCBI `datasets` CLI must be installed and on PATH.
Install instructions: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/
(e.g. `conda install -c conda-forge ncbi-datasets-cli`)

Usage:
    python prepare_genome.py                # prepare every registered organism
    python prepare_genome.py Human Rice      # only these organisms (by common_name_en)
    python prepare_genome.py --list          # report status only, don't download
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile

from crispr_designer.organisms import (
    REFERENCE_GENOME_ACCESSIONS,
    get_local_genome_dir,
    get_local_genome_path,
)

DATASETS_INSTALL_URL = (
    "https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/"
)


def _check_datasets_cli() -> bool:
    return shutil.which("datasets") is not None


def _extract_genomic_fasta(zip_path: str, dest_path: str) -> bool:
    """
    Extracts the genomic FASTA file out of an NCBI Datasets genome data package zip and
    writes it to dest_path. Returns True on success.

    The package's internal layout (ncbi_dataset/data/<accession>/<accession>_<assembly>.fna)
    includes an assembly-name suffix we don't want to hardcode, so we just glob for the .fna
    file instead of assuming an exact filename.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(tmp_dir)

        matches = glob.glob(os.path.join(tmp_dir, "ncbi_dataset", "data", "*", "*.fna"))
        if not matches:
            return False

        os.makedirs(os.path.dirname(dest_path) or ".", exist_ok=True)
        shutil.move(matches[0], dest_path)
        return True


def prepare_organism(organism: str, accession: str, download: bool) -> str:
    """
    Returns a short status code: 'present', 'downloaded', 'missing', 'failed'.
    """
    dest_path = get_local_genome_path(organism)
    dest_dir = get_local_genome_dir(organism)

    if os.path.exists(dest_path):
        print(f"  [OK] {organism}: {dest_path} already present.")
        return "present"

    if not download:
        print(f"  [ ] {organism}: missing -> would download {accession} to {dest_path}")
        return "missing"

    if not _check_datasets_cli():
        print(f"  [!] {organism}: the NCBI `datasets` CLI is not installed or not on PATH.")
        print(f"      Install instructions: {DATASETS_INSTALL_URL}")
        return "failed"

    os.makedirs(dest_dir, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp_dir:
        zip_path = os.path.join(tmp_dir, f"{accession}.zip")
        print(f"  Downloading {organism} ({accession}) via `datasets` CLI ...")
        try:
            subprocess.run(
                [
                    "datasets", "download", "genome", "accession", accession,
                    "--include", "genome",
                    "--filename", zip_path,
                ],
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as e:
            print(f"  [!] {organism}: `datasets` download failed.")
            print(f"      {e.stderr.strip() if e.stderr else e}")
            return "failed"

        if not _extract_genomic_fasta(zip_path, dest_path):
            print(f"  [!] {organism}: downloaded package did not contain a .fna file as expected.")
            return "failed"

    print(f"  [+] {organism}: saved -> {dest_path}")
    return "downloaded"


def main():
    parser = argparse.ArgumentParser(
        description="Download full reference genomes for crispr_designer's registered organisms."
    )
    parser.add_argument(
        "organisms", nargs="*",
        help="Specific organism(s) to process, by common_name_en (default: all). "
             f"Available: {', '.join(REFERENCE_GENOME_ACCESSIONS.keys())}"
    )
    parser.add_argument("--list", action="store_true", help="Only report status; don't download anything.")
    args = parser.parse_args()

    wanted = REFERENCE_GENOME_ACCESSIONS
    if args.organisms:
        wanted = {name: acc for name, acc in REFERENCE_GENOME_ACCESSIONS.items() if name in args.organisms}
        unknown = set(args.organisms) - set(REFERENCE_GENOME_ACCESSIONS.keys())
        if unknown:
            print(f"Unknown organism(s): {', '.join(unknown)}")
            print(f"Available: {', '.join(REFERENCE_GENOME_ACCESSIONS.keys())}")
            sys.exit(1)

    print("Preparing crispr_designer reference genomes...\n")
    results = {}
    for organism, accession in wanted.items():
        results[organism] = prepare_organism(organism, accession, download=not args.list)

    print()
    needs_attention = [name for name, status in results.items() if status not in ("present", "downloaded")]
    if needs_attention:
        print(f"{len(needs_attention)} organism(s) still need attention: {', '.join(needs_attention)}")
        sys.exit(1)
    print("All requested reference genomes are present.")


if __name__ == "__main__":
    main()
