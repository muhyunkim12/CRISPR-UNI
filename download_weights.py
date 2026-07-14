#!/usr/bin/env python
"""
download_weights.py
~~~~~~~~~~~~~~~~~~~~
Fetches the deep learning model weight files used by crispr_designer.predictors.

Model weights are intentionally excluded from git (see .gitignore) - they're too large
to commit safely. Rather than making every user read the README and manually place files
under ./weights/, this script automates that:

  1. Reads (model_name, weight_file) straight from the predictor classes in
     crispr_designer/predictors.py, so the expected filename/path always matches exactly
     what the running code looks for - one source of truth, no second hand-maintained list.
  2. Skips any weight file that already exists locally.
  3. For models with a known, real, direct-download URL (see KNOWN_DIRECT_URLS below),
     downloads the file automatically and saves it under the exact filename the code
     expects, regardless of what the source names it.
  4. For everything else, the predictor's own `download_url` (the same one shown in the
     app's ModelWeightsMissingError messages) is a project homepage / "more info" link, not
     a direct file - so those are reported with manual instructions instead of being
     silently mishandled as if they were downloadable.
  5. Optionally verifies a sha256 checksum, if one is configured in KNOWN_CHECKSUMS.

Usage:
    python download_weights.py                    # check/download every known model
    python download_weights.py DeepSpCas9 PRIDICT  # only these models (matched by model_name)
    python download_weights.py --list              # report status only, don't download
"""

import argparse
import hashlib
import os
import sys
import urllib.request

from crispr_designer.predictors import (
    DeepSpCas9Predictor,
    DeepCpf1Predictor,
    PRIDICTPredictor,
)

# Real, verified direct-download URLs, keyed by the exact filename crispr_designer/predictors.py
# expects under weights/. An entry here means "safe to auto-download". Anything missing here
# falls back to manual instructions using the predictor's own download_url (a project homepage,
# not a direct file link) - we don't guess at URLs we haven't confirmed actually serve the file.
KNOWN_DIRECT_URLS = {
    "DeepSpCas9_weights.h5": "https://raw.githubusercontent.com/myungjinkim/DeepSpCas9/master/DeepSpCas9_model.h5",
    "DeepCpf1_weights.h5": "https://raw.githubusercontent.com/myungjinkim/DeepCpf1/master/DeepCpf1_model.h5",
    "PRIDICT_model.pt": "https://zenodo.org/record/8208465/files/pridict_v2_model.pt",
}

# Optional sha256 checksums for integrity verification, keyed the same way as KNOWN_DIRECT_URLS.
# Left empty until real hashes are pinned down - add entries like "DeepSpCas9_weights.h5": "<hash>".
KNOWN_CHECKSUMS = {}

ALL_PREDICTOR_CLASSES = [
    DeepSpCas9Predictor,
    DeepCpf1Predictor,
    PRIDICTPredictor,
]


def _sha256_of(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _verify_checksum(path: str, filename: str) -> bool:
    expected = KNOWN_CHECKSUMS.get(filename)
    if not expected:
        return True  # nothing configured to check against - treat as fine
    actual = _sha256_of(path)
    if actual != expected:
        print(f"      [!] checksum mismatch: expected {expected[:12]}..., got {actual[:12]}...")
        return False
    return True


def check_and_download(predictor, download: bool) -> str:
    """
    Returns a short status code: 'present', 'downloaded', 'manual_required',
    'missing', 'failed', or 'checksum_mismatch'.
    """
    path = predictor.weight_file
    filename = os.path.basename(path)
    name = predictor.model_name

    if os.path.exists(path):
        if _verify_checksum(path, filename):
            print(f"  [OK] {name}: {path} already present.")
            return "present"
        print(f"  [!] {name}: {path} present but failed checksum verification.")
        return "checksum_mismatch"

    direct_url = KNOWN_DIRECT_URLS.get(filename)

    if not direct_url:
        print(f"  [ ] {name}: not found locally, and no known direct-download URL is configured.")
        print(f"      More info / source: {predictor.download_url}")
        print(f"      Please download the weights manually and save them as: {path}")
        return "manual_required"

    if not download:
        print(f"  [ ] {name}: missing -> would download from {direct_url}")
        return "missing"

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    print(f"  Downloading {name} from {direct_url} ...")
    try:
        req = urllib.request.Request(direct_url, headers={"User-Agent": "crispr-uni-weight-downloader"})
        with urllib.request.urlopen(req, timeout=120) as response, open(path, "wb") as out_file:
            out_file.write(response.read())
    except Exception as e:
        print(f"  [!] {name}: download failed ({e}).")
        print(f"      Please download manually from: {direct_url}")
        print(f"      and save it as: {path}")
        return "failed"

    if not _verify_checksum(path, filename):
        os.remove(path)
        print(f"  [!] {name}: removed downloaded file due to checksum mismatch.")
        return "failed"

    print(f"  [+] {name}: downloaded -> {path}")
    return "downloaded"


def main():
    parser = argparse.ArgumentParser(description="Check/download crispr_designer model weight files.")
    parser.add_argument("models", nargs="*", help="Specific model_name(s) to process (default: all).")
    parser.add_argument("--list", action="store_true", help="Only report status; don't download anything.")
    args = parser.parse_args()

    predictors = [cls() for cls in ALL_PREDICTOR_CLASSES]

    if args.models:
        wanted = {m.lower() for m in args.models}
        predictors = [p for p in predictors if p.model_name.lower() in wanted]
        if not predictors:
            print(f"No matching model(s) for: {', '.join(args.models)}")
            sys.exit(1)

    print("Checking crispr_designer model weights...\n")
    results = {}
    for predictor in predictors:
        results[predictor.model_name] = check_and_download(predictor, download=not args.list)

    print()
    needs_attention = [name for name, status in results.items() if status not in ("present", "downloaded")]
    if needs_attention:
        print(f"{len(needs_attention)} model(s) still need attention: {', '.join(needs_attention)}")
        sys.exit(1)
    print("All requested model weights are present.")


if __name__ == "__main__":
    main()
