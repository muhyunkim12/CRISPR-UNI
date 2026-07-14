<div align="right">
  <a href="README.md">🇰🇷 한국어</a> | 
  <a href="README.en.md">🇺🇸 English</a> | 
  <a href="README.jp.md">🇯🇵 日本語</a>
</div>

# CRISPR-UNI 
🧬 Next-Gen CRISPR & Prime Editing Multi-Model Design Pipeline (CLI)

This project is a high-performance CLI pipeline designed for laboratory workstations. It is engineered to comprehend the distinct molecular biology mechanisms of various CRISPR systems (e.g., SpCas9, Cas12a, Prime Editor) and automatically route them to the most optimized State-of-the-Art (SOTA) deep learning prediction models (e.g., DeepSpCas9, PRIDICT) to maximize genome editing efficiency.

## Core Architectural Features

- **AI Modeling (Factory Pattern):** Instead of relying on hardcoded conditional statements (if/else), an Object-Oriented Programming (OOP) design pattern is implemented. The `PredictorFactory` dynamically allocates dedicated deep learning inference objects tailored to the specific traits of the selected nuclease (DNA double-strand break, RNA targeting, reverse transcription, etc.).

- **Multi-Dimensional Optimization for Prime Editing:** To overcome the inherently low initial efficiency of Prime Editing (PE), the pipeline automatically generates a 3x3 combinatorial matrix of pegRNAs based on the nicking site—varying PBS lengths (10, 13, 15nt) and RT template lengths (12, 15, 20nt). Final scoring and optimization are executed via the Transformer-based PRIDICT architecture.

- **Defensive Programming:** Pipeline stability is fortified through rigorous exception handling. This includes strict input sequence validation (filtering non-ATGC characters) and verification logic to check the presence of local model weight files within the server infrastructure.

- **Off-target Analysis (optional):** Integrates [Cas-OFFinder](https://github.com/snugel/cas-offinder) against a locally prepared reference genome to check how closely each guide candidate matches other locations in the genome. See the "Off-target Analysis" section below for setup.

### PREREQUISITE(Important)

For the AI predictions of this pipeline to function correctly, the pre-trained weights for each model must be located in the `weights/` directory. 

The easiest way is to run the helper script below. It auto-downloads any model with a confirmed source URL and saves it under the exact filename the code expects.

```bash
python download_weights.py
```

To fetch them manually instead, use the commands below. Note: the filename after `-O` matters - it must match what the code expects, which differs from the upstream source's own filename.

```bash
# Create the directory for weight files and navigate into it.
mkdir -p weights && cd weights

# 1. DeepSpCas9 Weights (for SpCas9)
# Source: Yonsei Univ. Kim Lab Official GitHub
wget -O DeepSpCas9_weights.h5 "https://raw.githubusercontent.com/myungjinkim/DeepSpCas9/master/DeepSpCas9_model.h5"

# 2. DeepCpf1 Weights (for Cas12a)
# Source: Yonsei Univ. Kim Lab Official GitHub
wget -O DeepCpf1_weights.h5 "https://raw.githubusercontent.com/myungjinkim/DeepCpf1/master/DeepCpf1_model.h5"

# 3. PRIDICT Weights (for Prime Editor)
# Source: PRIDICT Official Zenodo Archive (Author's distribution link)
# Note: The PRIDICT model is large, so the download may take some time.
wget -O PRIDICT_model.pt "https://zenodo.org/record/8208465/files/pridict_v2_model.pt"

# Return to the parent directory once the downloads are complete.
cd ..
```

## Off-target Analysis (optional)

By default, the tool only performs an on-target scan within the sequence you provide. To also check genome-wide off-target risk for real experimental design, two extra setup steps are needed.

**Step 1: Prepare the reference genome (once per organism)**

Off-target search needs the organism's whole-genome FASTA available locally. `prepare_genome.py` fetches it automatically via the [NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/) (installable via `conda install -c conda-forge ncbi-datasets-cli`).

```bash
# Prepare genomes for every registered organism (Human, Mouse, Rice, Arabidopsis, Yeast, Ecoli)
python prepare_genome.py

# Or just one organism
python prepare_genome.py Human

# Check status only, without downloading
python prepare_genome.py --list
```

Downloaded genomes are saved under `genomes/<organism>/` and are already excluded via `.gitignore`.

**Step 2: Install Cas-OFFinder**

Rather than reimplementing genome-scale off-target matching, this relies on the established open-source tool [Cas-OFFinder](https://github.com/snugel/cas-offinder/releases). Download the binary and make sure `cas-offinder` is on your PATH.

**Step 3: Run with `--check-offtargets`**

```bash
python3 main.py --check-offtargets

# Adjust the allowed mismatch count (default: 3)
python3 main.py --check-offtargets --max-mismatches 2
```

If either prerequisite (reference genome or Cas-OFFinder) is missing, the tool won't crash - it prints a clear message and still shows on-target results normally.

> **Note:** For SpCas9 and Prime_Editor (both a 20nt spacer + 3' PAM), the results table shows the real published CFD (Cutting Frequency Determination, Doench et al. 2016) score - labeled "CFD". For SaCas9 (21nt spacer + 3' PAM), it shows the real published SaCas9 specificity model (Tycko et al. 2018, Nat Commun) - labeled "SaCas9-specificity". For every other system (Cas12a, Cas14a), no verified, portable scoring table is known to be available, so it automatically falls back to a simplified mismatch-count-only approximation labeled "heuristic" - which does not weight mismatches by their position relative to the PAM (seed region).

## How to Run

```bash
# 1. Environment Setup
pip install -r requirements.txt

# 2. Execute Pipeline
python3 main.py
```

```Result
==================================================
  CRISPR & PE Multi-Model Design Pipeline v1.0
==================================================
[1] SpCas9       (PAM: NGG    | Target: DNA)
[2] Cas12a       (PAM: TTTN   | Target: DNA)
[3] Prime_Editor (PAM: NGG    | Target: DNA)
...
Select CRISPR type (1-5): 3
Insert Target DNA Sequence: ATCG...
>> Starting [PRIDICT Transformer Model] based pegRNA Multi-Dimensional Optimization...
```

###  Manual Package Installation

`requirements.txt` only installs the one hard dependency (numpy). The deep learning frameworks (TensorFlow, PyTorch) are deliberately left out of the default install since they're heavy and you only need whichever one matches the CRISPR system(s) you actually plan to use. Install just what you need below. (Note: the tool still runs without any framework installed - on-target scanning and the GC-based mock scores work regardless; a framework is only required when you actually invoke a real deep learning predictor.)

```bash
# 1. Required dependency
pip install numpy

# 2. For SpCas9 and Cas12a predictors (DeepSpCas9, DeepCpf1)
# Requires TensorFlow.
pip install tensorflow

# 3. For the Prime Editor predictor (PRIDICT)
# Requires PyTorch.
pip install torch
```

> **Note:** Cas13d (RNA-targeting) is not supported by this tool. The real published Cas13d efficiency predictor (Wessels et al. 2020, Sanjana lab, [gitlab.com/sanjanalab/cas13](https://gitlab.com/sanjanalab/cas13)) is not a downloadable PyTorch weight file - it's a separate R + RNAhybrid + ViennaRNA pipeline, which doesn't fit this project's "download weights, then score instantly" architecture, so it was left out.
>
> **Note:** SaCas9 and Cas14a(Cas12f) support guide-finding, but have no deep learning efficiency predictor. A search turned up no real, publicly downloadable pretrained model under either name (the download URLs previously in this codebase pointed to repositories that don't exist). Their results table instead shows a GC-content based heuristic score, labeled "Heuristic (GC-based) Score".

### Feedback & Support

If you encounter any issues, have feature requests, or would like to provide feedback regarding this tool, please feel free to reach out to the contact below. 

* **Email:** kmh40283656@gmail.com
