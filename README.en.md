<div align="right">
  <a href="README.md">🇰🇷 한국어</a> | 
  <a href="README.en.md">🇺🇸 English</a> | 
  <a href="README.ja.md">🇯🇵 日本語</a>
</div>

# CRISPR-UNI 
🧬 Next-Gen CRISPR & Prime Editing Multi-Model Design Pipeline (CLI)

This project is a high-performance CLI pipeline designed for laboratory workstations. It is engineered to comprehend the distinct molecular biology mechanisms of various CRISPR systems (e.g., SpCas9, Cas12a, Prime Editor) and automatically route them to the most optimized State-of-the-Art (SOTA) deep learning prediction models (e.g., DeepSpCas9, PRIDICT) to maximize genome editing efficiency.

## Core Architectural Features

- **AI Modeling (Factory Pattern):** Instead of relying on hardcoded conditional statements (if/else), an Object-Oriented Programming (OOP) design pattern is implemented. The `PredictorFactory` dynamically allocates dedicated deep learning inference objects tailored to the specific traits of the selected nuclease (DNA double-strand break, RNA targeting, reverse transcription, etc.).

- **Multi-Dimensional Optimization for Prime Editing:** To overcome the inherently low initial efficiency of Prime Editing (PE), the pipeline automatically generates a 3x3 combinatorial matrix of pegRNAs based on the nicking site—varying PBS lengths (10, 13, 15nt) and RT template lengths (12, 15, 20nt). Final scoring and optimization are executed via the Transformer-based PRIDICT architecture.

- **Defensive Programming:** Pipeline stability is fortified through rigorous exception handling. This includes strict input sequence validation (filtering non-ATGC characters) and verification logic to check the presence of local model weight files within the server infrastructure.

### PREREQUISITE(Important)

For the AI predictions of this pipeline to function correctly, the pre-trained weights for each model must be located in the `weights/` directory. 

```bash
# Create the directory for weight files and navigate into it.
mkdir -p weights && cd weights

# 1. DeepSpCas9 Weights (for SpCas9)
# Source: Yonsei Univ. Kim Lab Official GitHub
wget -O DeepSpCas9_model.h5 "[https://raw.githubusercontent.com/myungjinkim/DeepSpCas9/master/DeepSpCas9_model.h5](https://raw.githubusercontent.com/myungjinkim/DeepSpCas9/master/DeepSpCas9_model.h5)"

# 2. DeepCpf1 Weights (for Cas12a)
# Source: Yonsei Univ. Kim Lab Official GitHub
wget -O DeepCpf1_model.h5 "[https://raw.githubusercontent.com/myungjinkim/DeepCpf1/master/DeepCpf1_model.h5](https://raw.githubusercontent.com/myungjinkim/DeepCpf1/master/DeepCpf1_model.h5)"

# 3. PRIDICT Weights (for Prime Editor)
# Source: PRIDICT Official Zenodo Archive (Author's distribution link)
# Note: The PRIDICT model is large, so the download may take some time.
wget -O PRIDICT_model.pt "[https://zenodo.org/record/8208465/files/pridict_v2_model.pt](https://zenodo.org/record/8208465/files/pridict_v2_model.pt)"

# Return to the parent directory once the downloads are complete.
cd ..
```

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
Select CRISPR type (1-6): 3
Insert Target DNA Sequence: ATCG...
>> Starting [PRIDICT Transformer Model] based pegRNA Multi-Dimensional Optimization...
```

###  Manual Package Installation

Depending on your local workstation environment, `pip install -r requirements.txt` might fail due to dependency conflicts. If this occurs, install the core dependencies first, and then manually install the specific AI frameworks required for your target CRISPR system.

```bash
# 1. Core Dependencies (Sequence processing & Data manipulation)
pip install biopython pandas numpy scikit-learn

# 2. For SpCas9, SaCas9, and Cas12a Models (DeepSpCas9, DeepCpf1)
# Requires TensorFlow for deep learning inference.
pip install tensorflow keras

# 3. For Prime Editor (PE) Models (PRIDICT)
# Requires PyTorch and Hugging Face Transformers.
pip install torch torchvision torchaudio transformers

# 4. Utility for downloading pre-trained weights
pip install gdown
```
