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

## How to Run

```bash
# 1. Environment Setup
pip install -r requirements.txt

# 2. Execute Pipeline
python3 main.py
