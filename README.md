<div align="right">
  <a href="README.md">🇰🇷 한국어</a> | 
  <a href="README.en.md">🇺🇸 English</a> | 
  <a href="README.jp.md">🇯🇵 日本語</a>
</div>

# CRISPR-UNI 
🧬 Next-Gen CRISPR & Prime Editing Multi-Model Design Pipeline (CLI)
본 프로젝트는 다양한 유전자 가위(SpCas9, Cas12a, Prime Editor 등)의 분자생물학적 메커니즘 차이를 이해하고, 각 시스템에 최적화된 최신 SOTA 딥러닝 예측 모델(DeepSpCas9, PRIDICT 등)을 자동으로 매칭하여 교정 효율을 극대화하는 실험실 워크스테이션용 고성능 CLI 파이프라인입니다.

# 핵심 아키텍처 특징

- AI Modeling: 사용자가 유전자 가위를 선택하면 조건문으로 코드를 기우는 대신, Predictor가 가위의 특성(DNA 절단, RNA 타겟팅, 역전사 등)에 맞는 전용 딥러닝 추론 객체를 동적으로 할당하는 객체 지향 디자인 패턴을 적용했습니다.

- Multi-Dimensional Optimization for Prime Editing: 초기 효율이 낮은 프라임 에디팅(PE)의 한계를 극복하기 위해, 절단 위치 기준 PBS(10, 13, 15nt)와 RT 템플릿(12, 15, 20nt)의 3x3 가이드 RNA 조합을 자동 생성하고, 트랜스포머(Transformer) 기반 PRIDICT 아키텍처를 통해 최종 스코어링을 수행합니다.

- Defensive Programming: 입력 서열의 유효성 검사(A, T, G, C 외 예외 처리) 및 서버 인프라 내 가중치 파일 유무를 체크하는 예외 처리 로직을 촘꼼히 반영하여 파이프라인의 안정성을 높였습니다.

- Off-target 분석 (선택 사항): 로컬에 준비된 참조 게놈을 대상으로 [Cas-OFFinder](https://github.com/snugel/cas-offinder)를 연동하여, 각 가이드 후보가 게놈 내 다른 위치와 얼마나 유사하게 매칭되는지 검사합니다. 자세한 사용법은 아래 "오프타깃 분석" 섹션을 참고하세요.

# PREREQUISITE (Important)

본 파이프라인의 AI 예측을 정상적으로 구동하기 위해서는 각 모델의 사전 학습된 가중치(Weights) 파일이 `weights/` 디렉토리 내에 존재해야 합니다.

**가중치 파일 외에 해당하는 딥러닝 프레임워크(TensorFlow 및/또는 PyTorch)도 별도로 설치해야 합니다.** `requirements.txt`는 용량 문제로 이 프레임워크들을 기본 설치에 포함하지 않으며, 사용하려는 시스템에 따라 둘 중 하나만 있으면 됩니다. 시스템별로 정확히 어떤 프레임워크가 필요한지는 아래 "패키지 수동 설치 가이드" 섹션을 참고하세요.

가중치 파일 자체를 받는 가장 쉬운 방법은 아래 스크립트를 실행하는 것입니다. 이미 확인된 다운로드 링크가 있는 모델은 자동으로 받아서 코드가 기대하는 정확한 파일명으로 저장해 줍니다.

```bash
python download_weights.py
```

수동으로 받고 싶다면 아래 명령어를 사용하세요. (주의: `-O` 뒤 파일명은 반드시 아래와 같이 지정해야 합니다. 소스 파일명과 크립스퍼 디자이너 코드가 기대하는 파일명이 다릅니다.)

```bash
# 가중치 파일이 저장될 폴더를 만들고 이동합니다.
mkdir -p weights && cd weights

# 1. DeepSpCas9 가중치 (SpCas9 용)
# 출처: Yonsei Univ. Kim Lab Official GitHub
wget -O DeepSpCas9_weights.h5 "https://raw.githubusercontent.com/myungjinkim/DeepSpCas9/master/DeepSpCas9_model.h5"

# 2. DeepCpf1 가중치 (Cas12a 용)
# 출처: Yonsei Univ. Kim Lab Official GitHub
wget -O DeepCpf1_weights.h5 "https://raw.githubusercontent.com/myungjinkim/DeepCpf1/master/DeepCpf1_model.h5"

# 3. PRIDICT 가중치 (Prime Editor 용)
# 출처: PRIDICT 공식 Zenodo 아카이브 (원작자 배포 링크)
# 주의: PRIDICT 모델은 용량이 크므로 다운로드에 시간이 걸릴 수 있습니다.
wget -O PRIDICT_model.pt "https://zenodo.org/record/8208465/files/pridict_v2_model.pt"

# 다운로드가 완료되면 다시 상위 폴더로 빠져나옵니다.
cd ..
```

# 오프타깃 분석 (선택 사항)

기본 실행에서는 입력한 서열 안에서 후보를 찾는 온타깃 스캔만 수행합니다. 실제 실험 설계를 위해 게놈 전체를 대상으로 오프타깃(off-target) 위험도까지 확인하려면 아래 두 단계가 추가로 필요합니다.

**1단계: 참조 게놈 준비 (조직당 1회만)**

오프타깃 검사는 조직의 전장 게놈 파일이 로컬에 있어야 합니다. `prepare_genome.py`가 [NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/)를 통해 이를 자동으로 받아줍니다 (`conda install -c conda-forge ncbi-datasets-cli`로 설치 가능).

```bash
# 등록된 모든 유기체(Human, Mouse, Rice, Arabidopsis, Yeast, Ecoli)의 게놈을 준비
python prepare_genome.py

# 특정 유기체만 준비하고 싶다면
python prepare_genome.py Human

# 다운로드 없이 상태만 확인
python prepare_genome.py --list
```

받은 게놈은 `genomes/<유기체>/` 폴더에 저장되며, `.gitignore`에 이미 포함되어 있어 커밋되지 않습니다.

**2단계: Cas-OFFinder 설치**

실제 오프타깃 매칭 검색은 직접 구현하는 대신 검증된 오픈소스 툴인 [Cas-OFFinder](https://github.com/snugel/cas-offinder/releases)를 사용합니다. 바이너리를 받아 `cas-offinder` 명령으로 실행 가능하도록 PATH에 등록해주세요.

**3단계: `--check-offtargets` 옵션으로 실행**

```bash
python3 main.py --check-offtargets

# 허용 미스매치 개수도 조절 가능 (기본값 3)
python3 main.py --check-offtargets --max-mismatches 2
```

두 준비물(참조 게놈, Cas-OFFinder) 중 하나라도 빠져 있으면 프로그램이 멈추는 대신 안내 메시지를 띄우고 온타깃 결과만 정상적으로 보여줍니다.

> **참고**: SpCas9와 Prime_Editor(둘 다 20nt 스페이서 + 3' PAM 구조)는 실제 논문에 발표된 CFD(Cutting Frequency Determination, Doench et al. 2016) 점수를, SaCas9(21nt 스페이서 + 3' PAM)는 실제 논문에 발표된 SaCas9 특이도 모델(Tycko et al. 2018, Nat Commun)을 그대로 계산해서 사용합니다 — 결과 테이블에 각각 "CFD", "SaCas9-specificity"로 표시됩니다. 그 외 시스템(Cas12a, Cas14a)은 접근 가능한 검증된 매트릭스가 없어 미스매치 개수만 반영하는 단순 근사 점수로 자동 대체되며 "heuristic"으로 표시됩니다. 이 경우엔 PAM에 가까운 미스매치가 더 치명적이라는 seed region 가중치가 반영되어 있지 않습니다.

# How to Run

```bash
# 1. 환경 세팅
pip install -r requirements.txt

# 2. 파이프라인 실행
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
가위 종류를 선택하세요 (1-5): 3
타겟 DNA 서열을 입력하세요: ATCG...
>> [PRIDICT Transformer Model] 기반 pegRNA 최적화 조합 연산을 시작합니다...
```

### 패키지 수동 설치 가이드 

`requirements.txt`는 필수 의존성(numpy) 하나만 설치합니다. 딥러닝 프레임워크(TensorFlow, PyTorch)는 무거운 데다 실제로 사용할 CRISPR 시스템에 따라 둘 중 하나만 있으면 되기 때문에 기본 설치에서는 빠져 있습니다. 사용하려는 시스템에 맞춰 아래에서 필요한 것만 설치해 주세요. (참고: 가중치 파일이 없어도 프로그램은 실행되며, 이 경우 온타깃 스캔/GC 기반 mock 점수까지만 동작하고 실제 딥러닝 예측 시에만 프레임워크 설치를 요구합니다.)

```bash
# 1. 필수 라이브러리
pip install numpy

# 2. SpCas9, Cas12a 예측 모델용 (DeepSpCas9, DeepCpf1)
# 텐서플로우(TensorFlow)가 필요합니다.
pip install tensorflow

# 3. Prime Editor 예측 모델용 (PRIDICT)
# 파이토치(PyTorch)가 필요합니다.
pip install torch
```

> **참고**: Cas13d(RNA 타겟팅, Cas13a/Cas13d 계열)는 이 툴에서 지원하지 않습니다. 실존하는 Cas13d 효율 예측 도구(Wessels et al. 2020, Sanjana lab, [gitlab.com/sanjanalab/cas13](https://gitlab.com/sanjanalab/cas13))는 PyTorch 가중치 파일이 아니라 R + RNAhybrid + ViennaRNA 기반의 별도 파이프라인이라, 이 프로젝트의 "가중치 파일 다운로드 후 즉시 스코어링" 구조와 맞지 않아 제외했습니다.
>
> **참고**: SaCas9와 Cas14a(Cas12f)는 가이드 탐색(scanning)은 지원하지만, 딥러닝 효율 예측 모델은 없습니다. 검색해본 결과 이 두 시스템 이름으로 실제 공개된 사전학습 가중치 파일을 찾을 수 없었습니다 (기존 코드의 다운로드 링크는 존재하지 않는 저장소를 가리키고 있었습니다). 대신 결과 테이블에는 GC 함량 기반의 근사(heuristic) 점수가 "Heuristic (GC-based) Score"로 표시됩니다.

### 피드백

툴 관련 불편사항이나 피드백 및 수정사항이 있으면 아래 연락처로 연락부탁드립니다.
* **Email:** kmh40283656@gmail.com

