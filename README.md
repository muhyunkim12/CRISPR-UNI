<div align="right">
  <a href="README.md">🇰🇷 한국어</a> | 
  <a href="README.en.md">🇺🇸 English</a> | 
  <a href="README.ja.md">🇯🇵 日本語</a>
</div>

# CRISPR-UNI 
🧬 Next-Gen CRISPR & Prime Editing Multi-Model Design Pipeline (CLI)
본 프로젝트는 다양한 유전자 가위(SpCas9, Cas12a, Prime Editor 등)의 분자생물학적 메커니즘 차이를 이해하고, 각 시스템에 최적화된 최신 SOTA 딥러닝 예측 모델(DeepSpCas9, PRIDICT 등)을 자동으로 매칭하여 교정 효율을 극대화하는 실험실 워크스테이션용 고성능 CLI 파이프라인입니다.

# 핵심 아키텍처 특징

- AI Modeling: 사용자가 유전자 가위를 선택하면 조건문으로 코드를 기우는 대신, Predictor가 가위의 특성(DNA 절단, RNA 타겟팅, 역전사 등)에 맞는 전용 딥러닝 추론 객체를 동적으로 할당하는 객체 지향 디자인 패턴을 적용했습니다.

- Multi-Dimensional Optimization for Prime Editing: 초기 효율이 낮은 프라임 에디팅(PE)의 한계를 극복하기 위해, 절단 위치 기준 PBS(10, 13, 15nt)와 RT 템플릿(12, 15, 20nt)의 3x3 가이드 RNA 조합을 자동 생성하고, 트랜스포머(Transformer) 기반 PRIDICT 아키텍처를 통해 최종 스코어링을 수행합니다.

- Defensive Programming: 입력 서열의 유효성 검사(A, T, G, C 외 예외 처리) 및 서버 인프라 내 가중치 파일 유무를 체크하는 예외 처리 로직을 촘꼼히 반영하여 파이프라인의 안정성을 높였습니다.

# PREREQUISITE (Important)

본 파이프라인의 AI 예측을 정상적으로 구동하기 위해서는 각 모델의 사전 학습된 가중치(Weights) 파일이 `weights/` 디렉토리 내에 존재해야 합니다.

```bash
# 가중치 파일이 저장될 폴더를 만들고 이동합니다.
mkdir -p weights && cd weights

# 1. DeepSpCas9 가중치 (SpCas9 용)
# 출처: Yonsei Univ. Kim Lab Official GitHub
wget -O DeepSpCas9_model.h5 "[https://raw.githubusercontent.com/myungjinkim/DeepSpCas9/master/DeepSpCas9_model.h5](https://raw.githubusercontent.com/myungjinkim/DeepSpCas9/master/DeepSpCas9_model.h5)"

# 2. DeepCpf1 가중치 (Cas12a 용)
# 출처: Yonsei Univ. Kim Lab Official GitHub
wget -O DeepCpf1_model.h5 "[https://raw.githubusercontent.com/myungjinkim/DeepCpf1/master/DeepCpf1_model.h5](https://raw.githubusercontent.com/myungjinkim/DeepCpf1/master/DeepCpf1_model.h5)"

# 3. PRIDICT 가중치 (Prime Editor 용)
# 출처: PRIDICT 공식 Zenodo 아카이브 (원작자 배포 링크)
# 주의: PRIDICT 모델은 용량이 크므로 다운로드에 시간이 걸릴 수 있습니다.
wget -O PRIDICT_model.pt "[https://zenodo.org/record/8208465/files/pridict_v2_model.pt](https://zenodo.org/record/8208465/files/pridict_v2_model.pt)"

# 다운로드가 완료되면 다시 상위 폴더로 빠져나옵니다.
cd ..
```

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
가위 종류를 선택하세요 (1-6): 3
타겟 DNA 서열을 입력하세요: ATCG...
>> [PRIDICT Transformer Model] 기반 pegRNA 최적화 조합 연산을 시작합니다...
```

### 패키지 수동 설치 가이드 

실험실 서버 환경에 따라 `pip install -r requirements.txt` 실행 시 버전 충돌이나 에러가 발생할 수 있습니다. 이 경우, 파이프라인의 공통 필수 패키지를 먼저 설치한 후, 사용하고자 하는 CRISPR 시스템에 맞춰 AI 딥러닝 프레임워크를 개별 설치해 주세요.

```bash
# 1. 공통 필수 라이브러리 (서열 처리 및 데이터 연산)
pip install biopython pandas numpy scikit-learn

# 2. SpCas9, SaCas9, Cas12a 예측 모델용 (DeepSpCas9, DeepCpf1 기반)
# 텐서플로우(TensorFlow) 연산이 필요합니다.
pip install tensorflow keras

# 3. Prime Editor (PE) 예측 모델용 (PRIDICT 기반)
# 파이토치(PyTorch) 및 트랜스포머(Transformers) 연산이 필요합니다.
pip install torch torchvision torchaudio transformers

# 4. 모델 가중치 다운로드 스크립트 실행용
pip install gdown
```
### 피드백

툴 관련 불편사항이나 피드백 및 수정사항이 있으면 아래 연락처로 연락부탁드립니다.
* **Email:** kmh40283656@gmail.com

