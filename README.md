<div align="right">
  <a href="README.md">🇰🇷 한국어</a> | 
  <a href="README.en.md">🇺🇸 English</a> | 
  <a href="README.ja.md">🇯🇵 日本語</a>
</div>

# CRISPR-UNI 
🧬 Next-Gen CRISPR & Prime Editing Multi-Model Design Pipeline (CLI)
본 프로젝트는 다양한 유전자 가위(SpCas9, Cas12a, Prime Editor 등)의 분자생물학적 메커니즘 차이를 이해하고, 각 시스템에 최적화된 최신 SOTA 딥러닝 예측 모델(DeepSpCas9, PRIDICT 등)을 자동으로 매칭하여 교정 효율을 극대화하는 실험실 워크스테이션용 고성능 CLI 파이프라인입니다.

# 핵심 아키텍처 특징

- AI Modeling: 사용자가 유전자 가위를 선택하면 조건문으로 코드를 기우는 대신, PredictorFactory가 가위의 특성(DNA 절단, RNA 타겟팅, 역전사 등)에 맞는 전용 딥러닝 추론 객체를 동적으로 할당하는 객체 지향 디자인 패턴을 적용했습니다.

- Multi-Dimensional Optimization for Prime Editing: 초기 효율이 낮은 프라임 에디팅(PE)의 한계를 극복하기 위해, 절단 위치 기준 PBS(10, 13, 15nt)와 RT 템플릿(12, 15, 20nt)의 3x3 가이드 RNA 조합을 자동 생성하고, 트랜스포머(Transformer) 기반 PRIDICT 아키텍처를 통해 최종 스코어링을 수행합니다.

- Defensive Programming: 입력 서열의 유효성 검사(A, T, G, C 외 예외 처리) 및 서버 인프라 내 가중치 파일 유무를 체크하는 예외 처리 로직을 촘꼼히 반영하여 파이프라인의 안정성을 높였습니다.

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
