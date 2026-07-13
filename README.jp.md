<div align="right">
  <a href="README.md">🇰🇷 한국어</a> | 
  <a href="README.en.md">🇺🇸 English</a> | 
  <a href="README.ja.md">🇯🇵 日本語</a>
</div>

# CRISPR-UNI 
🧬 Next-Gen CRISPR & Prime Editing Multi-Model Design Pipeline (CLI)

本プロジェクトは、多様なゲノム編集ツール（SpCas9、Cas12a、Prime Editorなど）の分子生物学的メカニズムの違いを理解し、各システムに最適化された最新のSOTAディープラーニング予測モデル（DeepSpCas9、PRIDICTなど）を自動的にマッチングさせ、編集効率を最大化する実験室ワークステーション用の高性能CLIパイプラインです。

## コアアーキテクチャの特徴

- **AIモデリング（Factory Pattern）：** ユーザーがゲノム編集ツールを選択した際、単純な条件分岐（if文）で処理するのではなく、オブジェクト指向のファクトリーパターンを採用しました。`PredictorFactory`が、各ツールの特性（DNA切断、RNAターゲティング、逆転写など）に合わせた専用のディープラーニング推論オブジェクトを動的に割り当てます。

- **プライムエディティング（PE）の多次元最適化：** 初期効率が低いプライムエディティングの限界を克服するため、切断位置を基準としたPBS（10, 13, 15nt）とRTテンプレート（12, 15, 20nt）の3×3のpegRNA組み合わせを自動生成し、TransformerベースのPRIDICTアーキテクチャを通じて最終的なスコアリングを実行します。

- **防御的プログラミング（Defensive Programming）：** 入力配列のバリデーション（A, T, G, C以外の例外処理）や、サーバーインフラ内の重み（Weight）ファイルの有無を確認する例外処理ロジックを緻密に実装し、パイプラインの安定性と堅牢性を高めました。

### PREREQUISITE(重要）

本パイプラインのAI予測を正常に動作させるためには、各モデルの事前学習済み重み付け（Weights）ファイルが `weights/` ディレクトリ内に存在している必要があります。

```bash
# 重み付けファイルを保存するディレクトリを作成し、移動します。
mkdir -p weights && cd weights

# 1. DeepSpCas9 の重み付け (SpCas9用)
# 出典: 延世大学 Kim Lab 公式 GitHub
wget -O DeepSpCas9_model.h5 "[https://raw.githubusercontent.com/myungjinkim/DeepSpCas9/master/DeepSpCas9_model.h5](https://raw.githubusercontent.com/myungjinkim/DeepSpCas9/master/DeepSpCas9_model.h5)"

# 2. DeepCpf1 の重み付け (Cas12a用)
# 出典: 延世大学 Kim Lab 公式 GitHub
wget -O DeepCpf1_model.h5 "[https://raw.githubusercontent.com/myungjinkim/DeepCpf1/master/DeepCpf1_model.h5](https://raw.githubusercontent.com/myungjinkim/DeepCpf1/master/DeepCpf1_model.h5)"

# 3. PRIDICT の重み付け (Prime Editor用)
# 出典: PRIDICT 公式 Zenodo アーカイブ (原著者の配布リンク)
# 注意: PRIDICTモデルはファイルサイズが大きいため、ダウンロードに時間がかかる場合があります。
wget -O PRIDICT_model.pt "[https://zenodo.org/record/8208465/files/pridict_v2_model.pt](https://zenodo.org/record/8208465/files/pridict_v2_model.pt)"

# ダウンロードが完了したら、親ディレクトリに戻ります。
cd ..
```

## How to Run

```bash
# 1. 環境セットアップ
pip install -r requirements.txt

# 2. パイプラインの実行
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
システムを選択してください (1-6): 3
ターゲットDNA配列を入力してください: ATCG...
>> [PRIDICT Transformer Model] に基づくpegRNA最適化組み合わせの計算を開始します...
```

###  パッケージの手動インストール 

実験室のサーバー環境によっては、`pip install -r requirements.txt` の実行時に依存関係の競合やエラーが発生する場合があります。その際は、共通の必須ライブラリをインストールした後、使用するCRISPRシステムに合わせてAIフレームワークを個別にインストールしてください。

```bash
# 1. 共通必須ライブラリ (配列処理およびデータ演算)
pip install biopython pandas numpy scikit-learn

# 2. SpCas9, SaCas9, Cas12a 予測モデル用 (DeepSpCas9, DeepCpf1)
# TensorFlow環境が必要です。
pip install tensorflow keras

# 3. Prime Editor (PE) 予測モデル用 (PRIDICT)
# PyTorchおよびTransformers環境が必要です。
pip install torch torchvision torchaudio transformers

# 4. モデル重み付け(Weight)ファイルダウンロード用
pip install gdown
```
