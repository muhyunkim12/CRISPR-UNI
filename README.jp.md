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

## 実行方法 (How to Run)

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
