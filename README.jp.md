<div align="right">
  <a href="README.md">🇰🇷 한국어</a> | 
  <a href="README.en.md">🇺🇸 English</a> | 
  <a href="README.jp.md">🇯🇵 日本語</a>
</div>

# CRISPR-UNI 
🧬 Next-Gen CRISPR & Prime Editing Multi-Model Design Pipeline (CLI)

本プロジェクトは、多様なゲノム編集ツール（SpCas9、Cas12a、Prime Editorなど）の分子生物学的メカニズムの違いを理解し、各システムに最適化された最新のSOTAディープラーニング予測モデル（DeepSpCas9、PRIDICTなど）を自動的にマッチングさせ、編集効率を最大化する実験室ワークステーション用の高性能CLIパイプラインです。

## コアアーキテクチャの特徴

- **AIモデリング（Factory Pattern）：** ユーザーがゲノム編集ツールを選択した際、単純な条件分岐（if文）で処理するのではなく、オブジェクト指向のファクトリーパターンを採用しました。`PredictorFactory`が、各ツールの特性（DNA切断、RNAターゲティング、逆転写など）に合わせた専用のディープラーニング推論オブジェクトを動的に割り当てます。

- **プライムエディティング（PE）の多次元最適化：** 初期効率が低いプライムエディティングの限界を克服するため、切断位置を基準としたPBS（10, 13, 15nt）とRTテンプレート（12, 15, 20nt）の3×3のpegRNA組み合わせを自動生成し、TransformerベースのPRIDICTアーキテクチャを通じて最終的なスコアリングを実行します。

- **防御的プログラミング（Defensive Programming）：** 入力配列のバリデーション（A, T, G, C以外の例外処理）や、サーバーインフラ内の重み（Weight）ファイルの有無を確認する例外処理ロジックを緻密に実装し、パイプラインの安定性と堅牢性を高めました。

- **オフターゲット分析（任意）：** ローカルに準備した参照ゲノムに対して[Cas-OFFinder](https://github.com/snugel/cas-offinder)を連携させ、各ガイド候補がゲノム内の他の位置とどれだけ類似しているかを検査します。設定方法は下記の「オフターゲット分析」セクションを参照してください。

### PREREQUISITE(重要）

本パイプラインのAI予測を正常に動作させるためには、各モデルの事前学習済み重み付け（Weights）ファイルが `weights/` ディレクトリ内に存在している必要があります。

**重みファイルに加えて、対応するディープラーニングフレームワーク（TensorFlowおよび/またはPyTorch）も別途インストールする必要があります。** `requirements.txt`は容量の都合上これらのフレームワークをデフォルトではインストールせず、使用するシステムに応じてどちらか一方があれば十分です。各システムにどちらのフレームワークが必要かは、下記の「パッケージの手動インストール」セクションを参照してください。

重みファイル自体を取得する一番簡単な方法は、下記のスクリプトを実行することです。ダウンロード元が確認できているモデルは自動的に取得し、コードが想定する正しいファイル名で保存します。

```bash
python download_weights.py
```

手動で取得する場合は下記のコマンドを使用してください。（注意: `-O` の後のファイル名は、配布元のファイル名とは異なり、コードが期待する名前に必ず合わせる必要があります。）

```bash
# 重み付けファイルを保存するディレクトリを作成し、移動します。
mkdir -p weights && cd weights

# 1. DeepSpCas9 の重み付け (SpCas9用)
# 出典: 延世大学 Kim Lab 公式 GitHub
wget -O DeepSpCas9_weights.h5 "https://raw.githubusercontent.com/myungjinkim/DeepSpCas9/master/DeepSpCas9_model.h5"

# 2. DeepCpf1 の重み付け (Cas12a用)
# 出典: 延世大学 Kim Lab 公式 GitHub
wget -O DeepCpf1_weights.h5 "https://raw.githubusercontent.com/myungjinkim/DeepCpf1/master/DeepCpf1_model.h5"

# 3. PRIDICT の重み付け (Prime Editor用)
# 出典: PRIDICT 公式 Zenodo アーカイブ (原著者の配布リンク)
# 注意: PRIDICTモデルはファイルサイズが大きいため、ダウンロードに時間がかかる場合があります。
wget -O PRIDICT_model.pt "https://zenodo.org/record/8208465/files/pridict_v2_model.pt"

# ダウンロードが完了したら、親ディレクトリに戻ります。
cd ..
```

## オフターゲット分析（任意）

デフォルトでは、入力した配列内でのオンターゲットスキャンのみが実行されます。実際の実験デザインのためにゲノム全体でのオフターゲットリスクも確認したい場合は、以下の2つの準備が追加で必要です。

**ステップ1: 参照ゲノムの準備（生物種ごとに1回のみ）**

オフターゲット検索には、対象生物種の全ゲノムFASTAがローカルに必要です。`prepare_genome.py`が[NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/)を通じて自動的に取得します（`conda install -c conda-forge ncbi-datasets-cli`でインストール可能）。

```bash
# 登録されている全生物種(Human, Mouse, Rice, Arabidopsis, Yeast, Ecoli)のゲノムを準備
python prepare_genome.py

# 特定の生物種のみ準備する場合
python prepare_genome.py Human

# ダウンロードせずステータスのみ確認
python prepare_genome.py --list
```

取得したゲノムは`genomes/<生物種>/`フォルダに保存され、`.gitignore`に既に含まれているためコミットされません。

**ステップ2: Cas-OFFinderのインストール**

ゲノム規模のオフターゲットマッチング検索を自前で実装する代わりに、実績のあるオープンソースツール[Cas-OFFinder](https://github.com/snugel/cas-offinder/releases)を使用します。バイナリをダウンロードし、`cas-offinder`コマンドとしてPATHに登録してください。

**ステップ3: `--check-offtargets`オプションで実行**

```bash
python3 main.py --check-offtargets

# 許容ミスマッチ数も調整可能(デフォルト3)
python3 main.py --check-offtargets --max-mismatches 2
```

参照ゲノムまたはCas-OFFinderのいずれかが準備できていない場合でも、プログラムは停止せず案内メッセージを表示し、オンターゲットの結果は通常どおり表示されます。

> **注記**: SpCas9とPrime_Editor(どちらも20ntスペーサー+3' PAM構造)については論文発表済みの実際のCFD(Cutting Frequency Determination、Doench et al. 2016)スコアを、SaCas9(21ntスペーサー+3' PAM)については論文発表済みの実際のSaCas9特異性モデル(Tycko et al. 2018, Nat Commun)をそのまま計算しており、結果テーブルにはそれぞれ「CFD」「SaCas9-specificity」と表示されます。それ以外のシステム(Cas12a、Cas14a)は、検証済みで移植可能なスコア表が確認できていないため、ミスマッチ数のみを反映した単純な近似スコアに自動的に切り替わり「heuristic」と表示されます。この場合、PAMに近いミスマッチほど影響が大きいというシード領域の重み付けは反映されていません。

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
システムを選択してください (1-5): 3
ターゲットDNA配列を入力してください: ATCG...
>> [PRIDICT Transformer Model] に基づくpegRNA最適化組み合わせの計算を開始します...
```

###  パッケージの手動インストール 

`requirements.txt`は必須依存関係(numpy)のみをインストールします。ディープラーニングフレームワーク(TensorFlow、PyTorch)は容量が大きく、実際に使用するCRISPRシステムによってどちらか一方だけあれば良いため、デフォルトのインストールには含めていません。使用したいシステムに合わせて必要なものだけインストールしてください。(注: フレームワークが未インストールでもツール自体は動作し、オンターゲットスキャンとGCベースのモックスコアはそのまま使えます。フレームワークが必要になるのは、実際にディープラーニング予測モデルを呼び出す時だけです。)

```bash
# 1. 必須ライブラリ
pip install numpy

# 2. SpCas9, Cas12a 予測モデル用 (DeepSpCas9, DeepCpf1)
# TensorFlowが必要です。
pip install tensorflow

# 3. Prime Editor 予測モデル用 (PRIDICT)
# PyTorchが必要です。
pip install torch
```

> **注記**: Cas13d(RNAターゲティング)はこのツールではサポートしていません。実在するCas13d効率予測ツール(Wessels et al. 2020, Sanjana lab, [gitlab.com/sanjanalab/cas13](https://gitlab.com/sanjanalab/cas13))はダウンロード可能なPyTorch重みファイルではなく、R + RNAhybrid + ViennaRNAベースの別パイプラインであり、本プロジェクトの「重みファイルをダウンロードして即座にスコアリングする」構造と合わないため、対象から除外しました。
>
> **注記**: SaCas9とCas14a(Cas12f)はガイド探索(scanning)には対応していますが、ディープラーニング効率予測モデルはありません。調査した限り、この2つのシステム名で実在する公開済みの学習済み重みファイルは見つかりませんでした(従来のコードのダウンロードリンクは存在しないリポジトリを指していました)。結果テーブルには代わりにGC含量に基づく近似(heuristic)スコアが「Heuristic (GC-based) Score」として表示されます。

###  フィードバック及びお問い合わせ

ツールの利用に関する不具合、フィードバック、または修正のご要望などがございましたら、以下の連絡先までお気軽にお問い合わせください。

* **Email:** kmh40283656@gmail.com
