# ヒト海馬 SGZ 加齢 Secretome 解析 — 実行指示書 完全版 v4

**プロジェクト名**: `hSGZ-aging-secretome`
**ハードウェア**: Ubuntu 22.04 / Core i9-14900 / RTX 4090 (24 GB VRAM) / **RAM 128 GB** / SSD 4 TB / HDD 20 TB
**支援AI**: Claude Code（対話的試行錯誤）、Codex CLI（定型処理・scaffolding）
**版**: v4 — v3 完全版に GSE160189 (regional axis reference) と GSE264692 (Tran 2025 Visium 本体) を追加。Figure 5 を二重 Visium + FISH で再構成し、regional specificity panel を新設

## v3 → v4 主要変更

1. **GSE160189 追加** (anterior/posterior axis snRNA-seq reference) → Figure 1 の細胞型 annotation reference に組み込み、Figure 5 で regional specificity panel を新設
2. **GSE264692 追加** (Tran 2025 Visium 本体) → v3 では GSE264624 単独で空間検証を記述していたが、Visium データの実体は GSE264692 側にある可能性が高く、これがないと cell2location が走らない重大な見落としを訂正
3. Section 1.1 / 1.2 / 3.2 / 4.5 / 8.2-8.4 を更新

---

## 目次

- [Section 0: プロジェクト全体方針](#section-0-プロジェクト全体方針)
- [Section 1: データセット戦略と差別化](#section-1-データセット戦略と差別化)
- [Section 2: Phase 0 — 環境構築](#section-2-phase-0--環境構築)
- [Section 3: Phase 1 — データ取得](#section-3-phase-1--データ取得)
- [Section 4: Figure 1 — 細胞アトラスと組成変化](#section-4-figure-1--細胞アトラスと組成変化)
- [Section 5: Figure 2 — Neurogenic lineage の精緻化](#section-5-figure-2--neurogenic-lineage-の精緻化)
- [Section 6: Figure 3 — Secretome DE と meta-analysis](#section-6-figure-3--secretome-de-と-meta-analysis)
- [Section 7: Figure 4 — Cell-cell communication + Regulatory](#section-7-figure-4--cell-cell-communication--regulatory)
- [Section 8: Figure 5 — 二重空間検証](#section-8-figure-5--二重空間検証)
- [Section 9: Figure 6 — クロス種比較](#section-9-figure-6--クロス種比較)
- [Section 10: Figure 7 — 統合解釈](#section-10-figure-7--統合解釈)
- [Section 11: 論文化フェーズ](#section-11-論文化フェーズ)
- [Section 12: 補遺 — AI 連携ワークフロー](#section-12-補遺--ai-連携ワークフロー)
- [Section 13: 査読対策の論理武装](#section-13-査読対策の論理武装)
- [Section 14: タイムライン & kickoff](#section-14-タイムライン--kickoff)

---

# Section 0: プロジェクト全体方針

## 0.1 設計思想

1. **Figure 駆動**: 各 phase の goal は「論文 figure 1 枚」に対応。それを再現できるスクリプトと中間データが残ることを成功条件とする。
2. **段階的 commit**: 1 figure 到達ごとに `git tag figureN-v1` を打つ。やり直しが容易になる。
3. **再現性 3 層**: ① conda environment.yml, ② numbered scripts または snakemake, ③ random seed 固定 (`np.random.seed(42)`, `set.seed(42)`)。
4. **データ非編集原則**: 生データ (raw/) には絶対に書き込まない。派生ファイルは processed/, results/ で別ディレクトリ。
5. **バックアップ自動化**: cron で 1 日 1 回 SSD → HDD rsync。Figure 完成時は手動 snapshot。
6. **計画より実データを優先**: Figure 1 で見える biology が Figure 2 以降の analytic decision を決める。事前計画で全体を縛らず、毎 figure 後に `DECISIONS.md` を更新して計画を修正する。

## 0.2 Claude Code / Codex CLI の使い分け

| AI | 適した場面 |
|---|---|
| **Claude Code** | データを見て次の判断をする場面: QC 閾値の決定、UMAP を見てクラスタ再分割、CellChat 出力の解釈、エラーデバッグなど。対話的で文脈が長い作業。 |
| **Codex CLI** | 決まった処理を一気に通す場面: 定型 script 量産、numbered script の scaffolding、繰り返し plot 生成、snakemake rule 追加など。 |

**共通ルール**: プロンプト前に `docs/CONTEXT.md` を最新化（テンプレートは Section 12）。

---

# Section 1: データセット戦略と差別化

## 1.1 確定データセット一覧

| GEO accession | 論文 | モダリティ | サンプル特性 | 役割 |
|---|---|---|---|---|
| **GSE268609** | Lazarov et al. 2024 "A roadmap to human hippocampal neurogenesis in adulthood, aging and AD" | snRNA-seq + **snATAC-seq (multiome)** | 若年成人, 認知正常高齢, MCI, AD | **★Primary discovery** |
| **GSE185277** | Zhou series — lifespan | snRNA-seq | 乳児〜高齢 | imGC 年齢勾配 |
| **GSE185553** | Zhou series — adult | snRNA-seq | 成人 | imGC marker 再現性 |
| **GSE198323** | Zhou et al. 2022 (Nature) | snRNA-seq | AD vs control | aging vs AD 分離 |
| **GSE186538** | Franjic et al. 2022 (Neuron) | snRNA-seq | ヒト・マカク・ブタ HPC + EC | **反証データ** |
| **GSE163737** | Aged macaque + human DG/HPC | snRNA-seq | ヒト 67/85/87/92 歳 + マカク連続 | 高齢 anchor + クロス種 |
| **GSE160189** | Human HPC anterior/posterior axis | snRNA-seq | 成人ヒト long-axis | **領域参照アトラス** + regional specificity 検証 |
| **GSE264624** | Tran et al. 2025 (Nat Neurosci) snRNA 部分 | snRNA-seq | 成人 10 例 前部 HPC | 空間検証用 reference snRNA |
| **GSE264692** | Tran et al. 2025 Visium 部分 (paired with GSE264624) | **Visium SRT** | 成人 10 例 前部 HPC | **空間検証 (1) — Visium 本体** |
| **GSE248545** | Spatial transcriptomics + FISH | spatial + FISH | 成人 DG | 空間検証 (2) — 直交検証 |
| **GSE199243** | Glial diversity in human HPC across lifespan and AD | snRNA-seq | lifespan, AD | ニッチ細胞補強 |
| **GSE325391** *要確認* | Immature neurons in aged human HPC, AD pathology, resilience | snRNA-seq | 高齢健常 + AD + resilience | **最新検証コホート** |
| (補) GSE95753 | Hochgerner 2018 マウス DG | scRNA-seq | mouse | クロス種 (mouse) |
| (補) SCP90 | Habib 2017 DroNc-seq | snRNA-seq | 4 GTEx donor | 古典 reference |

## 1.2 戦略マップ

```
            ┌─────────────────────────────────┐
            │   GSE268609 (Lazarov, multiome) │  ← Discovery
            │   NSC/neuroblast/imGC + niche   │
            └─────────────┬───────────────────┘
                          │ secretome DEGs を抽出
        ┌─────────────────┼─────────────────────┐
        ▼                 ▼                     ▼
  ┌──────────────┐  ┌──────────────┐    ┌──────────────┐
  │ Replication  │  │ Counter      │    │ Aged-specific│
  │ GSE185277/553│  │ GSE186538    │    │ GSE163737    │
  │ /198323      │  │ (Franjic)    │    │ GSE325391    │
  └──────┬───────┘  └──────┬───────┘    └──────┬───────┘
         │                 │                    │
         └────── meta-analysis (metafor) ───────┘
                          │
                          ▼
        ┌─────────────────────────────────────┐
        │  Regional reference (snRNA-seq)     │
        │  GSE160189 (anterior/posterior axis)│
        │  → SGZ 特異性 vs pan-HPC を判定     │
        └─────────────────┬───────────────────┘
                          │
                          ▼
        ┌─────────────────────────────────────┐
        │  Spatial validation (二重)          │
        │  GSE264624+GSE264692 (Visium pair)  │
        │  GSE248545 (spatial + FISH)         │
        └─────────────────┬───────────────────┘
                          │
                          ▼
        ┌─────────────────────────────────────┐
        │  Niche context: GSE199243 glia      │
        │  Cross-species: GSE163737 macaque   │
        │  + Hochgerner 2018 mouse            │
        └─────────────────────────────────────┘
```

## 1.3 Lazarov et al. 2024 (GSE268609) との差別化（最重要）

**先行研究**:
- 同一データ GSE268609 で multiomic 解析済
- CellChat と NeuronChat で aging/MCI/AD のリガンド・レセプター推定済
- 主要結論: NSC で最早期に老化シグナル出現、AD で immature neuron の翻訳機構が広範に停止
- 焦点: シナプス接着分子、神経伝達物質

**本研究の独自性 5 軸 + v4 で追加した第 6 軸**:

1. **Secretome 特化レンズ**: HPA × UniProt × MatrisomeDB の三重定義で前フィルタ。Lazarov は全シグナリング、本研究は**分泌タンパク質のみ**に絞り、cytokine/growth factor/matrisome/neuropeptide のカテゴリ別に解像。
2. **多データセット meta-analysis**: Lazarov は 1 cohort のみ。本研究は **5+ cohort** で metafor による effect size 統合 → reproducibility の数値化。
3. **二重空間検証**: GSE264624+GSE264692 (Visium) + GSE248545 (spatial+FISH) のクロス検証 → 単一プラットフォームのアーティファクトを排除。
4. **クロス種 ortholog 限定比較**: GSE163737 + Hochgerner 2018 で「**ヒト特異的に加齢で変化する secretome**」を抽出。
5. **snATAC × secretome 統合**: secretome 遺伝子の cis 調節領域 (peak) と TF motif 活性が加齢でどう変化するか → Lazarov は遺伝子発現中心で未カバー。
6. **Regional specificity (v4 新設)**: GSE160189 (anterior/posterior axis) と Visium 空間データを組み合わせ、aging secretome シグナルが **SGZ/DG 特異**か **pan-hippocampal** かを判定。「SGZ niche」というテーマの正当性を直接補強し、reviewer の典型質問「あなたの finding は SGZ 特有なのか」に明確に答える防御力。

**論文 Discussion の位置付け**: 「Lazarov 2024 と本研究は相補的: 彼らは neurogenic lineage を receiver 視点で AD 移行を、本研究は niche 細胞からの secretome を sender 視点で normal aging を見る」と整理する。

---

# Section 2: Phase 0 — 環境構築

## 2.1 OS と基本ツール

```bash
sudo apt update && sudo apt upgrade -y

sudo apt install -y build-essential git curl wget htop tmux \
    zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
    libxml2-dev libfontconfig1-dev libcairo2-dev libxt-dev \
    libhdf5-dev libgsl-dev libfftw3-dev pkg-config \
    aria2 pigz parallel rsync

# RAM 確認 (128 GB 想定)
free -h
```

128 GB RAM の場合、Bartosovic を含むような大規模統合でも全 in-memory 動作可。`backed='r'` モードは原則不要。

## 2.2 NVIDIA driver / CUDA

```bash
ubuntu-drivers devices
sudo ubuntu-drivers autoinstall
sudo reboot

# 確認
nvidia-smi   # CUDA 12.x が見える

# CUDA Toolkit はミニマルに。PyTorch 同梱の CUDA を使う方針で OK。
```

## 2.3 mamba と Python 環境

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3
echo 'source $HOME/miniforge3/etc/profile.d/conda.sh' >> ~/.bashrc
echo 'source $HOME/miniforge3/etc/profile.d/mamba.sh' >> ~/.bashrc
source ~/.bashrc
```

メイン Python 環境 (`env_python.yml`):

```yaml
name: sgz
channels: [rapidsai, nvidia, conda-forge, bioconda]
dependencies:
  - python=3.11
  - pip
  # core single-cell
  - scanpy>=1.10
  - anndata>=0.10
  - mudata
  - muon                      # multiome (RNA + ATAC)
  - scvi-tools>=1.1           # GPU integration
  - scvelo
  - cellrank>=2.0
  - palantir
  - scrublet
  - harmonypy
  - decoupler-py
  - pydeseq2
  - sccoda
  - pertpy
  # spatial
  - squidpy
  - spatialdata
  - spatialdata-io
  - cell2location             # GPU
  # multi-omics regulatory
  - snapatac2
  # plotting / utils
  - matplotlib
  - seaborn
  - plotly
  - pyarrow
  - polars
  - jupyterlab
  - pytables
  # PyTorch CUDA 12.1
  - pytorch>=2.2
  - pytorch-cuda=12.1
  # rapids-singlecell (GPU PCA/UMAP/Leiden)
  - rapids-singlecell
  - cuml
  - cudf
  - cugraph
  # workflow
  - snakemake>=8.0
  - dvc
  - pre-commit
  # GEO 取得
  - geoparse
  - pip:
      - cellbender            # GPU でアンビエント RNA 除去
      - solo-sc               # doublet detection
      - episcanpy             # snATAC
      - pycistopic            # multiome regulatory
      - multinichenetr        # 後で R から呼ぶ
```

```bash
mamba env create -f env_python.yml
mamba activate sgz
python -c "import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))"
# True, NVIDIA GeForce RTX 4090
```

## 2.4 R 環境

`env_r.yml`:

```yaml
name: sgz_r
channels: [conda-forge, bioconda]
dependencies:
  - r-base=4.4
  - r-essentials
  - r-seurat>=5.1
  - r-signac
  - r-harmony
  - r-presto
  - r-tidyverse
  - r-patchwork
  - r-ggrepel
  - bioconductor-deseq2
  - bioconductor-edger
  - bioconductor-limma
  - bioconductor-muscat        # pseudobulk DE
  - bioconductor-glmgampoi
  - bioconductor-singlecellexperiment
  - bioconductor-scran
  - bioconductor-scater
  - bioconductor-mofa2
  - bioconductor-jaspar2020
  - bioconductor-tfbstools
  - bioconductor-bsgenome.hsapiens.ucsc.hg38
  - bioconductor-spatiallibd
  - bioconductor-zellkonverter
  - r-seuratdisk
  - r-metafor
  - r-devtools
  - r-biocmanager
```

```bash
mamba env create -f env_r.yml
mamba activate sgz_r

Rscript -e 'devtools::install_github("jinworks/CellChat")'
Rscript -e 'devtools::install_github("saeyslab/nichenetr")'
Rscript -e 'devtools::install_github("saeyslab/multinichenetr")'
Rscript -e 'BiocManager::install("speckle")'
Rscript -e 'devtools::install_github("immunogenomics/presto")'
```

## 2.5 ディレクトリ構造（厳守）

```bash
export PROJ=$HOME/projects/hSGZ-aging-secretome
mkdir -p $PROJ/{raw,processed,results,figures,manuscript,scripts,envs,logs,refs,docs}
cd $PROJ
git init
```

```
hSGZ-aging-secretome/
├── raw/                 # read-only (絶対編集しない)
│   ├── GSE268609_lazarov2024/      # ★Discovery (multiome)
│   ├── GSE185277_zhou_lifespan/
│   ├── GSE185553_zhou_adult/
│   ├── GSE198323_zhou_AD/
│   ├── GSE186538_franjic2022/      # 反証
│   ├── GSE163737_aged_macaque_human/
│   ├── GSE160189_HPC_AP_axis/      # ★v4 追加 (regional reference)
│   ├── GSE264624_tran2025_snrna/   # Tran 2025 snRNA 部分
│   ├── GSE264692_tran2025_visium/  # ★v4 追加 (Tran 2025 Visium 本体)
│   ├── GSE248545_spatial_FISH/     # 空間 (2)
│   ├── GSE199243_glia_lifespan/
│   ├── GSE325391_immature_AD_resilience/  # 要確認
│   ├── GSE95753_mouse_DG/
│   └── refs/
├── processed/
│   ├── per_dataset/     # 各データセット個別 QC 後
│   ├── integrated/      # 統合済 anndata / Seurat
│   └── pseudobulk/      # 集計済 count
├── results/
│   ├── figure1/
│   ├── figure2/
│   └── ...
├── figures/             # PDF/SVG/PNG。Figure ごとサブフォルダ
├── scripts/
│   ├── 00_env/
│   ├── 01_download/
│   ├── 02_qc/
│   ├── 03_integrate/
│   ├── 04_figure1/
│   └── ...
├── envs/                # env_*.yml の snapshot
├── refs/                # secretome DB, ortholog table, marker gene
├── logs/
├── docs/
│   ├── CONTEXT.md       # AI に渡す現状サマリ
│   ├── DECISIONS.md     # QC 閾値・統合パラメタ等の判断ログ
│   └── DATA_PROVENANCE.md  # 出所・ライセンス
└── .gitignore
```

`.gitignore`:

```
raw/
processed/
results/**/*.h5ad
results/**/*.h5
figures/*.png
figures/*.pdf
.Rhistory
.RData
__pycache__/
.ipynb_checkpoints/
```

`processed/` の重要オブジェクトは DVC で版管理:

```bash
dvc init
dvc remote add -d hdd_backup /mnt/hdd_20tb/dvc_store
```

## 2.6 バックアップ自動化（最重要）

```bash
# crontab -e で追加
0 3 * * * rsync -aAX --delete --exclude='.snakemake' --exclude='.dvc/cache' \
    $HOME/projects/hSGZ-aging-secretome/ /mnt/hdd_20tb/backup/hSGZ-aging-secretome/ \
    >> $HOME/projects/hSGZ-aging-secretome/logs/rsync.log 2>&1
```

Figure 完成時の手動 snapshot:

```bash
# scripts/00_env/snapshot.sh
#!/bin/bash
DATE=$(date +%Y%m%d_%H%M)
TAG=$1   # 例: figure1-v1
rsync -aAX $PROJ/ /mnt/hdd_20tb/snapshots/${TAG}_${DATE}/
git tag $TAG
echo "Snapshot saved: ${TAG}_${DATE}"
```

---

# Section 3: Phase 1 — データ取得

## 3.1 ダウンロード基本テンプレート

```bash
# scripts/01_download/lib.sh
download_geo_suppl() {
    local GSE=$1
    local OUTDIR=$2
    mkdir -p "$OUTDIR" && cd "$OUTDIR"
    local PREFIX=$(echo $GSE | sed 's/...$/nnn/')
    # GEO supplementary を再帰取得
    wget -r -np -nH --cut-dirs=4 -R "index.html*" \
        "https://ftp.ncbi.nlm.nih.gov/geo/series/${PREFIX}/${GSE}/suppl/"
    cd -
}
```

## 3.2 各データセット取得

### 3.2.1 GSE268609 — Lazarov 2024 multiome (★Discovery)

```bash
source scripts/01_download/lib.sh
download_geo_suppl GSE268609 $PROJ/raw/GSE268609_lazarov2024
```

通常 multiome は以下を含む:
- `filtered_feature_bc_matrix.h5` (RNA + ATAC)
- `atac_fragments.tsv.gz` + `.tbi`
- `per_sample_metadata.csv`
- 著者公開の Seurat object (RDS) があれば最速

**Disk**: ATAC fragment が大きく 1 sample あたり数 GB。総量 30-100 GB。

著者公開コードがあれば取得:

```bash
# Lazarov 論文の Code Availability を確認
# 例: git clone https://github.com/lazarov-lab/<repo> $PROJ/refs/lazarov2024_code
```

### 3.2.2 Zhou シリーズ — GSE185277 / GSE185553 / GSE198323

```bash
for GSE in GSE185277 GSE185553 GSE198323; do
    download_geo_suppl $GSE $PROJ/raw/${GSE}_zhou
done
```

Zhou et al. 2022 Nature の GitHub (`songlab-cam` 系) には **imGC 同定 ML モデル**が公開されている可能性が高い:

```bash
# 著者の repository を要確認後
# git clone https://github.com/<authors>/<repo> $PROJ/refs/zhou_imGC_classifier
```

### 3.2.3 GSE186538 — Franjic 2022 (反証)

```bash
download_geo_suppl GSE186538 $PROJ/raw/GSE186538_franjic2022

# Sestan lab portal から処理済 Seurat object も取れる可能性
# http://resources.sestanlab.org/hippocampus/
```

### 3.2.4 GSE163737 — Aged macaque + human

```bash
download_geo_suppl GSE163737 $PROJ/raw/GSE163737_aged
```

ヒト 4 サンプル (67, 85, 87, 92 歳)。マカクは年齢勾配が連続的にあり、種内 aging trajectory のクロス種照合に最適。

### 3.2.5 GSE160189 — Anterior/posterior axis reference (v4 追加)

```bash
download_geo_suppl GSE160189 $PROJ/raw/GSE160189_HPC_AP_axis
```

ヒト海馬の long-axis (anterior–posterior) snRNA-seq atlas。**年齢比較に直接使うわけではないが**、以下の用途で本研究の質を上げる:

1. **細胞型 annotation の reference**: scANVI/scarches で label transfer の教師に使用。Lazarov の著者ラベルと突き合わせて annotation 品質を二重チェック。
2. **Regional specificity 解析の起点**: anterior vs posterior axis での baseline 細胞型分布を学習し、後の Figure 5 で aging secretome が axis 全体で一貫するか、特定セグメントに偏るかを判定。
3. **DG vs CA vs SUB の subfield 分類器**: snRNA 単独で region prediction model を学習し、Discovery (GSE268609) の各細胞に "予測 region" ラベルを付与可能。これによりニッチ細胞 (Astro, Micro) の DG 局在度合いが測れる。

### 3.2.6 GSE264624 + GSE264692 — Tran 2025 paired (snRNA + Visium)

**v4 重要訂正**: Tran 2025 (Nature Neurosci) は複数 GEO accession に分けて deposit されている。snRNA-seq 部分が GSE264624、Visium 空間データ本体が GSE264692 という構成。**両方を取得しないと Figure 5 の cell2location デコンが走らない**。

```bash
# snRNA-seq 部分
download_geo_suppl GSE264624 $PROJ/raw/GSE264624_tran2025_snrna

# Visium 部分 (空間検証の本体)
download_geo_suppl GSE264692 $PROJ/raw/GSE264692_tran2025_visium

# spatialLIBD R パッケージで配布されている可能性も確認
mamba activate sgz_r
Rscript -e '
  library(spatialLIBD)
  # 当該 dataset の fetch_data() 引数を vignette で確認
'
```

**ペアリング処理**: GSE264624 と GSE264692 のサンプルは同じドナー由来である可能性が高い。`metadata_master.csv` に `paired_with` カラムを追加し、snRNA と Visium のペア関係を明示記録。これにより cell2location の reference signature を Tran 2025 自身の snRNA で構築 → 同じドナーの Visium にデコンする self-consistency 確認も可能。

10 ドナー × 36 capture area、~150,000 spots。

### 3.2.7 GSE248545 — Spatial + FISH (直交検証)

```bash
download_geo_suppl GSE248545 $PROJ/raw/GSE248545_spatial_FISH
```

このデータは「ヒト DG では proliferative gene 発現細胞が極めて少なく、DCX+ 細胞は抑制性ニューロン様」という主張のソース。Figure 5 で Tran 2025 (GSE264692) と並置することで、SGZ secretome シグナルの空間整合性を相互検証。

### 3.2.8 GSE199243 — Glial atlas

```bash
download_geo_suppl GSE199243 $PROJ/raw/GSE199243_glia_lifespan
```

ニッチ secretome の主要 sender (astrocyte/microglia) のサブタイプ分類を本データで学習し、Discovery に label transfer すると粒度が上がる。

### 3.2.9 GSE325391 — 最新候補（要確認）

```bash
mamba activate sgz
python -c "
import urllib.request
url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325391'
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
print(urllib.request.urlopen(req).read().decode()[:5000])
"
# Status: Public ならダウンロード可、Private/embargoed なら待機
```

ダウンロード可能なら最新の AD + cognitive resilience コホートとして Figure 3 の meta-analysis に追加。

### 3.2.10 補助 reference

```bash
# Hochgerner 2018 マウス DG
download_geo_suppl GSE95753 $PROJ/raw/GSE95753_mouse_DG

# Habib 2017 DroNc-seq
# Single Cell Portal SCP90 で公開
# https://singlecell.broadinstitute.org/single_cell/study/SCP90
```

## 3.3 Reference databases

```bash
cd $PROJ/refs

# Human Protein Atlas - secretome
wget -O hpa_secretome.tsv \
    "https://www.proteinatlas.org/api/search_download.php?search=protein_class:Predicted+secreted+proteins&format=tsv&columns=g,gs,sc,upbp"

# UniProt - Secreted ヒト
wget -O uniprot_secreted_human.tsv \
    "https://rest.uniprot.org/uniprotkb/stream?query=organism_id%3A9606+AND+%28cc_subcellular_location%3Asecreted%29&format=tsv&fields=accession,gene_names,protein_name,sequence"

# MatrisomeDB (browser DL)
# http://matrisomedb.org/  → "Human Matrisome List" (Excel) を refs/ に保存
# refs/Hs_Matrisome_Masterlist.xlsx として配置

# NicheNet v2 priors
mkdir -p nichenet
wget -P nichenet/ "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"
wget -P nichenet/ "https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"
wget -P nichenet/ "https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"

# Ortholog (human-macaque-mouse 1:1) - BioMart
# scripts/01_download/get_orthologs_biomart.py で取得（再現性のためスクリプト化）

# SenMayo (senescence) gene set - Cell Stem Cell 2022 supplement (browser DL)
# SASP atlas - manual download
# CellChatDB は CellChat::CellChatDB.human で内蔵されている
```

## 3.4 メタデータ整備（横断比較に必須）

`refs/metadata_master.csv` 推奨カラム:

| 列 | 例 | 必須 |
|---|---|---|
| dataset | GSE268609 | ◯ |
| sample_id | GSM_XXXXXXX | ◯ |
| donor_id | D001 | ◯ |
| age_years | 67 | ◯ |
| age_group | 60-80 | ◯ |
| sex | M / F | ◯ |
| pmi_hours | 12.5 | ◯ |
| rin | 7.2 | △ |
| brain_region | anterior HPC, DG | ◯ |
| diagnosis | neurotypical / AD / MCI / epilepsy / resilient | ◯ |
| modality | snRNA / snATAC / multiome / Visium / FISH | ◯ |
| species | human / macaque / mouse | ◯ |
| platform | 10x v3.1 / 10x multiome / Drop-seq / Visium / Stereo-seq | ◯ |
| paired_with | GSM_YYYYYYY (snRNA ↔ Visium ペア) | △ |
| n_cells_raw | 12345 |  |
| n_cells_post_qc | 9876 |  |
| download_date | 2026-05-09 | ◯ |
| source_url | … | ◯ |
| reference_pmid | 35794479 | ◯ |

GEO metadata の自動抽出:

```python
# scripts/01_download/parse_geo_metadata.py
import GEOparse, pandas as pd

datasets = ['GSE268609','GSE185277','GSE185553','GSE198323','GSE186538',
            'GSE163737','GSE264624','GSE248545','GSE199243','GSE325391']

rows = []
for gse_id in datasets:
    try:
        gse = GEOparse.get_GEO(geo=gse_id, destdir='raw/refs/')
        for gsm_name, gsm in gse.gsms.items():
            rows.append({
                'dataset': gse_id,
                'sample_id': gsm_name,
                'title': gsm.metadata.get('title', [''])[0],
                'characteristics': '|'.join(gsm.metadata.get('characteristics_ch1', [])),
            })
    except Exception as e:
        print(f"Failed {gse_id}: {e}")
pd.DataFrame(rows).to_csv('refs/metadata_raw.csv', index=False)
```

抽出後、**手動レビューで age/sex/diagnosis を統一形式に正規化**（雑にすると全解析が壊れる重要工程）。

## 3.5 一次インベントリ

```python
# scripts/01_download/inventory.py
import hashlib, os, json, pathlib, pandas as pd
from datetime import datetime

ROOT = pathlib.Path(os.environ['PROJ']) / 'raw'
out = []
for ds in ROOT.iterdir():
    if not ds.is_dir(): continue
    for f in ds.rglob('*'):
        if f.is_file():
            with open(f, 'rb') as fh:
                h = hashlib.md5()
                for chunk in iter(lambda: fh.read(2**20), b''): h.update(chunk)
            out.append({'dataset': ds.name, 'path': str(f.relative_to(ROOT)),
                        'size_gb': f.stat().st_size / 1e9, 'md5': h.hexdigest()})
pd.DataFrame(out).to_csv(ROOT.parent / 'refs' / 'data_inventory.csv', index=False)
print(f'{datetime.now()} done; {len(out)} files')
```

## 3.6 Phase 1 完了条件

- [ ] 全 11 データセット (GSE325391 含めれば 12) の md5 が `data_inventory.csv` に記録
- [ ] `metadata_master.csv` に 100+ サンプル分の行が完成 (GSE264624 ↔ GSE264692 の `paired_with` 関係も記録)
- [ ] 各データセットの最小サンプル 1 件で AnnData / SeuratObject load 成功
- [ ] HDD バックアップ完了

→ `git tag phase1-complete` & snapshot.

---

# Section 4: Figure 1 — 細胞アトラスと組成変化

## 4.1 目的

GSE268609 (multiome, normal aging サブセット) を主軸に、海馬の主要細胞型 (excitatory/inhibitory neurons, astrocytes, microglia, OPC, oligodendrocytes, endothelial, VLMC, neurogenic lineage) が**加齢でどう量的に変動するか**を、再現性ある統一細胞型 annotation で示す。後続 figure 全ての出発点。

**戦略上の重要判断**: GSE268609 には AD/MCI も含まれるが、本研究は **normal aging に絞り**、AD/MCI は sensitivity analysis として別 figure (supplementary) で扱う。これにより Lazarov との論点の重ね合わせを避ける。

## 4.2 想定 panel

- **A**: workflow 模式図 (BioRender 別途、Discovery + 4 validation cohort の関係)
- **B**: GSE268609 multiome の integrated UMAP (snRNA 軸)
- **C**: snATAC peak の UMAP (同細胞)、両モダリティの一致を可視化
- **D**: 細胞型 marker dot plot (snRNA gene + snATAC gene activity)
- **E**: 年齢別細胞組成変化 (normal aging のみ)
- **F**: scCODA / propeller による組成変動の統計検定

## 4.3 multiome QC

```python
# scripts/04_figure1/01_multiome_load_qc.py
import scanpy as sc, anndata as ad, muon as mu
import numpy as np, pandas as pd
import scrublet as scr

mdata = mu.read_10x_h5('raw/GSE268609_lazarov2024/<sample>/filtered_feature_bc_matrix.h5')
# mdata.mod['rna'] と mdata.mod['atac'] が自動で分離される

# サンプル ID 付与
mdata.obs['dataset'] = 'GSE268609'
mdata.obs['sample_id'] = '<sample>'
meta = pd.read_csv('refs/metadata_master.csv').query('sample_id == "<sample>"').iloc[0]
for k in ['donor_id','age_years','age_group','sex','pmi_hours','diagnosis']:
    mdata.obs[k] = meta[k]

# RNA 側 QC
rna = mdata.mod['rna']
rna.var['mt'] = rna.var_names.str.startswith('MT-')
rna.var['ribo'] = rna.var_names.str.startswith(('RPS','RPL'))
sc.pp.calculate_qc_metrics(rna, qc_vars=['mt','ribo'], inplace=True, percent_top=None)
rna = rna[(rna.obs.n_genes_by_counts > 500) &
          (rna.obs.n_genes_by_counts < 10000) &
          (rna.obs.pct_counts_mt < 5)].copy()

# ATAC 側 QC (snapatac2 で nucleosome signal / TSSE)
atac = mdata.mod['atac']
sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=True)
# import snapatac2 as snap
# snap.pp.import_data('atac_fragments.tsv.gz', ...)

# 同じバーコードのみ保持
common = rna.obs_names.intersection(atac.obs_names)
rna = rna[common]; atac = atac[common]

# Doublet (RNA-based)
sd = scr.Scrublet(rna.X, expected_doublet_rate=0.06)
rna.obs['doublet_score'], rna.obs['predicted_doublet'] = sd.scrub_doublets()
keep = ~rna.obs['predicted_doublet']
rna = rna[keep]; atac = atac[keep]

# CellBender はオプションで raw_feature_bc_matrix から GPU で別途実施
# cellbender remove-background --input raw.h5 --output cb.h5 --cuda --epochs 150

mdata = mu.MuData({'rna': rna, 'atac': atac})
mdata.write('processed/per_dataset/GSE268609_<sample>_qc.h5mu')
```

**重要**: snRNA-seq では `pct_counts_mt` を 5% 以下に。Zhou など別データは protocol 差で閾値再調整必要。**閾値は dataset ごとに violin で確認してから決め、`docs/DECISIONS.md` に必ず記録**。

## 4.4 統合戦略

選択肢:
- **A: RNA で scVI 統合 → ATAC を knn transfer**: シンプル。本研究主軸は RNA secretome なので推奨。
- **B: MOFA2 / multiVI で multiomic 統合**: 細胞型 annotation が一段精緻。RTX 4090 + 128 GB なら現実的。

**推奨: A 軸で進め、Figure 4 の regulatory 解析対象細胞型に限り B も併走**。

```python
# scripts/04_figure1/02_integrate_RNA.py
import scanpy as sc, scvi, muon as mu, anndata as ad

samples = [...]  # GSE268609 全 normal aging samples
mdatas = [mu.read(f'processed/per_dataset/GSE268609_{s}_qc.h5mu') for s in samples]
rnas = [m.mod['rna'] for m in mdatas]
adata = ad.concat(rnas, join='inner', label='sample_id', keys=samples)
adata.var_names_make_unique()

# AD/MCI を除外 → normal aging 専用 anndata
adata_normal = adata[adata.obs['diagnosis'].isin(['Control','Normal'])].copy()
adata.write('processed/integrated/GSE268609_all.h5ad')      # AD/MCI 含む
adata_normal.write('processed/integrated/GSE268609_normal.h5ad')

# HVG (dataset 横断)
sc.pp.highly_variable_genes(adata_normal, n_top_genes=3000, flavor='seurat_v3',
                            batch_key='sample_id', subset=True)

# scVI on GPU
scvi.model.SCVI.setup_anndata(adata_normal, batch_key='sample_id',
                              categorical_covariate_keys=['sex'],
                              continuous_covariate_keys=['pmi_hours','pct_counts_mt'])
m = scvi.model.SCVI(adata_normal, n_layers=2, n_latent=30, gene_likelihood='nb')
m.train(max_epochs=400, accelerator='gpu', devices=1, batch_size=512,
        early_stopping=True, train_size=0.9)
adata_normal.obsm['X_scVI'] = m.get_latent_representation()

sc.pp.neighbors(adata_normal, use_rep='X_scVI', n_neighbors=30)
sc.tl.umap(adata_normal, min_dist=0.3)
sc.tl.leiden(adata_normal, resolution=1.0, flavor='igraph', n_iterations=2)

adata_normal.write('processed/integrated/GSE268609_normal_scvi.h5ad')
m.save('processed/integrated/scvi_268609_normal/', overwrite=True)
```

VRAM 使用 ~12-18 GB、所要時間 30-90 分。

## 4.5 細胞型 annotation

**3 方法併用** (v4 で 1 reference 追加):
1. **Lazarov 2024 著者ラベル** (GSE268609 supplement にあれば最初に label transfer)
2. **Franjic 2022 (Sestan lab) reference** を `scvi.model.SCANVI` または `scarches` で投影
3. **GSE160189 (anterior/posterior axis reference)** を併用 → DG/CA/SUB/Sub-region の region 予測ラベルも同時付与

GSE160189 は本研究で region predictor としても使う。具体的には:

```python
# scripts/04_figure1/03b_region_predictor.py
import scanpy as sc, scvi
ref = sc.read('processed/per_dataset/GSE160189_qc.h5ad')
# region label (anterior/posterior, DG/CA1/CA3/SUB) が著者ラベルにあること前提
scvi.model.SCANVI.setup_anndata(ref, batch_key='sample_id',
                                 labels_key='region_label',
                                 unlabeled_category='Unknown')
sc_model = scvi.model.SCANVI(ref, n_layers=2, n_latent=20)
sc_model.train(max_epochs=300, accelerator='gpu')
sc_model.save('processed/integrated/scanvi_region_predictor/', overwrite=True)

# Discovery (GSE268609) に投影
query = sc.read('processed/integrated/GSE268609_normal_scvi.h5ad')
scvi.model.SCANVI.prepare_query_anndata(query, sc_model)
q_model = scvi.model.SCANVI.load_query_data(query, sc_model)
q_model.train(max_epochs=100, plan_kwargs={'weight_decay':0.0})
query.obs['predicted_region'] = q_model.predict()
query.obs['region_confidence'] = q_model.predict(soft=True).max(axis=1)
query.write('processed/integrated/GSE268609_normal_scvi_with_region.h5ad')
```

これにより Discovery 細胞それぞれに **predicted region (anterior/posterior, DG/CA/SUB) ラベルが付き**、後続 Figure 3 / 5 で region-specific 解析が可能になる。

### Marker-based 検証

```python
# scripts/04_figure1/03_annotate.py
markers = {
  'NSC_RGL':     ['GFAP','VIM','HOPX','NES','SOX2','PAX6','SLC1A3'],
  'aNSC':        ['GFAP','SOX2','MKI67','TOP2A','ASCL1'],
  'NeuroBlast':  ['DCX','NEUROD1','EOMES','SOX11'],
  'imGC':        ['DCX','PROX1','CALB2','STMN1','BHLHE22','SEMA3C'],
  'GC_mature':   ['PROX1','CALB1','SLC17A7','GRIA2'],
  'ExN_CA1':     ['SLC17A7','MEF2C','POU3F1','WFS1'],
  'ExN_CA3':     ['SLC17A7','GRIK4','CHGB'],
  'InN_PV':      ['GAD1','GAD2','PVALB'],
  'InN_SST':     ['GAD1','SST'],
  'InN_VIP':     ['GAD1','VIP'],
  'Astro':       ['AQP4','GFAP','SLC1A2','GJA1','ALDH1L1','S100B'],
  'Astro_react': ['CHI3L1','GFAP','VIM','SERPINA3','C3'],     # 加齢/炎症型
  'Micro_homeo': ['CSF1R','C1QA','P2RY12','TMEM119','CX3CR1'],
  'Micro_act':   ['SPP1','APOE','CD83','TNF','IL1B','CCL3'],   # 加齢炎症型
  'OPC':         ['PDGFRA','CSPG4','SOX10'],
  'Oligo':       ['MOG','MBP','PLP1','MOBP'],
  'Endo':        ['CLDN5','FLT1','PECAM1','VWF'],
  'VLMC':        ['DCN','COL1A2','PDGFRB'],
  'Peri':        ['PDGFRB','RGS5','CSPG4'],
}
sc.pl.dotplot(adata_normal, markers, groupby='leiden',
              standard_scale='var', save='_panel_D_celltype.pdf')
```

**重要**: aged microglia (Micro_act) と aged astrocyte (Astro_react) のラベルを最初から立てておく。これらは secretome 解析の主要 sender。

annotation 完成後、`adata.obs['celltype_l1']`、`adata.obs['celltype_l2']` を作って保存。

## 4.6 組成解析（Panel E-F）

```python
# scripts/04_figure1/04_composition.py
import seaborn as sns, matplotlib.pyplot as plt, scipy.stats as sst, pandas as pd
import numpy as np

# Panel E: 年齢勾配 scatter + LOWESS
comp = (adata_normal.obs.groupby(['donor_id','celltype_l2']).size()
        .unstack(fill_value=0).pipe(lambda d: d.div(d.sum(axis=1), axis=0)))
donor_age = adata_normal.obs.groupby('donor_id')['age_years'].first()
comp = comp.loc[donor_age.sort_values().index]

target_ct = ['Astro','Astro_react','Micro_homeo','Micro_act','OPC','Oligo',
             'Endo','VLMC','NSC_RGL','aNSC','NeuroBlast','imGC']

fig, axes = plt.subplots(3, 4, figsize=(16,10))
for ax, ct in zip(axes.flat, target_ct):
    if ct not in comp.columns:
        ax.set_visible(False); continue
    df = pd.DataFrame({'age': donor_age.values, 'frac': comp[ct].values})
    sns.regplot(data=df, x='age', y='frac', lowess=True, ax=ax,
                scatter_kws={'s':30,'alpha':0.7})
    r, p = sst.spearmanr(df['age'], df['frac'])
    ax.set_title(f'{ct}\nρ={r:.2f}, p={p:.1e}', fontsize=10)
plt.tight_layout()
plt.savefig('figures/figure1/panel_E_age_composition.pdf', bbox_inches='tight')

# Panel F: scCODA Bayesian compositional
import sccoda.util.cell_composition_data as ccd
import sccoda.util.comp_ana as ca

df = comp.reset_index().merge(donor_age.reset_index(), on='donor_id')
data = ccd.from_pandas(df, covariate_columns=['age_years'])
mod = ca.CompositionalAnalysis(data, formula='age_years', reference_cell_type='Oligo')
res = mod.sample_hmc(num_results=20000)
print(res.summary())
res.summary_extended().to_csv('results/figure1/sccoda_age_composition.csv')

# Panel D 補助: stacked bar (donor 順)
ax = comp.plot.bar(stacked=True, figsize=(12,4), width=1.0, colormap='tab20')
plt.savefig('figures/figure1/panel_D_composition_bar.pdf', bbox_inches='tight')
```

ドナー単位の比率を扱うことで cell-level の偽の自由度膨張を避ける。

## 4.7 Validation checkpoint (Figure 1)

- [ ] UMAP で全主要細胞型が分離（dataset でも年齢でも均等に混ざる、scIB metrics で `kBET`, `iLISI` を一応算出）
- [ ] **NSC, aNSC, neuroblast, imGC が個別クラスタとして分離**（後の figure の前提）
- [ ] 各 cluster の marker gene が dot plot で一意に立ち上がる
- [ ] Lazarov 2024 既報の所見 (NSC で aging signal 早期出現、aged microglia 拡大、aged astrocyte 反応性 shift) が再現
- [ ] aged astrocyte 減少、OPC 減少、endothelial 減少が確認できる（hippocampus aging の標準所見）
- [ ] 図 PDF が `figures/figure1/` に F1A-F 出揃う
- [ ] `results/figure1/cell_composition.csv` 保存
- [ ] `docs/DECISIONS.md` に QC 閾値、scVI ハイパラ、annotation 根拠を記録

→ `git tag figure1-v1` & HDD snapshot.

---

# Section 5: Figure 2 — Neurogenic lineage の精緻化

## 5.1 目的

NSC (quiescent / activated) → neuroblast → imGC → mature GC の軌跡を再構成し、各サブタイプの**年齢別存在量と転写状態**を示す。後続 secretome 解析の receiver/sender を確定する step。

## 5.2 想定 panel

- **A**: neurogenic lineage 細胞のみの subset UMAP (GSE268609 normal aging)
- **B**: pseudotime trajectory (Palantir / Slingshot)
- **C**: 主要 marker (GFAP, SOX2, NES, MKI67, DCX, PROX1, CALB1) の pseudotime expression curve
- **D**: NSC subtype (qNSC1/qNSC2/aNSC/pNSC、eLife 2024 のラベル) の年齢別 fraction、bootstrap CI
- **E**: Zhou シリーズ (GSE185277/553/198323) での再現性 — lineage signature scores の cohort 間 Spearman correlation
- **F**: imGC ML score (Zhou 2022 model) の年齢勾配 + GSE325391 利用可なら resilience との関係

## 5.3 解析パイプライン

```python
# scripts/04_figure2/01_lineage_subset.py
import scanpy as sc, scvi, palantir
import numpy as np, pandas as pd, scipy.stats as sst

adata = sc.read('processed/integrated/GSE268609_normal_scvi.h5ad')

# Astro は連続性参照として残す
linc = adata[adata.obs['celltype_l2'].isin(
    ['NSC_RGL','aNSC','NeuroBlast','imGC','GC_mature','Astro'])].copy()

# Subset 用に scVI 再学習 (少数細胞でも収束しやすい)
scvi.model.SCVI.setup_anndata(linc, batch_key='sample_id',
                              categorical_covariate_keys=['sex'])
m2 = scvi.model.SCVI(linc, n_layers=2, n_latent=15)
m2.train(max_epochs=300, accelerator='gpu')
linc.obsm['X_scVI_lin'] = m2.get_latent_representation()

sc.pp.neighbors(linc, use_rep='X_scVI_lin', n_neighbors=20)
sc.tl.umap(linc, min_dist=0.2)
sc.tl.leiden(linc, resolution=0.6, key_added='leiden_lin')

# Pseudotime (Palantir)
linc_p = linc.copy()
sc.pp.normalize_total(linc_p, target_sum=1e4); sc.pp.log1p(linc_p)
linc_p.X = linc_p.X.toarray() if hasattr(linc_p.X, 'toarray') else linc_p.X

dm_res = palantir.utils.run_diffusion_maps(linc_p, n_components=10)
ms_data = palantir.utils.determine_multiscale_space(linc_p)

# Start cell: GFAP+/SOX2+/MKI67- の qNSC を選ぶ
gfap_idx = linc_p.var_names.get_loc('GFAP')
mki67_idx = linc_p.var_names.get_loc('MKI67') if 'MKI67' in linc_p.var_names else None
start_mask = (linc_p.obs['celltype_l2']=='NSC_RGL')
start_cell = linc_p.obs.index[start_mask][0]   # 仮、後で markers で精選

pr_res = palantir.core.run_palantir(ms_data, start_cell, num_waypoints=500)
linc.obs['pseudotime'] = pr_res.pseudotime
linc.obs['differentiation_potential'] = pr_res.entropy
linc.write('processed/integrated/GSE268609_lineage.h5ad')
```

## 5.4 重要な解析判断

**警告**: 成体ヒト snRNA-seq で NSC・neuroblast は数十〜数百細胞しか検出されないことが多い。Trajectory が unstable になりやすい。

**対策**:
- **Bootstrap stability**: `scanpy.tl.leiden` を 100 回 random seed で繰り返し ARI を出して安定性を panel S として残す
- **Cohort 横断検証**: GSE268609 と Zhou (GSE198323) を**それぞれ独立に解析**して同じ軌跡が出るかを確認
- **反証データでの実行**: Franjic GSE186538 で**同じパイプラインを走らせる** → 「ヒトでは neuroblast/imGC が他種比で少ない」という Franjic 報告が再現するか。これはネガコン的 figure として価値あり (`figures/figure2/franjic_lineage_comparison.pdf`)

## 5.5 Validation checkpoint (Figure 2)

- [ ] NSC marker (GFAP+SOX2+VIM+) と astrocyte marker (AQP4+SLC1A2+) が UMAP で連続だが分離可能
- [ ] pseudotime が GFAP→DCX→PROX1 の順を再現
- [ ] GSE268609 と Zhou で同様の軌跡（少なくとも順序）が出る
- [ ] Ma/Pei eLife 2024 の qNSC1/qNSC2/pNSC/aNSC ラベルとおおむね対応
- [ ] Bootstrap stability score が報告される

→ `git tag figure2-v1` & snapshot.

---

# Section 6: Figure 3 — Secretome DE と meta-analysis

**ここがプロジェクトの中核**。

## 6.1 目的

HPA / UniProt / MatrisomeDB で定義した secretome 遺伝子について、**主要ニッチ細胞型ごとに年齢依存的差分発現**を、ドナー pseudobulk で統計検定し、データセット間で meta-analyze する。

## 6.2 想定 panel

- **A**: secretome 遺伝子定義の UpSet plot (HPA × UniProt × MatrisomeDB × SignalP)
- **B**: 細胞型 × age での secretome DE 数の summary heatmap
- **C**: Discovery (GSE268609) の volcano plot — 細胞型 4 panel (Astro / Micro / Endo / NSC)
- **D**: Top 20 robust hits の forest plot (5+ cohort × cell type)
- **E**: Functional category 別 (cytokine / GF / matrisome / neuropeptide) の年齢相関分布
- **F**: 既知の老化マーカー (CHI3L1, CCL2, GDF15, CLU, APOE, SPP1 など) を highlight した age trajectory

## 6.3 Secretome 遺伝子セット構築

```python
# scripts/05_figure3/01_build_secretome_set.py
import pandas as pd
from collections import Counter
import upsetplot

hpa  = pd.read_csv('refs/hpa_secretome.tsv', sep='\t')['Gene'].dropna().unique()
uni  = pd.read_csv('refs/uniprot_secreted_human.tsv', sep='\t')
uni_genes = uni['Gene Names'].dropna().str.split().explode().unique()
matri = pd.read_excel('refs/Hs_Matrisome_Masterlist.xlsx',
                     sheet_name='Hs_Matrisome_Masterlist')
matri_genes = matri['Gene Symbol'].dropna().unique()

# 階層化セット
sec = {
  'HPA': set(hpa),
  'UniProt': set(uni_genes),
  'Matrisome': set(matri_genes),
}

# Core: 2 つ以上で支持されるもの (厳格)
c = Counter()
for s in sec.values():
    for g in s: c[g] += 1
core = {g for g,n in c.items() if n >= 2}

pd.DataFrame({'gene': sorted(core), 'category': ['core_secretome']*len(core)})\
  .to_csv('refs/secretome_core.csv', index=False)

# UpSet plot (Panel A)
import upsetplot, matplotlib.pyplot as plt
membership = upsetplot.from_contents(sec)
upsetplot.plot(membership, show_counts=True)
plt.savefig('figures/figure3/panel_A_upset.pdf', bbox_inches='tight')
```

**カテゴリ分類**を必ず付ける:
- Cytokine / chemokine (IL, CXCL, CCL, TNF…)
- Growth factor / neurotrophin (NGF, BDNF, GDNF, IGF, FGF, VEGF, PDGF, TGFB, BMP…)
- Neuropeptide (NPY, SST, CRH, VIP…)
- Matrisome → Core matrisome / Matrisome-associated (MatrisomeDB の subdivision を継承)
- Synaptic adhesion 分泌型 (SPARCL1, NRXN/NLGN secreted forms など)

## 6.4 Pseudobulk DE（主軸）

**重要原則**: 細胞単位 DE は pseudoreplication で偽陽性。**必ず donor × celltype の集計に戻して DESeq2 / edgeR**。

```python
# scripts/05_figure3/02_pseudobulk_de_per_cohort.py
import scanpy as sc, pandas as pd, numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

def run_pseudobulk_de(adata, dataset_name):
    sec = set(pd.read_csv('refs/secretome_core.csv')['gene'])
    results = []
    for ct in ['Astro','Astro_react','Micro_homeo','Micro_act','OPC','Oligo',
               'Endo','VLMC','NSC_RGL','aNSC','NeuroBlast','imGC','GC_mature']:
        if ct not in adata.obs['celltype_l2'].unique(): continue
        sub = adata[adata.obs['celltype_l2']==ct]

        # 最低 cell 数を満たす donor のみ
        n_cells = sub.obs.groupby('donor_id').size()
        keep_donors = n_cells[n_cells >= 10].index
        if len(keep_donors) < 5: continue   # 統計的に不可能
        sub = sub[sub.obs['donor_id'].isin(keep_donors)]

        # Pseudobulk (sum over cells per donor)
        pb = sc.get.aggregate(sub, by='donor_id', func='sum')
        counts = pb.X.toarray() if hasattr(pb.X, 'toarray') else pb.X
        counts_df = pd.DataFrame(counts, index=pb.obs_names, columns=pb.var_names).T

        meta = adata.obs.groupby('donor_id').agg({
            'age_years':'first','sex':'first','pmi_hours':'first',
            'dataset':'first'
        }).loc[counts_df.columns]
        meta['age_years'] = meta['age_years'].astype(float)

        # DESeq2 (pyDESeq2)
        dds = DeseqDataSet(counts=counts_df.T, metadata=meta,
                           design_factors=['sex','age_years'])
        dds.deseq2()
        stat = DeseqStats(dds, contrast=['age_years', None, None])
        stat.summary()
        res = stat.results_df.assign(celltype=ct, dataset=dataset_name)
        res['in_secretome'] = res.index.isin(sec)
        results.append(res.reset_index().rename(columns={'index':'gene'}))

    return pd.concat(results) if results else pd.DataFrame()

# 各 cohort で実行
for cohort, path in [
    ('GSE268609', 'processed/integrated/GSE268609_normal_scvi.h5ad'),
    ('GSE185277', 'processed/integrated/GSE185277_scvi.h5ad'),
    ('GSE185553', 'processed/integrated/GSE185553_scvi.h5ad'),
    ('GSE198323', 'processed/integrated/GSE198323_normal_scvi.h5ad'),
    ('GSE163737', 'processed/integrated/GSE163737_human_scvi.h5ad'),
]:
    a = sc.read(path)
    de = run_pseudobulk_de(a, cohort)
    de.to_csv(f'results/figure3/de_{cohort}_per_celltype.csv', index=False)
```

R `muscat` を併用すると cell-level 情報を保ったまま pseudobulk DE 可:

```r
# scripts/05_figure3/02_muscat_de.R
library(muscat); library(SingleCellExperiment); library(zellkonverter)
sce <- readH5AD("processed/integrated/GSE268609_normal_scvi.h5ad")
sce$age_group <- factor(sce$age_group, levels=c("20-40","40-60","60-80","80-100"))
sce <- prepSCE(sce, kid="celltype_l2", gid="age_group", sid="donor_id", drop=TRUE)
pb  <- aggregateData(sce, assay="counts", fun="sum",
                     by=c("cluster_id","sample_id"))
res <- pbDS(pb, design = ~ sex + age_group, coef = c("age_group80-100"))
saveRDS(res, "results/figure3/muscat_age80vs20.rds")
```

**3 つのモデル**を併走させ、結果一致のものを robust hit とする:
1. **連続モデル**: `~ sex + age_years`
2. **4 群比較**: `~ sex + age_group`、young (20-40) を reference
3. **マウス相同遺伝子の老化変化との一致** (Tabula Muris Senis 等)

## 6.5 Meta-analysis

```r
# scripts/05_figure3/03_meta_analysis.R
library(metafor); library(tidyverse)

cohorts <- c('GSE268609','GSE185277','GSE185553','GSE198323','GSE163737')
de_list <- map(cohorts, ~ read_csv(paste0('results/figure3/de_', .x, '_per_celltype.csv')))
names(de_list) <- cohorts

long <- bind_rows(de_list, .id='cohort')
sec <- read_csv('refs/secretome_core.csv') %>% pull(gene)
long_sec <- long %>% filter(gene %in% sec)

meta_res <- long_sec %>%
  filter(!is.na(log2FoldChange), !is.na(lfcSE)) %>%
  group_by(gene, celltype) %>%
  filter(n() >= 3) %>%   # 最低 3 cohort
  group_modify(~ {
    rma_fit <- tryCatch(
      metafor::rma(yi=.x$log2FoldChange, sei=.x$lfcSE, method='REML'),
      error=function(e) NULL)
    if (is.null(rma_fit)) return(tibble())
    tibble(meta_b = rma_fit$b[1,1],
           meta_se = rma_fit$se,
           meta_p = rma_fit$pval,
           meta_I2 = rma_fit$I2,
           n_cohort = nrow(.x))
  }) %>% ungroup() %>%
  mutate(meta_padj = p.adjust(meta_p, method='BH'))

write_csv(meta_res, 'results/figure3/secretome_meta_analysis.csv')
```

**Robust hit の定義**: `meta_padj < 0.05` AND `n_cohort >= 4` AND `I2 < 50%` (heterogeneity 低い) AND `符号一致 cohort 数 >= 4/5`。

## 6.6 Validation checkpoint (Figure 3)

- [ ] secretome 遺伝子セットの定義が再現可能 (`refs/secretome_core.csv` にハッシュつき)
- [ ] 各細胞型の DE 結果に明らかな biological sense (例: aged microglia で C1Q complement、IL-1β、CCL3/4 上昇)
- [ ] 既知の老化マーカー (CHI3L1/YKL-40 in astrocyte、GDF15 in senescence、APOE/SPP1 in microglia) が捉えられる
- [ ] GSE268609 と Zhou で同方向の hit が >50%
- [ ] Robust hit が細胞型あたり数十〜100 程度
- [ ] Franjic データでも astrocyte/microglia の方向は一致する (neurogenic lineage は乖離してOK)

→ `git tag figure3-v1` & snapshot.

---

# Section 7: Figure 4 — Cell-cell communication + Regulatory

## 7.1 目的

Figure 3 で同定した secretome 因子について、**どの細胞が send し、どの細胞が receive するか**、**通信が加齢でどう reshape するか**を CellChat v2 と NicheNet v2 で示す。さらに **multiome (snATAC) を活用して cis-regulatory landscape の年齢変化**を解析し、Lazarov 2024 にない方法論的 novelty を確立する。

## 7.2 想定 panel

- **A**: 全細胞型間の通信強度 chord (young vs old)
- **B**: 主要 secreted pathway (NOTCH, BMP, WNT, IL1, TGFβ, FGF, VEGF, IGF など) のシグナル強度の年齢別 bubble plot
- **C**: NSC/imGC を receiver にした NicheNet 上位リガンドランキングと推定下流ターゲット
- **D**: 加齢で**新規に出現/消失**する LR ペア (differential signaling)
- **E**: snATAC での secretome 遺伝子 peak の age correlation heatmap (multiome 独自軸)
- **F**: chromVAR で aging-associated TF motifs (Astro/Micro/Endo)

## 7.3 CellChat v2 (secretome restricted)

```r
# scripts/06_figure4/01_cellchat_secreted_only.R
library(CellChat); library(zellkonverter); library(SeuratDisk)

ad <- readH5AD("processed/integrated/GSE268609_normal_scvi.h5ad")
seu <- as.Seurat(ad)
seu$age_bin <- cut(seu$age_years, breaks=c(0,40,60,80,120),
                   labels=c("young","middle","old","oldest"))

# Secretome lens: secreted signaling category のみに subset
CellChatDB.use <- subsetDB(CellChatDB.human, search='Secreted Signaling')

run_cc <- function(seu_sub, group_lab){
  data.in <- GetAssayData(seu_sub, slot="data")
  meta <- seu_sub@meta.data
  cc <- createCellChat(object=data.in, meta=meta, group.by="celltype_l2")
  cc@DB <- CellChatDB.use
  cc <- subsetData(cc) |> identifyOverExpressedGenes() |>
        identifyOverExpressedInteractions() |> projectData(PPI.human)
  cc <- computeCommunProb(cc, type="triMean")
  cc <- filterCommunication(cc, min.cells=10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  saveRDS(cc, paste0("results/figure4/cellchat_", group_lab, ".rds"))
  cc
}

cc_y <- run_cc(subset(seu, age_bin=="young"),  "young")
cc_o <- run_cc(subset(seu, age_bin=="oldest"), "oldest")

cc_list <- list(young=cc_y, oldest=cc_o)
cc_merge <- mergeCellChat(cc_list, add.names=names(cc_list))
gg1 <- compareInteractions(cc_merge, show.legend=FALSE, group=c(1,2))
gg2 <- netVisual_diffInteraction(cc_merge, weight.scale=TRUE)
rankNet(cc_merge, mode="comparison", stacked=TRUE)
ggsave("figures/figure4/panel_B_pathway_rank.pdf", width=6, height=8)
```

## 7.4 NicheNet v2 (receiver 視点)

NSC, imGC, neuroblast を receiver、加齢で発現変化する遺伝子を target に、上流リガンドを推定:

```r
# scripts/06_figure4/02_nichenet_receiver.R
library(nichenetr); library(Seurat); library(tidyverse)

lr_network    <- readRDS("refs/nichenet/lr_network_human_21122021.rds")
ligand_target <- readRDS("refs/nichenet/ligand_target_matrix_nsga2r_final.rds")
weighted      <- readRDS("refs/nichenet/weighted_networks_nsga2r_final.rds")

seu <- readRDS("processed/integrated/GSE268609_seurat.rds")

run_nichenet <- function(seu, receiver, sender_cells, condition_col='age_bin',
                         test='oldest', ref='young'){
  seu_recv <- subset(seu, celltype_l2 == receiver)
  Idents(seu_recv) <- condition_col
  de_recv <- FindMarkers(seu_recv, ident.1=test, ident.2=ref,
                         min.pct=0.1, logfc.threshold=0.25) %>%
             rownames_to_column("gene")
  geneset_oi <- de_recv %>% filter(p_val_adj<0.05) %>% pull(gene)

  seu_test <- subset(seu, !!sym(condition_col) == test)
  expr_lig <- get_expressed_genes(sender_cells, seu_test, pct=0.05)
  ligands_oi <- intersect(expr_lig, lr_network$from)

  bg <- rownames(seu_recv)[Matrix::rowSums(seu_recv@assays$RNA@counts > 0) > 50]
  ligand_act <- predict_ligand_activities(
    geneset = geneset_oi,
    background_expressed_genes = bg,
    ligand_target_matrix = ligand_target,
    potential_ligands = ligands_oi)
  ligand_act %>% arrange(-aupr_corrected)
}

senders <- c("Astro","Astro_react","Micro_homeo","Micro_act","Endo","VLMC","OPC","ExN_CA1","ExN_CA3")
for (recv in c("NSC_RGL","NeuroBlast","imGC")){
  res <- run_nichenet(seu, receiver=recv, sender_cells=senders)
  write_csv(res, paste0("results/figure4/nichenet_", recv, "_old_vs_young.csv"))
}
```

## 7.5 Multiome regulatory analysis (v2 独自軸)

```python
# scripts/06_figure4/03_atac_secretome_peaks.py
import muon as mu, scanpy as sc, pandas as pd, numpy as np
import scipy.stats as sst, scipy.sparse as sp

mdata = mu.read('processed/integrated/GSE268609_normal_multi.h5mu')
atac = mdata.mod['atac']
sec = pd.read_csv('refs/secretome_core.csv')['gene'].tolist()

# secretome 遺伝子の近傍 peak を annotation
# muon の atac.tl.add_peak_annotation で gene linkage を構築済み前提
sec_peaks = atac.var.query('gene_name in @sec').index

results = []
for ct in ['Astro','Astro_react','Micro_homeo','Micro_act','Endo','NSC_RGL']:
    mask = atac.obs['celltype_l2'] == ct
    sub_atac = atac[mask, sec_peaks]
    age = sub_atac.obs['age_years'].values
    donors = sub_atac.obs['donor_id'].values
    X = sub_atac.X.toarray() if sp.issparse(sub_atac.X) else sub_atac.X
    pb = pd.DataFrame(X, index=donors, columns=sec_peaks)
    pb_sum = pb.groupby(level=0).sum()
    pb_age = pd.Series(age).groupby(donors).first()
    for peak in pb_sum.columns:
        r, p = sst.spearmanr(pb_age, pb_sum[peak])
        results.append({'celltype': ct, 'peak': peak,
                        'gene': atac.var.loc[peak, 'gene_name'],
                        'rho': r, 'p': p})
pd.DataFrame(results).to_csv('results/figure4/secretome_peak_age_corr.csv', index=False)
```

```r
# scripts/06_figure4/04_chromvar_motifs.R
library(Signac); library(JASPAR2020); library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)

seu <- readRDS('processed/integrated/GSE268609_seurat.rds')
DefaultAssay(seu) <- 'peaks'
pwms <- getMatrixSet(JASPAR2020, opts=list(species=9606, all_versions=FALSE))
seu  <- AddMotifs(seu, genome=BSgenome.Hsapiens.UCSC.hg38, pfm=pwms)
seu  <- RunChromVAR(seu, genome=BSgenome.Hsapiens.UCSC.hg38)
DefaultAssay(seu) <- 'chromvar'

# 年齢相関 (donor x motif の pseudobulk)
saveRDS(seu, 'processed/integrated/GSE268609_seurat_chromvar.rds')

# Aging-associated TF motif の同定:
# - donor 単位で motif activity を集計
# - 各 cell type で age との Spearman correlation
# - NFkB, AP-1, STAT3 (炎症), SOX/REST (neural identity 喪失) 等が出る想定
```

## 7.6 Validation checkpoint (Figure 4)

- [ ] CellChat の通信強度差が biological sense
- [ ] 加齢で増加する pathway: TGFβ, IL1, complement, SPP1 を期待
- [ ] 加齢で減少: BMP, NOTCH, IGF, BDNF, WNT
- [ ] **NicheNet トップリガンドが Figure 3 DE hit と >50% 重なる（Figure 3 ↔ 4 のクロスバリデーション）**
- [ ] snATAC peak の age correlation 方向と RNA 発現方向がおおむね一致
- [ ] aging-associated TF motif: NFkB, AP-1, STAT3, SOX/REST 等が出る
- [ ] AD/disease データを混ぜていないこと（純粋に「加齢」を見ている）の確認

→ `git tag figure4-v1` & snapshot.

---

# Section 8: Figure 5 — 二重空間検証 + Regional Specificity

## 8.1 目的

Figure 3-4 で見つけた secretome 因子が、**Visium ペア (GSE264624 snRNA + GSE264692 spatial) と spatial+FISH (GSE248545) の両プラットフォーム**で DG 顆粒層〜SGZ 近傍に局在しているかを確認。さらに **v4 新設パネル**として、aging secretome シグナルが SGZ/DG 特異か pan-hippocampal かを GSE160189 + GSE264692 の連携で判定する。

## 8.2 想定 panel (v4 改訂)

- **A**: GSE264692 Visium での DG/CA1/CA3/SUB の spatial domain (BayesSpace or NMF)
- **B**: GSE264624 snRNA を reference に GSE264692 Visium へ cell2location デコン → NSC/imGC の SGZ 局在
- **C**: 主要 secretome hit (3-5 genes) の Visium spatial expression map
- **D**: GSE248545 (FISH-based) で同じ遺伝子が DG に局在することを直交検証
- **E**: SpatialDM などで spatial-aware LR pair の同定
- **F**: **(v4 新設) Regional specificity panel** — secretome aging score の DG vs CA vs SUB 比較。GSE160189 で学習した region predictor を Discovery に適用 → region 別 secretome aging signature を計算 → Visium 空間で再投影して整合性確認

## 8.3 GSE264624 ↔ GSE264692 ペアでの cell2location

**v4 重要ポイント**: cell2location は reference snRNA-seq と target spatial が**同じ組織コンテキスト**である方が精度が高い。Tran 2025 の snRNA (GSE264624) を reference に、同じドナーの Visium (GSE264692) にデコンする self-consistency テストが可能 → これが通れば、別 reference (GSE268609) を使ったデコンの妥当性も担保される。

```python
# scripts/07_figure5/01_cell2location_paired.py
import cell2location as c2l, scanpy as sc, anndata as ad, numpy as np
import pandas as pd, glob, os

# Step 1: GSE264624 snRNA を reference signature に
ref_tran = sc.read('processed/integrated/GSE264624_tran_snrna_qc.h5ad')
ref_tran.X = ref_tran.layers['counts']

c2l.models.RegressionModel.setup_anndata(ref_tran, batch_key='sample_id',
                                          labels_key='celltype_l2')
ref_model_tran = c2l.models.RegressionModel(ref_tran)
ref_model_tran.train(max_epochs=300, accelerator='gpu', batch_size=2500)
ref_model_tran.export_posterior(ref_tran, sample_kwargs={'num_samples':1000})
inf_aver_tran = ref_tran.varm['means_per_cluster_mu_fg']

# Step 2: 並行して GSE268609 (Discovery) snRNA も reference signature に
ref_disc = sc.read('processed/integrated/GSE268609_normal_scvi.h5ad')
ref_disc.X = ref_disc.layers['counts']
c2l.models.RegressionModel.setup_anndata(ref_disc, batch_key='sample_id',
                                          labels_key='celltype_l2')
ref_model_disc = c2l.models.RegressionModel(ref_disc)
ref_model_disc.train(max_epochs=300, accelerator='gpu', batch_size=2500)
ref_model_disc.export_posterior(ref_disc, sample_kwargs={'num_samples':1000})
inf_aver_disc = ref_disc.varm['means_per_cluster_mu_fg']

# Step 3: GSE264692 Visium に対して両 reference でデコン → 一致度を確認
metadata = pd.read_csv('refs/metadata_master.csv')
visium_samples = metadata.query('dataset == "GSE264692"')['sample_id'].tolist()

for sample in visium_samples:
    visium_dir = f'raw/GSE264692_tran2025_visium/{sample}/'
    st = sc.read_visium(visium_dir)
    st.var_names_make_unique()

    for ref_name, inf_aver in [('tran', inf_aver_tran), ('disc', inf_aver_disc)]:
        st_copy = st.copy()
        c2l.models.Cell2location.setup_anndata(st_copy, batch_key='sample_id')
        st_model = c2l.models.Cell2location(st_copy, cell_state_df=inf_aver,
                                             N_cells_per_location=5,
                                             detection_alpha=20)
        st_model.train(max_epochs=30000, accelerator='gpu', batch_size=None)
        st_model.export_posterior(st_copy, sample_kwargs={'num_samples':1000})
        st_copy.write(f'processed/integrated/{sample}_c2l_{ref_name}.h5ad')

# Step 4: 両 reference のデコン結果を比較 (sanity check)
# 主要細胞型 (NSC, Astro, Micro) のデコン値が空間的に一致すれば self-consistent
```

両 reference でデコン結果が一致 → Discovery (GSE268609) の reference を使った解析の妥当性が空間レベルで担保される。これは**reviewer に対する強力な防御**になる。

## 8.4 Secretome hit の spatial mapping

```python
# scripts/07_figure5/02_secretome_spatial_map.py
import scanpy as sc, matplotlib.pyplot as plt, glob

top_hits = ['CHI3L1','CCL2','SERPINA3','SPP1','APOE','C1QA',
            'VEGFA','BMP4','NOTCH2','BDNF','IGF1','TGFB1']

for sample_path in glob.glob('processed/integrated/*_c2l_disc.h5ad'):
    st = sc.read(sample_path)
    fig, axes = plt.subplots(3, 4, figsize=(20,15))
    for ax, g in zip(axes.flat, top_hits):
        if g not in st.var_names: ax.set_visible(False); continue
        sc.pl.spatial(st, color=g, ax=ax, show=False, title=g)
    sample_name = sample_path.split('/')[-1].replace('_c2l_disc.h5ad','')
    plt.savefig(f'figures/figure5/panel_C_{sample_name}.pdf', bbox_inches='tight')
```

## 8.5 GSE248545 spatial+FISH 検証

GSE248545 のデータ形式は MERFISH/Xenium または Stereo-seq の可能性。`spatialdata-io` で統一形式に:

```python
# scripts/07_figure5/03_FISH_validation.py
import spatialdata_io as sdio
# 形式に応じて以下のいずれか:
# sdata = sdio.merscope('raw/GSE248545_spatial_FISH/<sample>/')
# sdata = sdio.xenium('raw/GSE248545_spatial_FISH/<sample>/')

# DG ROI を手動 / segmentation で定義
# 主要 hit が GSE264692 と GSE248545 両方の DG ROI で発現上昇しているか
# 二重 violin plot を panel D として作成
```

## 8.6 Regional specificity panel (v4 新設、Panel F)

「SGZ/DG 特異か pan-hippocampal か」を空間レベルで定量する。

```python
# scripts/07_figure5/04_regional_specificity.py
import scanpy as sc, pandas as pd, numpy as np, scipy.stats as sst
import matplotlib.pyplot as plt, seaborn as sns

# Step 1: meta-analysis robust hit からセルタイプ別 secretome aging signature を作る
meta = pd.read_csv('results/figure3/secretome_meta_analysis.csv')
robust = meta.query('meta_padj < 0.05 and n_cohort >= 4 and abs(meta_b) > 0.5')

signatures = {}
for ct in robust['celltype'].unique():
    ct_genes = robust.query(f'celltype == "{ct}"')
    signatures[f'{ct}_aging_up']   = ct_genes.query('meta_b > 0')['gene'].tolist()
    signatures[f'{ct}_aging_down'] = ct_genes.query('meta_b < 0')['gene'].tolist()

# Step 2: Discovery (GSE268609) に scanpy で signature score を付与
adata = sc.read('processed/integrated/GSE268609_normal_scvi_with_region.h5ad')
for sig_name, genes in signatures.items():
    sc.tl.score_genes(adata, gene_list=genes, score_name=sig_name)

# Step 3: predicted_region (DG/CA1/CA3/SUB) で signature score を比較
region_summary = []
for sig_name in signatures:
    for region in adata.obs['predicted_region'].unique():
        sub = adata[adata.obs['predicted_region'] == region]
        # ドナー単位の平均で donor n の独立性を保つ
        donor_means = sub.obs.groupby('donor_id')[sig_name].mean()
        donor_age = sub.obs.groupby('donor_id')['age_years'].first()
        r, p = sst.spearmanr(donor_age, donor_means)
        region_summary.append({'signature': sig_name, 'region': region,
                               'rho_age': r, 'p_age': p,
                               'n_donors': len(donor_means)})
region_df = pd.DataFrame(region_summary)
region_df.to_csv('results/figure5/regional_specificity.csv', index=False)

# Heatmap: signature × region の age correlation
pivot = region_df.pivot(index='signature', columns='region', values='rho_age')
fig, ax = plt.subplots(figsize=(6, max(4, len(pivot)*0.3)))
sns.heatmap(pivot, cmap='RdBu_r', center=0, annot=True, fmt='.2f', ax=ax,
            cbar_kws={'label': 'Spearman ρ (age vs signature score)'})
ax.set_title('Regional specificity of aging secretome signatures')
plt.savefig('figures/figure5/panel_F_regional_specificity.pdf', bbox_inches='tight')

# Step 4: GSE264692 Visium 空間で再投影 → 各 spot に signature score を計算
# DG ROI の spot で signature score 上昇 vs CA/SUB ROI で計算
for sample_path in glob.glob('processed/integrated/*_c2l_disc.h5ad'):
    st = sc.read(sample_path)
    for sig_name, genes in signatures.items():
        avail = [g for g in genes if g in st.var_names]
        if len(avail) >= 5:
            sc.tl.score_genes(st, gene_list=avail, score_name=sig_name)
    # 各 spot の region label は BayesSpace または NMF で別途付与済み前提
    # signature score を region で box plot
```

**期待される所見**:
- もし aging secretome signature が **DG で最強、SUB で弱い** → SGZ/DG niche 特異性が確認される (主張の正当化)
- もし **全領域で均一** → SGZ 特異性は弱く、論文ストーリーを「pan-hippocampal aging niche」に修正する必要
- 細胞型ごとに region 特異性が異なる可能性: 例えば microglia aging signature は pan-HPC、NSC niche aging は DG 特異、など → これ自体が興味深い知見

## 8.7 Spatial-aware LR (オプション)

```r
# scripts/07_figure5/05_spatial_cellchat.R
library(CellChat)
# CellChat v2 の spatial mode
# spatial coordinates を入力して、physically proximal cells のみで LR を計算
```

## 8.8 Validation checkpoint (Figure 5)

- [ ] GSE264692 Visium で DG 領域が PROX1+ で正しく検出される
- [ ] **GSE264624 snRNA reference と GSE268609 Discovery reference の両方でデコン結果が一致** (self-consistency check)
- [ ] NSC/imGC のデコンシグナルが SGZ 相当の granule layer 内側に集中
- [ ] Figure 3 の hit secretome 遺伝子が SGZ-relevant spot で発現 (少なくとも一部)
- [ ] **Tran 2025 と GSE248545 の主要結論が一致** (片方でしか出ない hit はリスト化して discussion で言及)
- [ ] FISH/MERFISH レベルで個別細胞での発現が確認される
- [ ] **Regional specificity panel: aging secretome signature の DG vs CA vs SUB 差が定量化される**

→ `git tag figure5-v1` & snapshot.

---

# Section 9: Figure 6 — クロス種比較

## 9.1 目的

ヒトで見つけた加齢 secretome 変化のうち、マカク・マウスで保存されているもの (種を超えた基本機構) と、ヒト固有のものを切り分ける。GSE163737 のヒト+マカク同一プラットフォームを活用。

## 9.2 想定 panel

- **A**: 種間統合 UMAP (GSE163737 + GSE268609 ヒト + Hochgerner マウス)
- **B**: 加齢 secretome の種間方向一致率 heatmap
- **C**: human-conserved aging secretome (3 種で同方向)
- **D**: human-specific aging secretome (ヒトのみ)
- **E**: positive control: METTL7B (Franjic 報告の primate-specific) が捕捉される
- **F**: ortholog-restricted GO enrichment

## 9.3 Ortholog 統合

```python
# scripts/08_figure6/01_ortholog_integration.py
import scanpy as sc, anndata as ad, pandas as pd, scvi

# 1:1 ortholog table (BioMart で取得済)
orth = pd.read_csv('refs/orthologs_hsa_mml_mmu.csv')
# columns: human_symbol, macaque_symbol, mouse_symbol

# 各種データを load
hsa = sc.read('processed/integrated/GSE268609_normal_scvi.h5ad')
mml = sc.read('processed/integrated/GSE163737_macaque.h5ad')
mmu = sc.read('processed/integrated/Hochgerner2018_mouse.h5ad')

# Macaque/Mouse の遺伝子名を ortholog 経由で human symbol に rename
def rename_to_human(adata, species_col, mapping):
    m = mapping.dropna(subset=['human_symbol', species_col])
    m = m.drop_duplicates(subset=species_col).set_index(species_col)['human_symbol']
    adata.var['human_symbol'] = adata.var_names.map(m)
    keep = adata.var['human_symbol'].notna() & ~adata.var['human_symbol'].duplicated()
    adata = adata[:, keep].copy()
    adata.var_names = adata.var['human_symbol'].values
    return adata

mml = rename_to_human(mml, 'macaque_symbol', orth)
mmu = rename_to_human(mmu, 'mouse_symbol', orth)

# 共通遺伝子で concat
common_genes = sorted(set(hsa.var_names) & set(mml.var_names) & set(mmu.var_names))
hsa = hsa[:, common_genes]; mml = mml[:, common_genes]; mmu = mmu[:, common_genes]

combined = ad.concat([hsa, mml, mmu], label='species', keys=['human','macaque','mouse'])
combined.var_names_make_unique()

# scVI で species を batch として
scvi.model.SCVI.setup_anndata(combined, batch_key='species',
                              categorical_covariate_keys=['sample_id'])
mx = scvi.model.SCVI(combined, n_layers=2, n_latent=30)
mx.train(max_epochs=400, accelerator='gpu')
combined.obsm['X_scVI_xspecies'] = mx.get_latent_representation()
combined.write('processed/integrated/cross_species.h5ad')
```

## 9.4 種間 DE 比較

```r
# scripts/08_figure6/02_compare_aging_directions.R
library(tidyverse)

hsa_de <- read_csv('results/figure3/secretome_meta_analysis.csv')
mml_de <- read_csv('results/figure6/macaque_aging_de.csv')   # GSE163737 macaque で別途実行
mmu_de <- read_csv('results/figure6/mouse_aging_de.csv')     # Hochgerner + Tabula Muris Senis

merged <- hsa_de %>%
  inner_join(mml_de, by='gene', suffix=c('_hsa','_mml')) %>%
  inner_join(mmu_de %>% rename(log2FC_mmu=log2FC), by='gene')

merged <- merged %>% mutate(
  conserved = sign(meta_b)==sign(log2FC_mml) & sign(meta_b)==sign(log2FC_mmu),
  human_specific = sign(meta_b) != sign(log2FC_mmu) & abs(log2FC_mmu) < 0.25
)
write_csv(merged, 'results/figure6/cross_species_secretome.csv')
```

## 9.5 注意点

- snRNA-seq vs scRNA-seq の差を補正するため、**protein-coding gene のみ、HVG を ortholog 共通遺伝子で取る**
- Hochgerner 2018 単独では年齢比較が限定的なので **Tabula Muris Senis** や **Harris 2021 (DG aging mouse)** を併用
- GSE163737 のマカクは年齢勾配が連続的 → マカクで強い "aged secretome signature" を抽出 → ヒトで test という方向の解析が可能

## 9.6 Validation checkpoint (Figure 6)

- [ ] 主要細胞型 (GFAP, AQP4, AIF1) が種を超えて整列
- [ ] METTL7B などの primate-specific positive control が再現
- [ ] aged microglia 炎症 (IL1B, TNF, CCL3) は種を超えて保存
- [ ] ヒト特異 hit が biological sense (POSTN, NRG3, primate-specific paralogs)
- [ ] Polzer et al. 2025 cross-analysis の知見と整合

→ `git tag figure6-v1` & snapshot.

---

# Section 10: Figure 7 — 統合解釈

## 10.1 目的

Figure 3-6 の hit を **pathway / TF activity / 既存プロテオミクスデータ** に投影し、機能的意義を提示。

## 10.2 想定 panel

- **A**: Cell type × age での pathway activity (decoupler の PROGENy / Dorothea)
- **B**: 上位 hit に対する CSF aging proteomics (Mayo, Seattle Aging cohort 等) との照合
- **C**: hit が含まれる SenMayo / SASP atlas との overlap
- **D**: ヒト海馬 GTEx bulk での年齢相関による orthogonal 検証
- **E**: GSE325391 利用可なら cognitive resilience 関連 secretome の特定
- **F**: Summary 模式図 (BioRender)

## 10.3 decoupler

```python
# scripts/09_figure7/01_pathway_activity.py
import decoupler as dc, scanpy as sc

adata = sc.read('processed/integrated/GSE268609_normal_scvi.h5ad')
progeny = dc.get_progeny(organism='human', top=500)
collectri = dc.get_collectri(organism='human')

dc.run_mlm(adata, net=progeny, source='source', target='target',
           weight='weight', use_raw=False)
dc.run_ulm(adata, net=collectri, source='source', target='target',
           weight='weight', use_raw=False)

# pathway activity を donor x cell type で集計、age との correlation
```

## 10.4 既存 proteomics との照合

- **Lehallier et al. 2019 (Nat Med)**: 血漿プロテオームの加齢変化
- **CSF aging proteomics** (Mayo, MSBB, ROSMAP CSF)
- **Carlyle et al. 2017 (Nat Neurosci)**: human brain proteome aging

公開 supplement の Excel から差分タンパク質リストを抽出し、自分の hit と Fisher's exact test で overlap を検定。

```python
# scripts/09_figure7/02_proteomics_overlap.py
import pandas as pd
from scipy.stats import fisher_exact

hits = pd.read_csv('results/figure3/secretome_meta_analysis.csv')
hits_robust = hits.query('meta_padj < 0.05 and n_cohort >= 4')['gene'].unique()

for prot_db, path in [
    ('Lehallier2019', 'refs/Lehallier2019_aging_plasma.xlsx'),
    ('Carlyle2017', 'refs/Carlyle2017_brain_aging.xlsx'),
]:
    prot = pd.read_excel(path)
    overlap = set(hits_robust) & set(prot['gene'])
    # Fisher exact test
    # ...
```

## 10.5 Validation checkpoint (Figure 7)

- [ ] pathway 活性差が biological sense (SASP 上昇、IGF 低下、TGFβ 上昇)
- [ ] 自分の hit が外部 proteomics で >25% 程度 overlap (snRNA / proteomics の technical gap を考慮)
- [ ] SenMayo / SASP atlas との overlap が有意
- [ ] Summary 模式図が論文の中心メッセージを 1 枚で表現

→ `git tag figure7-v1` & snapshot.

---

# Section 11: 論文化フェーズ

## 11.1 Figure 完成後のチェック

- [ ] すべての figure が独立スクリプトから再現できる (中間ファイル削除 → snakemake -F で完走)
- [ ] `results/` の数値 CSV が figure と完全一致
- [ ] `docs/DECISIONS.md` に解析判断、`docs/DATA_PROVENANCE.md` に出典が完備
- [ ] 主要 panel を 300 dpi PDF と編集可能 SVG の両方で出力

## 11.2 投稿先候補（v3 改訂）

| 候補 | 理由 |
|---|---|
| **Aging Cell** | テーマ完全一致、データセット規模十分。**第一候補** |
| **Nature Aging** | multiome regulatory + 二重空間検証で要件を満たせる場合の挑戦 |
| **Cell Reports** | open access、機構解析志向 |
| **Cell Reports Medicine** | 統合解析で機構提示まで完結する場合 |
| **Neurobiology of Aging** | コンサバな選択 |
| **Nature Neuroscience** | 機能的 validation (オルガノイド or 動物実験) が追加で入る場合 |
| **eLife** | open review、reproducibility 評価される |
| **Scientific Data** | 統合データセット resource paper 路線 |

## 11.3 公開準備

- 全コード GitHub (投稿時 private、accept 後 public)
- 中間 h5ad を Zenodo に DOI 付きで deposit
- インタラクティブ可視化 (cellxgene VIP / ShinyCell) の URL を論文に記載
- 解析パイプラインを Snakemake workflow として公開（将来の再利用を考慮）

---

# Section 12: 補遺 — AI 連携ワークフロー

## 12.1 Claude Code 利用テンプレート

各 figure の作業前に `docs/CONTEXT.md` を最新化:

```markdown
# CONTEXT (last updated: YYYY-MM-DD)

## 現在の作業
Figure 3 — secretome 遺伝子の加齢 DE。
GSE268609, Zhou (3 シリーズ) は QC 済み。Franjic は QC 中。

## 直前の判断
- pct_counts_mt 閾値: GSE268609 5%, Zhou 5%, Franjic 8% (DECISIONS.md L42 参照)
- 統合: scVI, n_latent=30, batch=sample_id, covariate=[sex, pmi_hours]

## 今やりたいこと
DESeq2 pseudobulk を per-celltype × per-cohort で走らせ、結果を metafor で統合する R スクリプトを書く。

## 制約
- 入力 AnnData: processed/integrated/GSE268609_normal_scvi.h5ad (22 GB)
- celltype_l2 列を groupby
- donor 単位 pseudobulk
- 共変量: sex, pmi_hours
- target: age_years (連続)

## 期待する出力
results/figure3/de_age_per_celltype.csv
(columns: gene, celltype, cohort, log2FoldChange, lfcSE, pvalue, padj, in_secretome)
```

セッション先頭でこのファイルを読ませると、文脈を毎回説明し直さずに済む。

## 12.2 Codex CLI 利用テンプレート

```bash
codex "create scripts/05_figure3/03_volcano_plots.py that reads results/figure3/de_GSE268609_per_celltype.csv, makes a volcano plot per celltype, highlights genes in_secretome with red, saves PDFs to figures/figure3/volcano_<celltype>.pdf. Use matplotlib, no seaborn, dpi=300."
```

定型処理は `scripts/` 配下に numbered で揃えて、Codex CLI に「次の番号のスクリプトを以下の仕様で生成して」と頼むと scaffold が早い。

## 12.3 安全策

- 大規模スクリプトは **必ず dry-run / 小サンプルテスト** してから本番投入
- GPU 長時間ジョブは `tmux` 内で `nohup` 起動、`logs/<date>_<job>.log` に残す
- snakemake は `--use-conda` で env を rule 単位に分離
- すべての figure で乱数 seed を 42 に固定し、scripts に明記

## 12.4 トラブルシューティング Tips

| 症状 | 対処 |
|---|---|
| OOM (CPU) | anndata の `backed='r'` モードで読む（128GB あれば通常不要） |
| OOM (GPU) | scVI の `batch_size` を 512 → 256 に、`n_latent` を 30 → 20 に |
| scVI が dataset を分離しすぎる | `categorical_covariate_keys` から `dataset` を外し、`batch_key=sample_id` のみに / scANVI で label を入れる |
| CellChat が遅い | `triMean` を使い `min.cells=10` で thin out |
| rare cell type が安定しない | bootstrap 100 回の中央値報告 + `n_donors` を panel に注記 |
| Visium デコンが収束しない | `max_epochs` 増、`detection_alpha` 調整、reference の HVG を増やす |
| metafor で SE 不足 | 各 cohort で `lfcSE` を出力するよう DE スクリプトを修正 |

---

# Section 13: 査読対策の論理武装

| 想定質問 | 対応 |
|---|---|
| **GSE268609 を使うなら Lazarov 2024 の重複では？** | 同一 raw data からも切り口で独自性を出せる。①secretome 特化 lens、②multi-cohort meta-analysis、③二重空間検証、④cross-species ortholog、⑤multiome regulatory、の 5 軸で異なる結論を導出（Section 1.3）。 |
| **snRNA-seq では分泌タンパク質 mRNA を取り逃がすのでは？** | caveat として明記。proteomic data (CSF aging) との照合 (Figure 7B)、in situ hybridization (GSE248545、Figure 5D) で核外 mRNA も含めた発現を確認、と多段防御。 |
| **ヒトの neurogenesis 自体が controversial では？** | 本研究の主軸は **niche secretome** であり、神経新生の存在/非存在に依らない結論を立てる。Franjic データでも astrocyte/microglia 老化 secretome は再現することを示す（Figure 3）。 |
| **AD と aging を分離できているか？** | GSE268609 normal aging サブセットを主軸、AD/MCI は sensitivity analysis として supplementary figure で扱う。共変量に diagnosis を入れたモデルでも main hit が robust であることを示す。 |
| **年齢以外の confounder (sex, PMI, RIN, batch) はどう補正？** | 全 DE モデルに共変量として明示的に組み込む。`docs/DECISIONS.md` で設計判断を逐一記録。 |
| **rare cell type (NSC, neuroblast) の hit はどう信頼するか** | bootstrap stability 解析、cohort 横断再現性、空間検証 (Figure 5) の 3 段で防御。 |
| **マウス研究との連続性は？** | Cross-species panel (Figure 6) で保存されている軸を明示。ヒト特異性も同時に提示。 |
| **あなたの finding は SGZ 特有なのか pan-hippocampal aging か？** | **Figure 5 panel F (Regional specificity) で定量化**。GSE160189 で学習した region predictor で Discovery 細胞に DG/CA/SUB ラベルを付与し、aging secretome signature score の region 差を Spearman correlation で評価。GSE264692 Visium 空間でも再投影。 |

---

# Section 14: タイムライン & kickoff

## 14.1 タイムライン

| Phase | 期間目安 | GPU 重要度 |
|---|---|---|
| Phase 0 環境構築 | 2-3 日 | 低 |
| Phase 1 データ取得 | 1-2 週間 | 低 |
| Figure 1 (multiome 統合) | 2-3 週間 | 高 (scVI) |
| Figure 2 lineage | 1-2 週間 | 中 |
| Figure 3 meta-analysis (5+ cohort) | 3-4 週間 | 低 (DE は CPU) |
| Figure 4 communication + regulatory | 2-3 週間 | 中 (chromVAR) |
| Figure 5 二重空間検証 | 3-4 週間 | 高 (cell2location) |
| Figure 6 cross-species | 1-2 週間 | 中 |
| Figure 7 integration | 1-2 週間 | 低 |
| 論文化 | 4-6 週間 | — |

合計 **5-7 ヶ月**（1 人での進行を想定）。

## 14.2 最初の 3 日 kickoff チェックリスト (v4)

- [ ] Ubuntu / nvidia driver / CUDA 動作確認 (`nvidia-smi`)
- [ ] mamba で `sgz` (Python) と `sgz_r` (R) 環境構築、PyTorch GPU 動作確認
- [ ] `$PROJ` ディレクトリ作成、git init、.gitignore、DVC init
- [ ] HDD への cron rsync 設定 + `snapshot.sh` テスト
- [ ] **GSE268609 の最新 metadata を NCBI GEO で確認** (Lazarov 論文の data availability)
- [ ] **GSE264624 と GSE264692 のペア関係を Tran 2025 論文 supplement で確認** (どちらがどのドナーの snRNA / Visium かをマッピング)
- [ ] **GSE160189 の region label (anterior/posterior, DG/CA/SUB) の品質を確認** — region predictor として使えるか
- [ ] **GSE325391 の status を確認** (Public か embargoed か)
- [ ] Lazarov 2024 論文 (PMC11160907) の Methods と Code Availability を精読 → 著者公開コードがあれば取得
- [ ] GSE268609 のサンプル 1 件をダウンロードし、muon で multiome load 成功確認
- [ ] Zhou シリーズ (GSE198323) のサンプル 1 件で snRNA load 成功確認
- [ ] **GSE264692 のサンプル 1 件を `sc.read_visium()` で load 成功確認** (Visium 標準ディレクトリ構造になっているかが鍵)
- [ ] **GSE160189 のサンプル 1 件で region label が AnnData の obs にあるか確認**
- [ ] `refs/metadata_master.csv` 空テンプレ + `paired_with` カラム + GEOparse で raw metadata 抽出パイプライン構築
- [ ] secretome reference (HPA / UniProt / MatrisomeDB) ダウンロード
- [ ] NicheNet v2 priors ダウンロード
- [ ] `docs/CONTEXT.md` を Claude Code に渡せる状態に整備

## 14.3 最重要原則の再掲

**Figure 1 が出る前に Figure 7 を縛らない**。

この protocol は研究の骨子だが、**実際には Figure 1 で見えた biology が次の analytic decision を決める**。毎 figure 後に `DECISIONS.md` を更新し、必要なら計画を修正してください。事前計画の固定化より、**実データを見て判断する原則**を貫くことが論文の質を決めます。

---

**完全版 v3 終わり**

各 phase で実装が止まったり、Figure 1 の中間結果が出た時点で具体相談したいことが出たら、Section 12.1 の `CONTEXT.md` 形式でまとめて持ってきてください。文脈を取り直さず直接議論に入れます。

Good hunting.
