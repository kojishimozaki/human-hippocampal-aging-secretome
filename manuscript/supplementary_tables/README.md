# Supplementary Tables S1, S2, S(R1)

原稿 `manuscript/manuscript.md` の「Supplementary tables」節に対応する提出用テーブル。
**新たな計算は一切していない。** すべてコミット済みの結果 CSV から組み立てただけであり、
各セルは下記の出典ファイルに1対1で遡れる。

再生成: `make supplementary-tables`（生成器 `scripts/16_supplementary_tables/build_supplementary_tables.py`）

| ファイル | 行×列 | 出典 |
|---|---|---|
| `Table_S1.csv` | 10 × 21 | `results/audit_reruns/de_calibration/donor_label_permutation_summary.csv` |
| `Table_S2.csv` | 10 × 21 | `results/de/atac_gene_vs_peak_level_matched_null.csv` ＋ `atac_linkedpeak_enrichment.csv` |
| `Table_SR1.csv` | 62 × 21 | `results/regulatory/motif_tf_rna_concordance_{HAvYA_primary,ADvHA_projection}_pairs.csv` |
| `Supplementary_Tables.xlsx` | README＋3シート | 上記3表を1ブックに |

> ファイル名を `Table_SR1` としたのは、括弧が Makefile とシェルを壊すため（CLAUDE.md 鉄則10）。
> 原稿中の表番号は **Table S(R1)** のままである。

---

## Table S1 — ドナーラベル置換によるDE件数の較正

コホート × ニッチ細胞型ごとに、ドナー疑似バルクモデルをラベル置換して100回フルに再フィットした結果。

| 列 | 意味 |
|---|---|
| `status` | `OK` / `GATE_FAIL`（最小ドナー数ゲート不通過） |
| `n_donors_YA` / `n_donors_HA` | 当該細胞型で実際に寄与したドナー数 |
| `n_distinct_label_splits` | ラベル割当の総数。GSE268609 の 8対9 は **24,310 通り**で、100置換は汲み尽くさない |
| `obs_secretome_fdr01` | 観測されたセクレトーム件数（補正 p < 0.1） |
| `null_median_*` / `null_p95_*` / `null_max_*` | 置換帰無分布の中央値・95パーセンタイル・最大 |
| `emp_p_*` | 経験的 p。**下限は 1/101 = 0.0099** |
| `null_frac_ge1_*` / `null_frac_ge10_*` | 帰無下で BH 有意遺伝子が1件以上／10件以上出る割合 |

**GSE278576 の内皮は `GATE_FAIL`（YA 6 対 HA 2）として行を残してある**——省略すると欠測が見えなくなるため。

## Table S2 — アクセシビリティ・マッチ帰無：ピーク単位 対 遺伝子単位

細胞型 × 遺伝子集合（凍結60ヒット連結／全セクレトーム連結）ごとに、同一の帰無を
**ピークを単位**にした場合と**遺伝子を単位**にした場合で並べたもの。

| 列 | 意味 |
|---|---|
| `peaks_per_gene` | 遺伝子あたりピーク数（**3.17–6.15**）。遺伝子集約が必要な理由 |
| `top_gene` / `top_gene_share` | 単一遺伝子が集合のピークに占める割合（**2.5%–68.4%**。Astro=TAFA1 33.1%、OPC=FGF13 68.4%） |
| `obs_peak_level` / `p_peak_level` | ピーク単位の観測値と経験的 p（**比較用。推論単位ではない**） |
| `obs_gene_level` / `p_gene_level` | 遺伝子単位（**こちらが推論単位**） |
| `p_gene_signperm` / `p_gene_LOO` | 遺伝子水準の符号置換 ／ leave-one-gene-out |
| `null_sd_ratio_gene_over_peak` | 遺伝子帰無 ÷ ピーク帰無の標準偏差比（**1.23–1.85**＝ピーク単位帰無の過小分散） |
| `n_peaks_fdr_significant` | per-peak FDR を通過したピーク数。**10行すべてで 0** |
| `prespecified_gonogo_pass` | 事前規定 go/no-go の合否。**10行すべて False**（第3成分が構成上到達不能） |

## Table S(R1) — モチーフ活性 対 転写因子RNA の一致（TF × 細胞型ごと）

chromVAR で加齢動的（名目 p < 0.05）だったモチーフの凍結前景31ペアについて、
**HA対YA（一次）と AD対HA（記述的投影）を `contrast` 列で積み上げた**もの（31 × 2 = 62行）。

| 列 | 意味 |
|---|---|
| `contrast` | `HAvYA_primary` / `ADvHA_projection` |
| `dz` / `cv_p` / `cv_padj` | chromVAR のモチーフ活性変化と名目 p・細胞型内 BH q |
| `rna_lfc` / `rna_se` / `rna_p` / `rna_q` | 当該 TF 自身の RNA。**元の細胞型別DE宇宙から引いており、前景で再計算していない** |
| `testable` | 当該細胞型で TF が頑健に発現しているか。**31件中9件が `False`** |
| `concordant` | 符号一致（RNA有意性とは独立に判定） |
| `category` | `A`（一致かつFDR有意、**0件**）／`B_concordant_q>=0.05` 11件／`C_discordant` 11件／`D_not_testable` 9件 |

検定可能22件のうち一致は11件（50%）で、マッチ背景の 52.4% を上回らない。

---

## 原稿の記述との照合（生成時に検算済み）

| 原稿の主張 | 実データ |
|---|---|
| GSE278576 内皮は GATE_FAIL、6 YA 対 2 HA | ✓ 一致 |
| 8対9 は 24,310 通りのラベル割当 | ✓ 一致 |
| 経験的 p の下限 1/101 = 0.0099 | ✓ 最小値 0.009901 |
| peaks per gene 3.2–6.2 | ✓ 3.17–6.15 |
| top_gene_share 2.5%–68%、TAFA1 33%、FGF13 68% | ✓ 0.025–0.684、33.1%、68.4% |
| null_sd_ratio 1.23–1.85 | ✓ 一致 |
| per-peak FDR は10行すべて0 | ✓ 一致 |
| 31 TF 中 9 が NOT_TESTABLE | ✓ 一致 |
| 検定可能22件中11件が一致 | ✓ 一致 |
