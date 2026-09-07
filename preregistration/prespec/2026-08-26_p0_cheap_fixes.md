# 事前固定した判定基準 — P0-1 が残した STILL-OPEN のうち「安い3件」

作成: 2026-08-26。**計算を1本も走らせる前に書き、この時点で commit する**（CLAUDE.md 鉄則5）。
由来: `audit_log/2026-08-25_p0_1_unread_outputs/RESOLUTION.md` §7 の推奨。

対象は `F-5-006` / `F-5-007` / `F-3b-012`。いずれも解析の再設計ではなく、
監査キット側スクリプトの小さな欠陥で出力が得られていなかったもの。

---

## 0. 共通規則

- 閾値は凍結シグネチャと同じ `padj < 0.1` を使う。**別の閾値で再判定しない**（鉄則7）。
- 生物学的反復の単位は **donor**。cell 単位の p を推論に使わない（鉄則1）。
- 数値はすべて出力 CSV の `(file, row, column)` に遡って報告する（鉄則2）。
- 修正は「数値を動かさない実装の修正」に限る。**手法は置き換えない**（鉄則4）。
  置き換えが必要だと判明したら STOP して報告する。
- 判定は **CLOSED-OK / CLOSED-CHANGE / STILL-OPEN** の3値（P0-1 と同じ定義）。

---

## 1. `F-5-006` — LFC 収縮

### 事前に確定させた事実（実測済み、2026-08-26）

- `~grp` の design 行列の係数は **`grp[T.YA]`** の1本のみ（HA が参照レベル）。
  スクリプトが要求している `grp[T.HA]` は存在せず KeyError になる。
  一巡目ログ `audit/pass6/_logs/prior_pass5_05_covariates.log` に
  `The available LFC coeffs are Index(['grp[T.YA]'])` として5回記録されている。
- 合成データ（真の DE 80 遺伝子、UP +2 / DOWN −2）で検証した結果:
  `lfc_shrink(coeff="grp[T.YA]")` は **YA vs HA** を返す。
  `corr(MLE_{HA vs YA}, shrunk) = −1.0000`、`corr(MLE, −shrunk) = +1.0000`。
  → **符号を反転して初めて `lfc_MLE` 列（contrast HA vs YA）と同じ向きになる。**

### 判定基準

`stats/pass5_lfc_shrinkage_frozen60.csv` の 60 行について:

1. `lfc_shrunk` の NaN が **0/60** になること。ならなければ **STILL-OPEN** のまま。
2. **符号の健全性ゲート**: `sign(lfc_shrunk) == sign(lfc_MLE)` が **60/60**。
   1件でも外れたら符号反転の向きを誤っているので **STOP して報告する**（自動で先へ進まない）。
3. 報告する量: frozen 60 の `median |lfc_MLE|` と `median |lfc_shrunk|`、
   `|lfc_shrunk| > 1` の件数、`shrink_factor = |shrunk|/|MLE|` の分布。
4. **SERPINE1（Micro, baseMean 10.686, lfc_MLE +4.125, lfcSE 1.038, padj 0.0299）**について:
   - `|lfc_shrunk| < 1.0` → **CLOSED-CHANGE**。抄録・本文の見出しから "+4.1" を外し、
     収縮値を主語にする。
   - `1.0 ≤ |lfc_shrunk|` かつ `shrink_factor < 0.5` → **CLOSED-CHANGE**。両方を併記する。
   - `shrink_factor ≥ 0.5` → **CLOSED-OK**。収縮しても主張は保たれる旨を記す。

**注**: 収縮値は `~grp`（一次 design）に対してのみ計算する。他の design には拡張しない。

---

## 2. `F-5-007` — PMI を含む同時調整

### 事前に確定させた事実（実測済み、2026-08-26）

- 一巡目で `~pmi+grp` / `~arm+pmi+grp` / `~sex+arm+pmi+grp` が
  **5細胞型すべてで SKIPPED**（ログに 15 行）。`summary.csv` にも1行も無い。
- 原因: `pmi_src = refs/gse268609_group_pmi_age.csv` は **群単位の要約**
  （列 `Group,n_donor,pmi_mean,pmi_sd,...`）で、donor 列も素の `pmi` 列も無い。
  そのため `pmi_map = {}` になり `md["pmi"]` が全 NaN → SKIP。
- **per-donor PMI は、既に `--sexcsv` で渡している
  `audit/confounding/GSE268609_recovered_sex.csv` の `pmi` 列に存在する**
  （YA/HA 17 donor、NaN 0 件、`donor == SampleNumber`）。
- その per-donor 値を群集計すると `refs/gse268609_group_pmi_age.csv` を**完全に再現**する
  （HA n=9 mean 7.67 sd 2.08 min 4.4 max 10.0 / YA n=8 mean 5.84 sd 0.96 min 4.8 max 7.0）。
  → **同一の SOFT 由来 PMI であり、新たな導出ではない**（鉄則4 を満たす）。

### 判定基準

1. **回帰ゲート**: 既に出ている 15 行（`~grp` / `~arm+grp` / `~logncell+grp` × 5細胞型）は
   `n_donor` `n_evaluable` `sign_retained` `padj01_retained` `lfc_corr` が
   **一巡目と完全一致**すること。1つでもずれたら PMI の追加が既存 fit を汚しているので **STOP**。
2. 新規3 design が 5細胞型で SKIPPED にならないこと。SKIP が残るなら理由を明記して STILL-OPEN。
3. **原稿 Limitations 7 に載せる数値は `~arm + pmi + grp` の `padj01_retained` の5細胞型合計**
   とあらかじめ決める。`~pmi+grp` と `~sex+arm+pmi+grp` は併記するが、
   **結果を見てから「都合のよい design」を主役に選び直さない**（鉄則5・鉄則7）。
4. 判定:
   - `~arm+pmi+grp` の合計が既報の `~arm+grp`（42/60）と比べて
     **±20% 以内**（すなわち 34–50/60）→ **CLOSED-CHANGE**（数値を Limitations に追加するだけ）
   - それを外れる → **CLOSED-CHANGE** かつ「PMI 調整で結論が動く」ことを本文で開示する
   - どちらにせよ **符号保持（sign_retained）を必ず併記**する。

### 適用範囲の制限（プロジェクト側の既知の注意を引き継ぐ）

`scripts/03_de/04_aging_robustness.py:L75-77` が明記するとおり、SOFT の PMI/age は
**YA/HA についてのみ群整合が検証済み**（8/8・9/9）。AD/MCI/SA には
RNA+ATAC の GSM 間で不整合がある。`pass5_05` は YA/HA しか使わないので該当しないが、
**この結果を疾患軸の PMI 調整に流用してはならない**旨を RESOLUTION に明記する。

---

## 3. `F-3b-012` — ATAC 側のドナーあたり最低細胞数

### 事前に確定させた事実

- `pass3b_03_atac_multiplet_upperbound.py` の感度ループは
  `rule ∈ {none, p99, p99.5, p99.9, 5x_median}` のみで、
  **ドナーあたり最低細胞数の門を課した variant が存在しない**。
  現状はある細胞型で 1 細胞しか持たないドナーも pseudobulk 列として等価に寄与する。
- 既存のドナー門は `if ha.sum() < 3 or ya.sum() < 3: continue`（ドナー数のみ）。

### 判定基準

`min_cells ∈ {0, 10}` の第2軸を追加し、出力に `min_cells` 列を足す。

1. **回帰ゲート**: `min_cells = 0` の 25 行が、既存の
   `pass3b_atac_multiplet_sensitivity.csv` と **すべての列で完全一致**すること。
   ずれたら実装が既存挙動を変えているので **STOP**。
2. 判定（`min_cells = 10` の行について）:
   - `concordance` と `A2_obs_mean_oriented` が、既に報告されている multiplet 規則の
     振れ幅（Fig 4A 40–43/60、A2 obs はベースラインの 87–101%）の**内側**に収まる
     → **CLOSED-OK**（門を課しても結論不変。Methods にドナー包含基準を1文）
   - 外side に出る → **CLOSED-CHANGE**（Supplementary に感度表を載せ、本文に注記）
   - 門によって `ha.sum() < 3 or ya.sum() < 3` で落ちる細胞型が出たら、
     その細胞型は **検定不能**と明記する（「問題なし」に丸めない）。

---

## 4. やらないと事前に決めたこと

- 原稿（`manuscript/manuscript.md` および `_ja.md`）・図・図スクリプトは**編集しない**。
  成果物は数値と変更案までとする（P0-1 と同じ扱い）。
- P0-2〜P0-5、および P0-1 が STILL-OPEN に残した他の3件
  （`F-3b-009` 行順、`F-1-021`/`F-8-025` ドナー重複、`F-3b-002` の OPC 肢）には着手しない。
- 介入1・7（疎行列化）の独立検証は行わない。
- `raw/` `processed/` には書き込まない（読み取り専用のまま）。

## 5. 出力先

`audit_log/2026-08-26_p0_cheap_fixes/RESOLUTION.md`
キット側スクリプトへの修正は `~/audit-runs/20260821-225834/compute/_PATCHES.md` §8 に開示する。
