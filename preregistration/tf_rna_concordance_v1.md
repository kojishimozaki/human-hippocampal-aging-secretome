# PRE-REGISTRATION — Motif↔TF-RNA same-direction cross-modal concordance (HAvYA primary; ADvHA projection; GSE278576 RNA-direction robustness)

**Status: LOCKED pre-analysis. Written and committed BEFORE any analysis code is written or run.**
日付 / Date: 2026-06-19  ·  著者 / Author: kj.shimozaki  ·  設計 / Design: locked via `/grill-me` (7 branches, all PI-confirmed)
登録コミット / Registration commit: _TBD — step H.3 で commit 後に追記_  ·  annotated tag: `tf-rna-concordance-prereg`

> 本書は Fig 3B（chromVAR aging motif-activity, HAvYA, nominal）でアクセシビリティが加齢変動する
> motif に対応する **TF 遺伝子の RNA が、同方向に加齢変動するか**を、適格な経験背景に対する
> **same-direction cross-modal concordance の濃縮**として記述的に検定する計画である。
> これは既存の三角測量（PHASE_REPORT Layers A–D = motif濃縮 / NicheNet / footprint / L-R、
> convergence の3層 = 濃縮・chromVAR方向・標的数）が**一度も検定していない "TF発現脚"** にあたる。

---

## 0. 解釈の天井（branch A — 最上流の固定）
- 位置づけ = **記述的 / 補完的（descriptive / supporting）**。cis-TF 機構の結論は動かさない。
- cis-TF 機構は既に **footprint（高検出力, n=8/9）・chromVAR・motif濃縮すべてで FDR-null** ＝
  「honest null / 機構未解決」に reframe 済み（`results/regulatory/PHASE_REPORT.md`）。
- **陽性が出ても cis 機構を復活させない**（鉄則9：reframe 済み主張の silent revival 禁止）。
  footprint = 「結合は加齢で増えない（nominal はむしろ浅くなる）」が高検出力の上位読み出しであり、
  本解析（TF mRNA 方向）が陽性でも、それは「TF発現が motif活性変化と同方向に並ぶ」という
  **記述**にとどまり、TF活性・結合・占有・因果の証明にはしない。

## 1. Headline question（HAvYA primary, GSE268609 単独）
> 「健康な加齢（HA vs YA）に伴う motif-associated chromatin accessibility 変化は、
>  対応する TF の RNA の **同方向変化**を、**適格な matched empirical background より高頻度に**伴うか。」

- primary inference は **GSE268609 内の HAvYA same-direction cross-modal concordance enrichment のみ**。
- 用語は **same-direction cross-modal concordance の濃縮**で統一。"coherence" / "TF activity" /
  機序的 coupling は用いない。

### 述べない主張（locked non-claims）
TF タンパク量・結合・活性化・motif 占有・因果が証明された、とは一切述べない。activator/repressor
分類だけで符号関係を機械的に反転させない。caveat：repressor / pioneer factor / motif-family
ambiguity / 翻訳後制御。

---

## 2. 入力データ（凍結 — SHA256 で固定）
| role | file | SHA256 |
|---|---|---|
| chromVAR (HAvYA & ADvHA, donor-pseudobulk t, per-celltype BH) | `results/de/chromvar_motifs_strengthened.csv` | `8cbe09cb…52609` |
| RNA DE HAvYA (donor-pseudobulk, `~[sex+]grp`, HA vs YA) | `results/de/de_GSE268609_per_celltype.csv` | `950173e4…67c2b` |
| RNA DE ADvHA (donor-pseudobulk, AD vs HA) | `results/de/de_GSE268609_ADvHA.csv` | `62043233…760e0` |
| RNA DE GSE278576 (sex-adjusted, old vs young) | `results/de/de_GSE278576_sexadj_per_celltype.csv` | `206e4da6…b5195` |

完全な SHA256 は `…_run_manifest.json` に記録。モデル差（268609 は `~[sex+]grp`／278576 は
sex-adjusted、donor 数も異なる）は出力表に明記する。

### 符号規約（コードで検証必須）
両モダリティ・両コホートとも **「高齢 − 若齢」が正（aging-UP = +）**。
検証済み（pre-analysis）: RNA log2FC>0 が aging-UP（IL15/Micro +3.11, TNC/Astro +1.35, COL21A1 +1.91 @268609;
IL15/Micro +3.23 @278576）、chromVAR は `contrast=HAvYA`（dz>0 = 老齢で高アクセス）。
スクリプトは join 前に **assert で符号規約を再検証**する。

---

## 3. Foreground 定義（branch B — RNA を見る前に凍結）
- **単位 = TF × cell type**。
- foreground = chromVAR **HAvYA** で **nominal p < 0.05** の (TF × cell type)。
  （FDR<0.1 は全細胞型で 0。3B は p昇順 top-N 表示の viz であり、ここでは閾値で原理的に凍結。）
- motif→TF は §4 の direct mapping。RNA 発現閾値を満たし age coefficient が推定可能な
  foreground のみ headline 分母に入れる（= **22**）。残り 9 は **NOT_TESTABLE_LOW_EXPRESSION**
  として暗黙除外せず別掲。
- 同一 (TF gene, cell type) を複数 motif が指す場合は canonical 1 件に collapse（代表 = chromVAR
  nominal-p 最小、同点は非二量体/非variant 優先）。**本 foreground は衝突 0（確認済み）**。
- **ADvHA・GSE278576 の結果に基づく foreground の追加・削除・再選択は一切行わない。**

### 凍結 foreground 一覧（31 TF×ct; dz/p は chromVAR HAvYA; gene = direct map）
| TF (motif) | celltype | dz | p | mapped gene | test@268609 | test@278576 |
|---|---|---|---|---|---|---|
| ETV1 | Astro | -0.251 | 0.0013 | ETV1 | YES | YES |
| ELF3 | Astro | -0.218 | 0.0027 | ELF3 | **NO(low-expr)** | YES |
| GABPA | Astro | -0.260 | 0.0037 | GABPA | YES | YES |
| EHF | Astro | -0.262 | 0.0039 | EHF | YES | YES |
| SPIB | Astro | -0.239 | 0.0047 | SPIB | YES | low-expr |
| ELF1 | Astro | -0.238 | 0.0055 | ELF1 | YES | YES |
| ETV4 | Astro | -0.227 | 0.0057 | ETV4 | YES | YES |
| ETV2 | Astro | -0.185 | 0.0258 | ETV2 | **NO(low-expr)** | low-expr |
| ETV6 | Astro | -0.223 | 0.0261 | ETV6 | YES | YES |
| ZNF143 | Astro | -0.204 | 0.0339 | ZNF143 | YES | YES |
| ELF5 | Astro | -0.184 | 0.0422 | ELF5 | **NO(low-expr)** | low-expr |
| RFX1 | Astro | +0.367 | 0.0470 | RFX1 | YES | YES |
| SOX9 | Endo | -0.325 | 0.0135 | SOX9 | YES | celltype-mismatch |
| SMAD5 | Endo | -0.203 | 0.0142 | SMAD5 | YES | celltype-mismatch |
| TCF7L1 | Endo | -0.381 | 0.0222 | TCF7L1 | YES | celltype-mismatch |
| TCF7L2 | Endo | -0.393 | 0.0243 | TCF7L2 | YES | celltype-mismatch |
| BARHL2 | Endo | -0.295 | 0.0273 | BARHL2 | **NO(low-expr)** | celltype-mismatch |
| SOX8 | Endo | -0.259 | 0.0283 | SOX8 | YES | celltype-mismatch |
| EGR3 | Endo | -0.365 | 0.0314 | EGR3 | YES | celltype-mismatch |
| NR2C1 | Endo | -0.196 | 0.0473 | NR2C1 | YES | celltype-mismatch |
| E2F3 | Micro | +0.174 | 0.0334 | E2F3 | YES | YES |
| TP63 | Micro | +0.163 | 0.0344 | TP63 | YES | YES |
| TP73 | Micro | +0.215 | 0.0417 | TP73 | **NO(low-expr)** | low-expr |
| HOXD9 | Micro | +0.249 | 0.0455 | HOXD9 | **NO(low-expr)** | low-expr |
| ZBTB18 | OPC | -0.342 | 0.0116 | ZBTB18 | YES | YES |
| TWIST1 | OPC | -0.269 | 0.0245 | TWIST1 | **NO(low-expr)** | low-expr |
| MYF5 | OPC | -0.447 | 0.0295 | MYF5 | **NO(low-expr)** | low-expr |
| ELK1 | Oligo | -0.272 | 0.0150 | ELK1 | YES | YES |
| FEV | Oligo | -0.219 | 0.0274 | FEV | **NO(low-expr)** | low-expr |
| ELK3 | Oligo | -0.240 | 0.0276 | ELK3 | YES | YES |
| ETV3 | Oligo | -0.224 | 0.0388 | ETV3 | YES | YES |

集計: **268609 testable = 22 / 31**（headline 分母）, NOT_TESTABLE_LOW_EXPRESSION = 9。
278576 robustness: YES=15, celltype-mismatch(Endo)=8, low-expr=8, **両コホート交差 testable = 14**。
構造的特徴（記録のみ、選択には未使用）: 31 中 26 が dz<0（加齢で motif 活性*喪失*; ETS/SOX/TCF/bHLH）、
正は RFX1/Astro と Micro の E2F3/TP63/TP73/HOXD9 のみ。

---

## 4. motif→TF 写像（branch C — direct primary / family は分離した sensitivity）
- **Primary = 直接名一致のみ**。motif 名 → 同名 gene。二量体は構成遺伝子へ分割、`(var.N)` は除去。
  （本 foreground は全て single 名のため分割不要。規則は背景・将来のため明文化。）
- 同 celltype の RNA DE に gene が無い場合は **NOT_TESTABLE_LOW_EXPRESSION** として理由付きで残す
  （暗黙 drop 禁止）。
- **Family sensitivity は完全に別レイヤー・別表・別検定 family**（§7・§9）。primary の数値に混ぜない。
  手動 family グルーピング（ETS / bHLH / STAT / IRF / FOX / SOX / TCF-LEF / AP-1 / KLF-SP …、外部 DB 不要）で
  「発現している family パラログのいずれかが motif 方向に動くか」を探索的に報告。

---

## 5. Concordance の操作的定義（branch D-1）
同一 contrast・同一符号規約で:
- **concordant**: `sign(RNA age log2FC) == sign(chromVAR age dz)`
- **discordant**: 符号不一致
- **RNA padj は concordance 判定の条件にしない**（非有意を「変化なし」と扱わず、方向と統計的支持を分離）。
- 各 pair を 4 カテゴリで併記: **A** concordant & RNA BH q<0.05 / **B** concordant & q≥0.05 /
  **C** discordant / **D** not testable。
- 符号関係は activator 仮定（同符号=concordant）で機械的に評価し、activator/repressor で反転しない。
  repressor / pioneer / family ambiguity / PTM を caveat に明記。

---

## 6. Empirical background（branch D-2 — matched/stratified が primary）
**Primary 背景** = 以下を全て満たす **非 foreground** TF×cell type:
同一 mapping ルール / chromVAR age effect 推定可能 / RNA が同一発現閾値で age coef 推定可能 /
同一 celltype 解析・共変量モデル / foreground 自身は背景から除外。

**Primary 検定 = matched / stratified empirical background に対する foreground の concordance 濃縮。**
層化（match）変数（最低限）:
1. cell type
2. chromVAR age effect の符号（up / down）
3. RNA baseMean（celltype 内 log(baseMean) tertile; quartile を sensitivity）

**stratified permutation**: foreground の (celltype × chromVAR符号 × expr-tertile) 同時分布を固定し、
背景から各層 foreground と同数を抽出 → concordance率 / risk difference を計算 → 反復。
`seed=42`, `n=10000`, `p=(1+#{perm stat ≥ obs})/(n+1)`。背景層が foreground 層より薄い時のみ復元抽出し
**ログに記録**。反復数・seed・p 計算法を manifest に保存。

**併記する量**: foreground concordance 数 / testable 数、率と 95% CI、matched background 率、
absolute risk difference、odds ratio と 95% CI、背景較正した exact conditional または stratified
permutation の p。

### Sensitivity（primary とは分離）
1. 全適格背景 pool の 2×2 Fisher exact
2. 50% 帰無の単純 binomial は **参考値のみ**（headline にしない）
3. `|chromVAR effect|`（または順位 bin）を合わせた背景 ← foreground は chromVAR-p で選抜のため重要
4. motif family を合わせた背景
5. RNA log2FC と chromVAR effect の Spearman 相関 / 符号付き連続 score

**pool-Fisher と matched-background が食い違う場合は matched を優先**し、差を celltype 構成・発現量
構成で説明する。

---

## 7. ADvHA secondary（branch F — α: descriptive projection、3群再実行なし）
- **HAvYA で凍結した同一 foreground (motif–TF×celltype)** をそのまま ADvHA で評価する **projection**。
  AD 側で foreground を選び直さない（selection bias 回避）。
- 「健康加齢で同定した motif–TF RNA 関係が AD で維持 / 追加変化 / 反転を示すか」を調べる。
- **HAvYA と ADvHA は pool しない**。検定・BH・表・図・記述を contrast ごとに**完全分離**。
- coherence enrichment は contrast ごとに独立計算（HAvYA = primary, ADvHA = secondary projection）。
  **contrast 間の formal な差の検定（3群 interaction / second-difference）は行わない。**
- **HA が共有参照**ゆえ 2 contrast を独立標本扱いせず、**p値・q値を直接比較しない**。
- cross-contrast 記述カテゴリ（per pair）:
  1 age-only concordant / 2 disease-only concordant / 3 shared same-direction /
  4 shared opposite-direction / 5 age-concordant & disease-discordant / 6 not testable。
- **用語制限**（効果量・符号・区間で記述; 有意/非有意の切替で断定しない）:
  - 許可: *same-direction continuation* / *opposite-direction secondary shift* /
    *no clear additional shift* / *inconclusive*
  - 「維持/増強/消失/反転」は群平均・効果量から明確な時のみ。
  - ADvHA 非有意だけで「維持」としない／同符号だけで「増強」としない／逆符号だけで
    「消失・正常化」としない。opposite sign は原則 *trajectory reversal / opposite-direction shift* と書き、
    AD が YA を越えて反転した意味では使わない。
  - **age-specific / AD-specific / significantly amplified|reversed in AD とは書かない**（formal 検定なし）。
- 3群モデルは将来 cross-contrast difference 自体を独立仮説として検証する際の追加候補として記録。

---

## 8. GSE278576 RNA-direction robustness（branch G — G-b; **replication ではない**）
- 厳密に **external RNA-direction robustness check**。joint cross-modal concordance の replication /
  validation / external reproduction とは**呼ばない**（278576 に対応する外部 chromVAR が無い）。
- 対象 = **HAvYA で凍結した foreground のみ**。278576 側で TF/celltype を結果に基づき追加・再選択しない。
  family パラログは本 robustness の primary 対象に含めない。
- **cell type 1:1 対応（凍結）**: Astro↔Astro, Micro↔Micro, OPC↔OPC, Oligo↔Oligo。
  **Endo は 278576 に対応無し → NOT_TESTABLE_CELLTYPE_MISMATCH**（曖昧統合・1対多・事後 regroup 禁止）。
- 278576 も「高齢−若齢」符号に統一し、268609 と方向規約一致を**コードで検証**。モデル差を表に明記。
- 判定: `sign(268609 RNA age effect) == sign(278576 RNA age effect)` = cross-cohort RNA-direction agreement。
  RNA padj<0.05 を一致の必要条件にしない。
- per pair 出力: 両コホートの RNA effect / SE or 95%CI / exact p / BH q / baseMean / agreement / testability。
- set-level（記述のみ）: testable pair 数、agreeing 数、agreement率 + exact 95%CI、effect-size の Spearman、
  scatter 用データ。50% 基準の単純有意性主張はしない。**両コホート交差 testable = 14**（事前確認値）。
- （任意・別 family）適格 TF×celltype から celltype・発現量を合わせた背景で foreground の cross-cohort
  方向一致率が一般的 RNA 方向安定性より高いかを exploratory 評価。既存ファイルだけで安定に作れない場合は
  検定せず agreement率・CI・効果量相関の記述に留める。
- 許可表現: *RNA component showed cross-cohort directional robustness* /
  *direction of age-associated TF RNA changes was broadly preserved* / *limited or mixed RNA-direction robustness*。
  禁止: external replication of cross-modal concordance / replicated motif–TF mechanism /
  independently validated chromatin–RNA coupling / externally confirmed concordance enrichment。
- 結果が弱くても省略せず報告し、268609 primary・frozen foreground は変更しない。

---

## 9. 多重検定（branch E）
- **事前固定の single global concordance test を単一 headline** として扱う。
- per-TF RNA DE の BH q は **元の contrast の全 DE testing universe** に基づく値を使う。
  **結果を見た後に 22 TF だけで BH をやり直して primary q としない。**
- 別 family として独立に BH 補正:（i）HAvYA primary concordance（ii）ADvHA projection concordance
  （iii）family sensitivity（iv）GSE278576 robustness（v）cross-cohort exploratory 背景。
- **2 contrast を合算して BH を再計算しない。** cross-contrast / cross-cohort 比較を複数 TF で行う場合は
  その比較群に別途 BH。

---

## 10. 解釈の限界（locked）
結論は「foreground motif では対応 TF RNA の同方向変化が経験背景より**濃縮していた／濃縮の証拠は
得られなかった**」までとする。同方向でも TF タンパク量・結合・活性化・motif 占有・因果の証明にしない。
NULL も結果（検出力を添えて "no evidence for"）。limitation に残す:（i）joint concordance は単一コホート評価、
（ii）外部で chromVAR leg 未再構築、（iii）278576 は RNA component のみの robustness、（iv）n=8/9 の検出力、
（v）motif↔TF の family ambiguity（9/31 は命名 TF がその細胞で非発現）。

---

## 11. 出力（branch I — contrast 完全分離がファイル名で判る）
`results/regulatory/` 配下:
- `motif_tf_rna_concordance_frozen_foreground.csv`（§3 の凍結集合 + mapping/testability）
- `motif_tf_rna_concordance_mapping_audit.csv`（motif→TF 写像の全件監査）
- `motif_tf_rna_concordance_HAvYA_primary_pairs.csv`
- `motif_tf_rna_concordance_HAvYA_primary_enrichment.csv`
- `motif_tf_rna_concordance_ADvHA_projection_pairs.csv`
- `motif_tf_rna_concordance_ADvHA_projection_enrichment.csv`
- `motif_tf_rna_concordance_HAvYA_family_sensitivity.csv`
- `motif_tf_rna_concordance_ADvHA_family_sensitivity.csv`
- `motif_tf_rna_concordance_crosscontrast_descriptive.csv`
- `motif_tf_rna_concordance_GSE278576_RNA_robustness.csv`
- `motif_tf_rna_concordance_run_manifest.json`

### run_manifest.json（最低限）
git commit hash / 実行日時 / Python・package versions / seed / 入力パス + SHA256 /
foreground 件数 / motif→TF mapping version / cell type mapping version / 発現閾値 /
BH family の定義 / 出力ファイル一覧。

---

## 12. Process（branch H — 監査証跡）
1. 本書（A–G + H + I）を commit。**解析コード変更・実行より前**。
2. remote へ push（off-machine 耐久性）。**※現状 remote 未設定 → PI 判断待ち（下記 Open）。**
3. commit hash・日時・入力一覧を本書に追記（§登録コミット）。
4. 解析実行（`scripts/05_regulatory/06_motif_tf_rna_concordance.py`; bio env, seed=42, PROJ, cache→/tmp;
   seed は permutation/sampling 部のみ有効と明記）。
5. **図・本文作成前に STOP & go/no-go report**（実施 / 出力パス / 赤旗 / 未解決 / 推奨1つ）。
6. `RESOLUTION.md` に判断・逸脱を記録。
7. その後に descriptive trajectory plot（既存 donor-level 値の YA/HA/AD 群推定; 新規 inference なし）
   および原稿作業。
- pre-registration commit に **annotated tag**。事後に仕様変更が要る場合は**原文を上書きせず
  amendment**（日時・理由・影響範囲を追記）。

---

## OPEN（PI 判断待ち — commit/push 前に解決）
- **remote 未設定**（`git remote -v` 空）。step 2 の push 不可。
  → (a) remote を追加して push する／(b) ローカル commit + annotated tag のみで進める（timestamp は
  ローカルで確保; remote ができ次第 push）。
- commit ブランチ（推奨: 専用 `analysis/tf-rna-concordance`、既存トピックブランチ規約に合わせる）。
