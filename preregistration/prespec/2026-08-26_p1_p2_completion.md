# 事前登録 — P1 / P2 の完了

作成: 2026-08-26。**原稿を編集する前**に commit（CLAUDE.md 鉄則5）。
BLOCKER と P0 は既に完了（`audit_log/2026-08-26_manuscript_finalisation/`）。本 run は残る
**P1 20 項目**と **P2**（数値訂正・図・引用・Declarations）を対象とする。

## 0. 重要な前提の確認（着手前に確定させた）

**P1-16 と P1-19 は再計算を要しない。** いずれも監査が既に計算済みで、出力が実在する:
- `P1-19`: `results/de/atac_linkedpeak_enrichment.csv` の `pass` 列が **10/10 すべて False**。
  原因は全行で `n_fdr_concordant = 0`（ピーク単位で FDR を通るものが皆無）。
  すなわち事前指定 GO 基準（`obs>0 かつ p<0.05 かつ nfdr>=3`）は**この標本サイズでは構造上到達不能**。
- `P1-16`: `$AUDIT/stats/pass5_composition_compositional.csv`。**microgliosis は手法依存**
  （割合 MWU padj 0.0021 / **CLR padj 0.484** / ALR vs mOli 0.036 / ALR vs Astro 0.156 /
  ALR vs mGC 0.963 / beta-binom 0.0030）。NSC は全手法で生存（CLR 0.0011、beta-binom 0.0030）。

したがって**本 run に新規の重い計算は無い**。すべて既存出力の開示と表記の訂正である。
計算が必要と判明したら STOP して報告する（鉄則4）。

## 1. 判定ゲート（各項目に適用）

- **G-P1 出典**: 書く数値はすべて `(file, row, column)` から引き直す。
  **監査レポートの本文からは引かない**（実 CSV を正とする）。食い違いは記録する。
- **G-P2 二重言語**: EN を直したら JA・図凡例・`cover_letter.md`・`submission_checklist.md` を grep。
- **G-P3 主張の強度**: 開示を足した結果、既存の主張が過大になっていないか確認する。
  過大なら同時に弱める（開示だけ足して主張を残さない）。
- **G-P4 監査の誤りを再生産しない**: ピアが記録した監査側の誤り5件
  （`F-8-001` `F-2-001` `F-6b-010` `B-5` の機構 `F-6b-001` の集合）は**修正対象にしない**。
- **G-P5 最終確認**: 完了後に PDF を再生成し、**ページ画像を目視**する
  （2026-08-26 の PDF 化で、テキスト抽出では見えない誤りが2件見つかったため必須とする）。

## 2. 対象と扱い

### 2-1. 原稿に反映する P1（14 項目）

`P1-1` 原著 supplementary の共変量（性別・APOE・BRAAK・LATE・海馬硬化）が公開されている事実の開示
／`P1-4` 実効ドナー数（Micro・Endo は n=7/9）／`P1-5` GSE325391・GSE264692 の原著引用と既報性
／`P1-7` "for the first time" の優先権主張の撤回／`P1-8` protein-level 3 参照の実体の正確な記述
／`P1-9` QC・doublet・ambient・統合の Methods 追記／`P1-10` 死後脳の交絡を Limitations に
／`P1-11` GSE199243 v1/v2 と Leiden パラメータ・seed／`P1-12` 未使用データセット 2 件の理由
／`P1-14` 多重補正のスコープ／**`P1-16` 組成解析の手法依存性**／`P1-17` ATAC の QC 指標
／`P1-18` 事前登録の DEVIATIONS 節／**`P1-19` A2 の GO 基準 0/10**／`P1-20` データセット選択基準。

### 2-2. 変更案の提示にとどめる P1（2 項目）

`P1-3`（原著著者への照会）と `P1-13`（Zenodo deposit）は**外部作業**。
原稿には既に開示済み（Limitation 8／Data availability）なので、追加は行わない。

### 2-3. P2

`FINAL_REPORT.md` §7 の P2 節に列挙された数値訂正・図の修正・引用の修正・Declarations を、
**1件ずつ出典に当たって**処理する。`$AUDIT/citation_check.tsv` の非 OK 29 件も同様に確認する。

## 3. 判定の三値

各項目を **APPLIED**（原稿を直した）/ **ALREADY-OK**（確認したが変更不要）/
**DEFERRED**（外部作業が要る。理由と必要な作業を明記）に分類し、RESOLUTION に全件記録する。
**「確認していないもの」を ALREADY-OK に丸めない。**

## 4. やらないこと

- 新規の重い計算（要ると判明したら STOP）
- `figures/` の図そのものの再生成（ただし本文と食い違う凡例テキストは直す）
- `raw/` `processed/` への書き込み
- 監査側の誤り5件の「修正」

## 5. 出力先

`audit_log/2026-08-26_p1_p2_completion/RESOLUTION.md`
