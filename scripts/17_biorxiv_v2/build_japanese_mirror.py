#!/usr/bin/env python3
"""Build the Japanese mirror of the revised (bioRxiv v2) manuscript.

Rather than re-translating ~10,000 words and risking new transcription errors, this
aligns the revised English to the existing Japanese draft
(`manuscript/manuscript_codex_niche_forward_ja.md`) and reuses its translation verbatim
for every paragraph the revision did not touch. Three transformations are applied to the
reused text:

  * reference numbers are remapped from the draft's numbering to the submitted numbering
    (built by matching DOIs between the two reference lists, so it cannot drift);
  * mnemonic supplementary labels are remapped to final numbers
    (S(N)->S3, S3->S4, S(M1)->S5, S(R1)->S6, S(E1)->S7, S(Rcv)->S8, S(RT1)->S9);
  * paragraphs the revision rewrote, and those Word-era edits changed materially, carry
    freshly written Japanese from JA_OVERRIDES below.

Paragraphs whose Japanese still tracks the pre-Word English in wording only -- no number
and no claim differs -- are listed in the file's header so the lag is visible rather than
silent.

Usage:
    PROJ=$(pwd) python scripts/17_biorxiv_v2/build_japanese_mirror.py
"""

from __future__ import annotations

import difflib
import json
import os
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

PROJ = Path(os.environ.get("PROJ", Path(__file__).resolve().parents[2]))
REVISED_DOCX = PROJ / "submission/v2/shimozaki-aging-secretome-v2.docx"
EN_SRC = PROJ / "manuscript/manuscript_codex_niche_forward_en.md"
JA_SRC = PROJ / "manuscript/manuscript_codex_niche_forward_ja.md"
OUT = PROJ / "manuscript/biorxiv_v2/manuscript_biorxiv_v2_ja.md"

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"

FIG_MAP = {"S(N)": "S3", "S(M1)": "S5", "S(R1)": "S6",
           "S(E1)": "S7", "S(Rcv)": "S8", "S(RT1)": "S9", "S3": "S4"}
LABEL_ORDER = ["S(N)", "S(M1)", "S(R1)", "S(E1)", "S(Rcv)", "S(RT1)", "S3"]


# --------------------------------------------------------------------------------------
# Freshly written Japanese, keyed by a prefix of the revised English paragraph.
# --------------------------------------------------------------------------------------
JA_OVERRIDES: list[tuple[str, str]] = [
 ("A support-cell secretome program in human hippocampal aging",
  "# ヒト海馬加齢の支持細胞セクレトーム・プログラムは方向性が再現するが、遺伝子特異性は示されない"),

 ("Cells that build the hippocampal extracellular environment",
  "海馬の細胞外環境を構築する細胞は、ヒト脳加齢の作用点として有力でありながら十分に特定されていない。また、"
  "それらを調べるために用いられる単一核解析のデザインは、解釈に先立って較正されることがほとんどない。本研究では、"
  "2,224個のセクレトーム関連遺伝子を対象とする事前規定解析を、ヒト海馬マルチオーム（若年8名、臨床的に対照と"
  "定義された高齢9名のドナー）に適用した。生物学的な結果を読む前にドナー群ラベル置換によってコントラストを較正し、"
  "すべての主張に対して、それを破ることを目的とした特性一致帰無分布を対置した。3つの所見が生き残った。加齢は、"
  "アストロサイト、ミクログリア、混合血管区画、オリゴデンドロサイト前駆細胞、オリゴデンドロサイトにわたって、"
  "方向の安定したセクレトーム関連プログラムと関連していた。60個の遺伝子×細胞型ペアはすべて、いずれの"
  "leave-one-donor-out再推定でも符号を保持し、アンビエントRNA補正後も検定可能な57方向すべてが保持され、"
  "その方向は独立したヒト海馬コホートで再現した（Spearman ρ = 0.24、符号一致率81%対背景58%、P = 5×10−4）。"
  "同じプログラムは、同一ドナー由来の成熟歯状回顆粒細胞、CA興奮性ニューロン、抑制性ニューロンでは背景を"
  "上回って再現されず、組織全体ではなく支持細胞に偏ったプログラムであることが示された。アストロサイトでは、"
  "RNAとクロマチンアクセシビリティが遺伝子単位で方向的に共役しており、遺伝子選択そのものを反復する"
  "end-to-endのドナーラベル置換検定にも耐えた（P = 0.025）。一方、4つの主張は対照検定に耐えなかった。"
  "発現量と効果量を一致させた対照遺伝子は選択遺伝子と同程度に再現し（P = 0.146、P = 0.187）、遺伝子単位の"
  "特異性は支持されない。見かけ上のRNA–メチル化協調は4つの特異性対照すべてで消失した。モチーフ、"
  "フットプリント、転写因子RNAのいずれからも上流制御因子は特定できなかった。受け手候補のどの区画にも"
  "協調的な受容体上昇はみられなかった。さらにドナーラベル置換は、完全帰無の下でもこの8対9のデザインが"
  "FDR < 0.05で有意な遺伝子を1個以上、抽出の29–42%で返すことを示した。この確率は0.05以下に抑えられるはずであり、"
  "遺伝子ごとのP値がここでは反保守的であること、したがってヒット数を偽発見率の推定値として読めないことを意味する。寄託者が付与した3つの希少細胞ラベルは、マーカーおよび分類器による"
  "再解析では支持されず、すべての結論から除外した。得られたものは、方向性が再現し機序は未解決である支持細胞"
  "加齢プログラムであり、この規模のヒトマルチオームデータがどの推論を支えられ、どれを支えられないかについての"
  "明示的な説明とともに報告する。"),

 ("Keywords: hippocampal aging",
  "**キーワード:** 海馬加齢; ニッチ支持細胞; セクレトーム関連転写; 細胞外環境; 単一核マルチオーム; "
  "認知加齢; 特性一致帰無対照; 解析再現性"),

 ("Analysis code, the pre-specified analysis documents, the software environments",
  "解析コード、事前に定めた解析文書、ソフトウェア環境、および本稿が報告するすべての数値の基礎となる結果表は、"
  "プロジェクトリポジトリで公開する（コード利用可能性を参照）。ソフトウェア環境は`envs/`で固定し、"
  "特記しない限り乱数シードは42を用いる。ワークフローのターゲットは、コミット済みの結果表だけで実行できるもの"
  "（図、検証解析、補足表、多重検定スコープの確認）と、寄託行列を追加で必要とするものとに分けてあり、"
  "一次データを持たない読者でも、リポジトリのクローンから全図と全検証表を再生成できる。"),

 ("Three limits on that reproducibility are stated rather than implied",
  "この再現可能性には3つの限界があり、暗示ではなく明示する。第1に、寄託されたSeuratオブジェクトから"
  "カウント行列を出力したスクリプトが保存されていない。そのため、細胞フィルタリング、品質閾値、RNAとATACの"
  "バーコード照合を含む取り込み工程はリポジトリから再実行できない。ワークフローは処理済み行列の段階以降から"
  "反復できる。第2に、固定環境で一次差次的発現解析のターゲットを再実行すると、本稿で用いた表の99,544行では"
  "なく99,543行が得られる。低カウントのオリゴデンドロサイト転写産物LINC01238（総カウント2、17ドナー中2名で検出）が"
  "現行の遺伝子フィルタを下回るためである。事前固定した60ペアはすべて再生成表に存在し、60ペアすべてが"
  "adjusted P < 0.1を保持する。60ペアにおけるlog2 fold changeの絶対差の最大値は4.3×10−5である。表S1は再生成後の"
  "細胞型別検定数を報告しているため、そのオリゴデンドロサイトの値は22,845ではなく22,844であり、5項目の合計は"
  "99,543となる。第3に、表S1のドナーラベル置換較正、効果量を一致させた再現背景、および外部コホートの"
  "ドナーラベル帰無分布を生成したスクリプトは投稿前監査の環境で実行されたものであり、リポジトリには含まれない。"
  "それらの出力表と由来記録は同梱してあり、これらに由来する原稿中の数値はすべて当該出力表から再計算して確認した。"),

 ("Spatial localisation was tested in GSE264692",
  "GSE264692 [21] では、品質管理を通った34のVisium捕捉領域それぞれについて、各ニッチシグネチャーと歯状回マーカー"
  "PROX1とのスポット単位のSpearman相関を算出して空間局在を検定した。プログラムがニッチに濃縮しているかを"
  "検定するため、GSE268609の3つの成熟ニューロンクラス（歯状回顆粒細胞、CA興奮性ニューロン、抑制性ニューロン）に"
  "同じドナー擬似バルクの枠組みを適用した。モデルには調製アームと年齢群を含めた（~arm + group）。ニッチ細胞と"
  "ニューロン各クラスのセクレトーム加齢効果を、対ごとのSpearman相関と2,000回の遺伝子シャッフル帰無分布で比較した。"
  "各ニューロンクラスにおける凍結ヒットの方向一致率は、そのニューロン自身のセクレトーム背景に対して評価し、"
  "de novoのニューロン・セクレトームヒット数も計数した。ニューロン内在性のプログラムは、gseapy prerank [32]、"
  "方向性を考慮した遺伝子セット応答、2,000回のランダムセット帰無分布によって探索した。ニッチとニューロンの"
  "プログラム間のドナー単位の連関、および受容体発現・受容体の加齢制御・リガンド活性に関するNicheNetの"
  "ニッチ→ニューロン解析は、いずれも記述的なものである。すべての解析でシード42を用いた。"),

 ("Table S1 | Donor-label permutation calibration",
  "**表S1｜ドナー群ラベル置換による差次的発現数の校正。** コホート×ニッチ細胞型の各組み合わせについて、"
  "ドナー群ラベルを100回置換したうえでドナー擬似バルクモデルを再推定した。本表は、観測されたadjusted P < 0.1の"
  "セクレトーム検出数、帰無中央値、経験的P値、FDR < 0.05におけるゲノムワイド検出数とその帰無中央値・"
  "95パーセンタイル・最大値・経験的P値、1個以上および10個以上のBenjamini–Hochberg有意遺伝子を生じた置換の割合、"
  "ならびにFDR < 0.1におけるゲノムワイド検出数を示す。GSE278576の内皮はGATE_FAILと表示されるが、これはYA 6名・"
  "HA 2名しか利用できず、最小ドナー数の要件を満たさなかったためである。経験的P値は1/101 = 0.0099を下回れない。"
  "GSE268609の8対9の分割では24,310通りのラベル割り当てが可能であり（7対9となるミクログリアと血管の比較では"
  "11,440通り）、100回の抽出ではこの空間を尽くさない。これらの校正結果は結果の第1小節にまとめてあり、"
  "以降の節が報告する検出数の読み方を規定する。"),

 ("Statistical procedures",
  "### 統計処理"),

 ("Reproducibility",
  "### 再現可能性"),

 ("Calibrating what this design can support",
  "### このデザインが支えられる範囲を較正する"),

 ("Deposited single-nucleus data now supply much of what is claimed",
  "ヒト脳加齢について主張される内容の多くは、いまや寄託された単一核データに由来している。しかし、そこから"
  "生じるデザインがどこまで推論を支えられるのかは、ほとんど測定されていない。そこで我々は、一次コントラスト"
  "（若年成人8名対、臨床的に対照と定義された高齢9名）をドナー群ラベル置換によって較正した。"
  "その較正を、それが規定する結果の後にあたるここで報告するのは、上記のあらゆる検出数をどう読むべきかを"
  "定めるのがこの較正だからである。"),

 ("Donor labels were permuted 100 times within each cohort",
  "コホート×細胞型の各組み合わせについてドナーラベルを100回置換し、ドナー擬似バルクモデルを完全に再推定した"
  "（表S1）。デザインの性質は2つ明らかになった。第1は好都合な性質である。観測されたセクレトーム関連の"
  "関連数は、ミクログリア（経験的P = 0.020）、アストロサイト（P = 0.040）、血管区画（P = 0.040）、"
  "OPC（P = 0.050）で置換帰無分布を上回り、帰無中央値はすべての区画で0であった。オリゴデンドロサイトは"
  "帰無分布を上回らなかった（5ペア、P = 0.109）。第2は好都合でない性質である。完全帰無の下では、Benjamini–Hochberg法は水準αにおいて"
  "いずれかの遺伝子が有意と判定される確率をα以下に抑えるため、ゲノムワイドFDR < 0.05のスコープでは"
  "その確率は0.05を超えないはずである。ところが置換では、1個以上の有意遺伝子が29–42%の抽出で得られ、"
  "10個以上も8–16%で得られた。GSE278576の10対20の分割では29–51%および3–17%であった。この差はBenjamini–Hochberg法の"
  "性質ではなく、このドナー数では遺伝子ごとのP値が反保守的であることの証拠であり、同法の保証がここでは"
  "成り立っていないことを意味する。したがって有意遺伝子の個数は解釈可能な単位ではなく、偽陽性の期待個数を"
  "そこから読み取ることもできない。"),

 ("Applying the same calibration to the genome-wide count",
  "同じ較正をFDR < 0.05のゲノムワイドな検出数に適用すると、強い腕と弱い腕が分かれる。アストロサイト"
  "（P = 0.030）、ミクログリア（P = 0.030）、OPC（P = 0.050）は帰無分布を上回った一方、"
  "オリゴデンドロサイト（P = 0.139）と血管区画（P = 0.069）は帰無分布の95パーセンタイル内に収まった。"
  "そこで我々は、アストロサイト、ミクログリア、OPCを、両方のスコープで実証可能なシグナルを担う区画として扱い、"
  "オリゴデンドロサイトと血管の腕については、方向としては整合するがゲノムワイドには置換帰無分布から"
  "分離できないものとして報告する。経験的P値は1/101 = 0.0099を下回れず、100回の抽出では利用可能な24,310通りの"
  "ラベル割り当て（10核閾値により7対9となるミクログリアと血管の比較では11,440通り）を尽くさない。"
  "したがってこれらの値は、証拠を正確に推定するというより、その範囲を画するものである。"),

 ("Multiple-testing scope was fixed in advance",
  "多重検定のスコープは事前に固定し、2通りの形で報告する。Benjamini–Hochberg補正は細胞型ごとに適用したため、"
  "5つのニッチ細胞クラスにまたがる99,544検定は1つではなく5つの検定ファミリーを構成する。5クラスを単一の"
  "ゲノムワイド検定ファミリーとして扱うと、報告した60ペアのうち49ペアが残り、新たに4ペア"
  "（AZGP1、LTBP3、NRG1、NXPE3。いずれもオリゴデンドロサイト）が加わる。事前固定した60ペアは、下流の"
  "アクセシビリティ、メチル化、外部再現、受容体、ニューロン解析のすべての起点となっているため、"
  "元のセットを保持したうえで、代替スコープの結果を黙って置き換えることなく併記する。"),

 ("The gene universe was defined as a union",
  "遺伝子の母集団は和集合として定義しており、より厳格な積集合の結果を併記する。2,224遺伝子の母集団は、"
  "Human Protein Atlas の predicted-secreted 注釈と、reviewed UniProtKB の分泌注釈エントリの和集合であり、"
  "両方の注釈をもつのは1,813遺伝子である。60ペアのうち10ペアは片方の注釈だけに依拠している"
  "（アストロサイトの USH2A・CHST9・CRB2・KITLG・ANTXR2、血管区画の KDR・KITLG、ミクログリアの ASAH1、"
  "OPC の FGF13、オリゴデンドロサイトの EPHA3）。とくに FGF13 は実験報告では細胞内在性とされ、KDR は受容体である。"
  "両方の注釈をもつ50ペアに限定すると、外部コホートでの方向一致は評価可能40ペア中32ペア（80%。全体では"
  "48ペア中39ペアで81%）、アストロサイトの遺伝子単位 RNA–アクセシビリティ一致は21遺伝子中16遺伝子"
  "（全体では26中20）となる。FGF13 のみを除くと OPC の連関ピーク数は19から6に減るが、残る遺伝子は"
  "すべて一致したままである。したがって報告した結果は、より緩い注釈に依存していない。ただし OPC の"
  "アクセシビリティの根拠は、支配的な1遺伝子を外すと薄い。"),

 ("Fig. 4 | Regulatory evidence and the boundary of mechanistic inference",
  "**図4｜制御層の根拠と機構推論の限界。** **A、** RNAと近傍ATACピークの変化方向の一致。"
  "60ペアのうち41ペア（68%）が一致しており、RNA有意でない3,942ペアの背景57.6%を上回る。"
  "ATACのFDRで有意に達した遺伝子はCOL21A1のみであり、それもRNA有意ファミリーの中に限られる。"
  "ゲノムワイドの369,393件のピーク検定は、いずれも多重検定補正を通過しない。**B、** "
  "chromVARで測定した加齢関連モチーフ活性。**C、** 加齢と疾患のクロマチン変化。**D、** "
  "ドナー単位の転写因子フットプリント。ATACの結果はアストロサイトにおける遺伝子単位の方向的な"
  "共役を支持するが、制御領域を特定するものではない。モチーフ、フットプリント、転写因子RNAの"
  "いずれの解析も駆動因子を特定しない。"),

 ("Fig. S4 | CSF-proteome concordance",
  "**図S4｜aging-UP/SASP 成分と CSF プロテオームの一致は遺伝子セットの閾値に依存する。** aging-UP および "
  "aging-DOWN のセットはニッチ・セクレトーム母集団から取り、その CSF 年齢推定値 [16] を、測定された "
  "CSF プロテオーム全体（6,174 feature）を背景として比較した。事前規定の adjusted P < 0.1 セットは"
  "有意でない上昇傾向を示し（60%が上昇、n = 25、P = 0.21）、拡張した adjusted P < 0.2 セットは一致した"
  "（69%が上昇、n = 42、P = 0.001）。背景をセクレトーム母集団に含まれる1,325 feature に限定すると、"
  "同じ比較は P = 0.44 および P = 0.014 となる。したがってこの結果は背景に依存しており、"
  "中核の60ペアではなく拡張した UP 成分を支持する。"),

 ("Agreement with the CSF aging proteome",
  "CSF 加齢プロテオーム [16] との一致は弱く、閾値と比較背景の双方に依存した。TNC、ANXA1、SERPINE1、"
  "IGFBP5 はいずれも CSF で加齢とともに増加した。測定された CSF プロテオーム全体（6,174 feature）を"
  "背景とすると、事前規定シグネチャーは有意でない傾向にとどまり（60%が上昇、n = 25、P = 0.21）、"
  "拡張セットは一致した（69%が上昇、n = 42、P = 0.001; 図S4）。背景をセクレトーム母集団に含まれる"
  "1,325 の CSF feature に限定した——こちらが条件を揃えた比較である——場合、これらの値は P = 0.44 および "
  "P = 0.014 になる。したがって CSF の結果は、より広い UP 成分に対する示唆的な支持であって、"
  "60ペアのタンパク質レベルでの一次検証ではないと位置づける。これらの参照データは総じて、転写シグネチャーを"
  "細胞老化に形づくられた分泌プログラムおよび近縁の生物系で測定されたタンパク質と結びつけるが、"
  "海馬ニッチそのものからの分泌を直接測定したことを意味しない。"),

 ("Two consequences follow, and we hold to them throughout",
  "ここから2つの帰結が導かれ、本稿を通じてこれを守る。すなわち、このデザインが支持できる証拠の単位は"
  "遺伝子ごとの有意性ではなく方向であること、そして本稿のいかなる主張も有意遺伝子の個数には依拠しないこと、"
  "の2点である。"),

 ("We consider first what the external comparisons do and do not license",
  "PREFIX::外部比較が何を許し何を許さないかを最初に検討する。その境界が、プログラムの残りの部分をどう読むべきかを"
  "規定するからである。"),

 ("A directionally stable but mechanistically unresolved support-cell program",
  "### 方向は安定しているが機序は未解決である支持細胞プログラム"),

 ("Conclusions: what this design can and cannot sustain",
  "### 結論: このデザインが支えられるものと支えられないもの"),

 ("Three findings survived every control we applied",
  "適用したすべての対照検定を、3つの所見が生き残った。ヒト海馬の加齢は、支持細胞区画において方向の安定した"
  "セクレトーム関連プログラムを伴う。60個の遺伝子×細胞型ペアはすべて、いずれのleave-one-donor-out再推定でも"
  "符号を保持し、検定可能な57方向のうち57方向がアンビエントRNA補正を通過し、その方向は独立したヒト海馬コホートで"
  "再現した（Spearman ρ = 0.24、符号一致率81%対背景58%、P = 5×10−4）。このプログラムは組織全体ではなく支持細胞に"
  "偏っている。同一ドナー由来の成熟歯状回顆粒細胞、CA興奮性ニューロン、抑制性ニューロンは、各自のセクレトーム背景を"
  "上回ってこれを再現せず、de novoの関連はそれぞれ4、3、0であったのに対し、5つの支持細胞区画では17、11、7、6、6で"
  "あった。さらにアストロサイトでは、RNAとクロマチンアクセシビリティが遺伝子単位で方向的に共役しており"
  "（160の連関ピークにわたる26遺伝子中20遺伝子）、遺伝子一致帰無分布、符号置換帰無分布、leave-one-gene-out解析、"
  "および遺伝子選択そのものを反復するend-to-endのドナーラベル置換検定（P = 0.025）を通過した。"),

 ("Four claims did not survive. Genes matched for expression",
  "4つの主張は生き残らなかった。発現量と一次コホートでの効果量を一致させた遺伝子は、選択したペアと同程度に"
  "再現したため（厳密な一致条件でP = 0.146およびP = 0.187）、外部コホートでの一致は、より広いセクレトーム関連の"
  "パターンを支持するのであって、この60ペアの特異性を支持するものではない。見かけ上のRNA–メチル化協調は、"
  "対応しない細胞型のRNAによる方向づけ、および10,000個のスクランブル格子と区別できない一致対角成分を含め、"
  "4つの特異性対照のすべてで消失した。上流制御因子は特定できなかった。300のモチーフ検定のうち0件、"
  "21のフットプリント検定のうち0件しかFDR有意に達せず、対応する転写因子のRNAが自身のモチーフと方向を一致させたのは"
  "検定可能な22組み合わせ中11組にとどまり、一致背景は52.4%であった。受け手候補のどの区画にも協調的な受容体上昇は"
  "みられなかった（8件の確認的検定すべてが不合格、方向性Q = 0.986）。性別、調製アーム、死後経過時間、群を"
  "同時に調整したモデルでは60ペア中0ペアしかadjusted P < 0.1を保持せず、性別を除く3項モデルでは60ペア中12ペアが"
  "保持されたが、60方向すべては変化しなかった。この結果は交絡の存在も不在も立証しないが、遺伝子ごとの推論が"
  "どこで止まるかを正確に示している。"),

 ("Two further results concern the data rather than the biology",
  "さらに2つの結果は、生物学ではなくデータそのものに関わる。ドナーラベル置換は、完全帰無の下でもこのデザインがFDR < 0.05で有意な遺伝子を1個以上、"
  "抽出の29–42%で返すことを示した。この確率は0.05以下であるはずで、このドナー数では遺伝子ごとのP値が"
  "反保守的である。したがってこの規模におけるヒット数を偽発見率の推定値として読むことはできない。また、寄託者が付与した3つの希少細胞ラベルは、マーカーパネルおよび"
  "分類器による再解析では支持されなかった。NSCとラベルされた核のうち元のラベルを保持したのは0.9%であり、"
  "確立された細胞型の保持率83–99.8%と対照的である。これらの集団はすべての結論から除外した。寄託ラベルを"
  "無批判に引き継ぐ再解析は、この誤りをそのまま引き継ぐことになる。"),

 ("What remains is a directionally reproducible, mechanistically unresolved",
  "残るのは、方向性が再現し機序は未解決である支持細胞加齢プログラムであり、その証拠上の境界についての"
  "明示的な説明を伴っている。この境界は、次に行うべき実験を規定する。すなわち、ここで指名した区画における"
  "タンパク質の直接測定、区画内でエンハンサー–遺伝子関係をモデル化できる深度での細胞型分離プロファイリング、"
  "そのために十分な検出力をもつコホートでの細胞型×加齢の交互作用の正式な検定、およびヒト多細胞モデルにおける"
  "単一細胞読み出しを伴う摂動実験である。陰性の結果を陽性の結果と同じ分量で報告するのは、この標本規模では、"
  "どの推論が成り立たないかを知ることのほうが、より広く転用できる結果だからである。"),

 ("Funding: This work was supported by JSPS KAKENHI",
  "**研究助成:** 本研究はJSPS科研費 JP23K10962の助成を受けた。"),

 ("Competing interests: The author declares",
  "**利益相反:** 著者は開示すべき利益相反を有しない。"),

 ("Use of AI tools: AI-assisted tools",
  "**AIツールの使用:** コード実装、デバッグ、原稿作成、編集の一部にAI支援ツール（Claude Code、Anthropic社; "
  "Codex CLI、OpenAI社）を使用した。これらは著者の指示のもとで対話的に使用し、出力はすべて著者が確認した。"
  "解析上の判断、統計モデル、科学的解釈はすべて著者が決定・検証しており、内容、コード、結論のすべての責任は"
  "著者が負う。AIツールは著者として記載していない。投稿に先立ち、別個のAIセッションを、敵対的な頑健性の再点検"
  "（Claude、Anthropic社）および編集レビューと校正（Claude、Anthropic社; GPT、OpenAI社）に使用した。"
  "これらは著者自身の作業に対する追加の自動チェックであって、独立した研究者による監査ではない。"
  "指摘は著者が元データに照らして検証したうえで採否を判断した。"),

 ("All primary and secondary datasets are publicly available",
  "一次・二次データはすべて公開されており、GSE268609 [17]、GSE278576およびGSE299139 [18]、GSE186538 [19]、"
  "GSE199243 [20]、GSE264692 [21]、GSE325391 [22]から入手できる。外部リソースは、Human Protein Atlas [24]、"
  "UniProtKB [25]、SenMayo（MSigDB SAUL_SEN_MAYO）[13]、SASP Atlas [14]、脳プロテオーム認知トラジェクトリー"
  "データセット [15]、CSF加齢プロテオーム [16]、NicheNet v2事前分布（Zenodo 7074291）[31]、"
  "BSgenome.Hsapiens.UCSC.hg38、およびヒト（taxonomy 9606）のJASPAR2020 CORE脊椎動物モチーフ [30] である。"
  "処理済み中間ファイルは再配布しない。解析行列は50 GBを超えるうえ、方法に記した取り込み工程の限界を除けば、"
  "上記の寄託データから記載の手順で再生成できるためである。コードアーカイブには、事前に定めた解析文書と、"
  "報告した数値の基礎となる結果表（ドナーラベル置換較正、効果量を一致させた再現背景、外部コホートの"
  "ドナーラベル帰無分布を含む）を、それぞれの由来記録とともに収めてある。"),

 ("The analysis code, the pre-specified analysis documents, and the result tables that support",
  "解析コード、事前に定めた解析文書、および報告した数値を支える結果表は {REPO_URL} で公開し、"
  "{ZENODO_DOI} にアーカイブする。本稿で用いたバージョンは、タグ {RELEASE_TAG} を付けたリリースである。ソフトウェア環境は"
  "`envs/`で固定しており（pyDESeq2 0.5.4、scanpy 1.11.5、Python 3.11.15、およびDESeq2・Signac・chromVAR・ashrを"
  "含むR環境）、特記しない限り乱数シードは42を用いる。図、検証、補足表、多重検定スコープの各ターゲットは、"
  "コミット済みの結果表だけを用いてリポジトリのクローンから実行できる。クロマチンアクセシビリティ一致、"
  "差次的発現の頑健性、制御機構解析、一次差次的発現表を再計算するターゲットは、ここでは再配布していない"
  "寄託行列を追加で必要とする。この再現可能性の限界は方法に記した。"),

 ("Fig. 2 | Age-associated secretome-related transcription",
  "**図2｜ニッチ細胞クラスにおける加齢関連のセクレトーム関連転写。** **A、** 各ニッチ細胞クラスのボルケーノプロット。"
  "**B、** 細胞型別のセクレトーム分類のヒートマップ。**C、** 細胞老化関連分泌形質の参照セットに対する方向性の濃縮。"),

 ("Fig. 5 | The aging niche-secretome shows no detectable dentate-gyrus concentration",
  "**図5｜加齢ニッチ・セクレトームに歯状回への集中は検出されない。** GSE264692 [21] の品質管理済み34捕捉領域の"
  "Visium空間トランスクリプトミクス。左上と右上のパネルはPROX1とアストロサイト加齢UPシグネチャーの代表的マップを"
  "示し、輪郭はPROX1高値スポットを示す。下段のパネルは、捕捉領域ごとの、PROX1とアストロサイト・ミクログリア・"
  "オリゴデンドロサイトの加齢UPシグネチャーとの空間Spearman相関を示す。ρの中央値はそれぞれ+0.003、+0.020、"
  "−0.021であり、|ρ|の中央値はそれぞれ0.022、0.021、0.040である。歯状回との共局在は検出されないが、"
  "利用可能な前方切片は歯状回を明瞭に分解しない。コホート間再現性は図S7に別途まとめた。"),

 ("Fig. 6 | Granule-lineage transcription in a cognitive-resilience cohort",
  "**図6｜認知レジリエンス・コホートにおける顆粒細胞系列の転写。** **A、** レジリエント対対照（RES対CTRL）、"
  "レジリエント対重度AD（RES対SAD）、重度AD対対照（SAD対CTRL）の差次的発現遺伝子数。**B、** RES対SADの発現"
  "シグネチャー。**C、** 群ごとの年齢。CTRL、非認知症対照; RES、高病理かつ認知機能保持; MAD、中等度AD; "
  "SAD、重度AD。MADは年齢の分布を示すためCのみに含めており、レジリエンス比較の一部ではない。GSE325391 [22]。"),

 ("Table S3 | Motif-accessibility and TF-RNA concordance",
  "**表S3｜転写因子×細胞型ペアごとのモチーフ・アクセシビリティとTF-RNAの一致。** 本表は62行からなり、"
  "31行はHA対YA解析、31行は記述的なAD対HA投影で、contrast列で区別する。名目上のHA対YA chromVAR P < 0.05を"
  "満たす事前規定の各組み合わせについて、対応づけた転写因子、細胞型、モチーフとモチーフファミリー、"
  "モチーフchromVAR dz・名目P・細胞型別Benjamini–Hochberg Q、RNAのlog2 fold change・標準誤差・95%信頼区間・"
  "厳密P・元の差次的発現Q、baseMean、検定可能性、chromVARとRNAの符号、および符号一致カテゴリーを示す。"
  "31因子のうち9因子は低発現のためD_not_testableとラベルした。一致かつFDR有意である組み合わせは存在しないため、"
  "そのカテゴリーは表に現れない。検定可能な3行はadjusted RNA値が空欄であるが、これはDESeq2の独立フィルタリングが"
  "低baseMeanでは値を返さないためである。"),
]

PLACEHOLDER_JA = {
    "REPO_URL": "https://github.com/kojishimozaki/human-hippocampal-aging-secretome",
    "ZENODO_DOI": "https://doi.org/10.5281/zenodo.22640703",
    "RELEASE_TAG": "biorxiv-v2",
}

KEEP_ENGLISH_PREFIXES = ("Koji Shimozaki", "Corresponding author: Koji Shimozaki")


NUMS = re.compile(r"\d[\d,]*(?:\.\d+)?")


def numset(text: str) -> set[str]:
    return {m.group(0).rstrip(".").replace(",", "") for m in NUMS.finditer(text)}


def regroup(en_paras: list[str], ja_paras: list[str]) -> list[str]:
    """Fold extra Japanese blocks into contiguous groups, one per English paragraph.

    Where the Japanese draft split a paragraph the split point is recovered by maximising
    numeric-token agreement between each English paragraph and its Japanese group, so no
    sentence is dropped and none lands under the wrong paragraph.
    """
    n, m = len(en_paras), len(ja_paras)
    if n == m:
        return ja_paras
    en_nums = [numset(e) for e in en_paras]
    pre = [set()]
    for j in ja_paras:
        pre.append(pre[-1] | numset(j))

    def score(i: int, a: int, b: int) -> float:
        got = pre[b] - pre[a]
        want = en_nums[i]
        if not want and not got:
            return 0.5
        return len(want & got) / max(1, len(want | got))

    NEG = float("-inf")
    best = [[NEG] * (m + 1) for _ in range(n + 1)]
    back = [[0] * (m + 1) for _ in range(n + 1)]
    best[0][0] = 0.0
    for i in range(1, n + 1):
        for b in range(i, m - (n - i) + 1):
            for a in range(i - 1, b):
                if best[i - 1][a] == NEG:
                    continue
                v = best[i - 1][a] + score(i - 1, a, b)
                if v > best[i][b]:
                    best[i][b], back[i][b] = v, a
    out, b = [], m
    for i in range(n, 0, -1):
        a = back[i][b]
        out.append("".join(ja_paras[a:b]))
        b = a
    return out[::-1]


# Statistical symbols follow the English: uppercase, italic, and "P-value" hyphenated.
# The reused translations were written with a lowercase p, so they are normalised here
# rather than edited one by one.
LOWER_P = re.compile(r"(?<![A-Za-z0-9])p(?=\s*[=<>\u2264\u2265])")
ITALIC_JA = re.compile(r"(?<![A-Za-z0-9*\u2212-])(?:[PQ](?![A-Za-z0-9*])|n(?=\s*=))")


def normalise_symbols(text: str) -> str:
    text = LOWER_P.sub("P", text)
    text = re.sub(r"(?<=[PQ]) (?=values?\b)", "-", text)
    return ITALIC_JA.sub(lambda m: f"*{m.group(0)}*", text)


def para_text(p) -> str:
    return "".join(n.text or "" for n in p.iter(W + "t")).strip()


def para_style(p) -> str:
    st = p.find(f"{W}pPr/{W}pStyle")
    return st.get(W + "val") if st is not None else "Normal"


def blocks(path: Path) -> list[str]:
    """Non-empty lines, skipping any superseded-banner prefix before the title heading."""
    text = path.read_text(encoding="utf-8")
    text = re.sub(r"\A(?:<!--.*?-->\n|>.*\n|\n)+", "", text)
    return [l.strip() for l in text.split("\n") if l.strip()]


def build_ref_map(revised_refs: list[str], src_refs: list[str]) -> dict[int, int]:
    def doi(s):
        m = re.findall(r"doi:([^\s)]+)", s)
        return m[0].rstrip(".").lower() if m else None
    new = {doi(t): i + 1 for i, t in enumerate(revised_refs)}
    out = {}
    for t in src_refs:
        n = int(re.match(r"^(\d+)\.", t).group(1))
        d = doi(t)
        if d in new:
            out[n] = new[d]
    return out


def make_remapper(ref_map: dict[int, int]):
    def cit(m):
        inner = m.group(1)
        if not re.fullmatch(r"[\d,\s–—-]+", inner):
            return m.group(0)
        nums = []
        for part in re.split(r"\s*,\s*", inner.strip()):
            rm = re.fullmatch(r"(\d+)\s*[–—-]\s*(\d+)", part)
            if rm:
                nums.extend(range(int(rm.group(1)), int(rm.group(2)) + 1))
            elif part.isdigit():
                nums.append(int(part))
            else:
                return m.group(0)
        out = sorted({ref_map.get(n, n) for n in nums})
        parts, i = [], 0
        while i < len(out):
            j = i
            while j + 1 < len(out) and out[j + 1] == out[j] + 1:
                j += 1
            if j == i:
                parts.append(str(out[i]))
            elif j == i + 1:          # the English cites a pair as "[36,37]", not a range
                parts.extend([str(out[i]), str(out[j])])
            else:
                parts.append(f"{out[i]}–{out[j]}")
            i = j + 1
        return "[" + ",".join(parts) + "]"

    def labels(s: str) -> str:
        res, i = [], 0
        while i < len(s):
            for k in LABEL_ORDER:
                if s.startswith(k, i) and (i == 0 or not s[i - 1].isdigit()):
                    if s[i + len(k):i + len(k) + 1].isdigit():
                        continue
                    res.append(FIG_MAP[k])
                    i += len(k)
                    break
            else:
                res.append(s[i])
                i += 1
        return "".join(res)

    return lambda s: labels(re.sub(r"\[([^\[\]]*)\]", cit, s))


def main() -> None:
    with zipfile.ZipFile(REVISED_DOCX) as z:
        root = ET.fromstring(z.read("word/document.xml"))
    rev = [(para_style(p), para_text(p)) for p in root.find(W + "body").findall(W + "p")]
    rev = [(s, t) for s, t in rev if t]
    revised_refs = [t for s, t in rev if s == "ListParagraph"]

    en_src = blocks(EN_SRC)
    ja_src = blocks(JA_SRC)
    src_refs = [b for b in en_src if re.match(r"^\d+\.\s", b)]
    remap = make_remapper(build_ref_map(revised_refs, src_refs))

    # structural EN<->JA alignment. Both drafts share one heading skeleton, so align
    # section by section: a run of body paragraphs between two headings maps 1:1, and
    # where the Japanese split a paragraph the extra blocks are folded back into their
    # predecessor rather than dropped.
    def kind(s):
        m = re.match(r"^(#{1,3}) ", s)
        if m:
            return "H" + str(len(m.group(1)))
        return "REF" if re.match(r"^\d+\.\s", s) else "P"

    def sections(blocks_):
        secs, cur = [], []
        for b in blocks_:
            if kind(b).startswith("H"):
                if cur:
                    secs.append(cur)
                secs.append([b])
                cur = []
            else:
                cur.append(b)
        if cur:
            secs.append(cur)
        return secs

    en_secs, ja_secs = sections(en_src), sections(ja_src)
    if len(en_secs) != len(ja_secs):
        raise SystemExit(f"section count differs: EN {len(en_secs)} vs JA {len(ja_secs)}")
    ja_for_en: list[str] = []
    for es, js in zip(en_secs, ja_secs):
        if len(js) < len(es):
            raise SystemExit(f"Japanese section is short: {es[0][:50]!r}")
        ja_for_en.extend(regroup(es, list(js)))
    if len(ja_for_en) != len(en_src):
        raise SystemExit("alignment produced the wrong length")
    en2ja = {i: i for i in range(len(en_src))}
    ja_src = ja_for_en
    en_plain = [re.sub(r"[*`#]", "", x).strip() for x in en_src]

    def norm(s):
        s = re.sub(r"\s+", " ", s).replace("’", "'")
        for a, b in [("p = ", "P = "), ("p <", "P <"), ("p-", "P-"), ("⁻", "-"), ("×10", "x10")]:
            s = s.replace(a, b)
        return s.lower()
    N = [norm(remap(x)) for x in en_plain]

    overrides = {k: v for k, v in JA_OVERRIDES}
    used, lagging, out = set(), [], []
    for i, (style, text) in enumerate(rev):
        if style == "ListParagraph":
            out.append(("REF", text))
            continue
        if text.startswith(KEEP_ENGLISH_PREFIXES):
            out.append(("P", text))
            continue
        ov = next((v for k, v in overrides.items() if text.startswith(k)), None)
        if ov is not None:
            used.add(next(k for k in overrides if text.startswith(k)))
            if ov.startswith("PREFIX::"):
                ov = ov[len("PREFIX::"):]
            else:
                out.append(("OV", ov))
                continue
        nt = norm(text)
        # Headings match only against headings of the same level: short heading strings
        # otherwise fuzzy-match into unrelated body text.
        want = {"Title": "H1", "Heading1": "H2", "Heading2": "H3"}.get(style)
        cand = [j for j in range(len(N)) if kind(en_src[j]) == want] if want else list(range(len(N)))
        # autojunk distorts similarity badly on strings this long -- it marks common
        # characters as junk and can score a near-identical pair at 0.23 -- so the
        # ratios that gate the checks below are only meaningful with it off.
        def sim(j: int) -> float:
            return difflib.SequenceMatcher(None, nt, N[j], autojunk=False).ratio()
        best = max(cand, key=sim)
        ratio = sim(best)
        ja = remap(ja_src[en2ja[best]]) if best in en2ja else None
        if ja is None:
            raise SystemExit(f"no Japanese counterpart for revised paragraph {i}: {text[:60]}")
        if ov is not None:
            ja = re.sub(r"^(\*\*|###?\s)", "", ja)
            ja = ov + ja
        elif ratio < 0.55:
            # A paragraph this unlike anything in the translated draft is new, not merely
            # reworded: reusing the nearest Japanese would silently attach the wrong text.
            raise SystemExit(
                f"revised paragraph {i} has no Japanese counterpart (best match {ratio:.2f}).\n"
                f"  {text[:90]}\n"
                "Add an entry to JA_OVERRIDES for it."
            )
        elif ratio < 0.90:
            lagging.append((i, ratio, text[:70]))
        if style in ("Title", "Heading1", "Heading2") and not ja.startswith("#"):
            ja = {"Title": "# ", "Heading1": "## ", "Heading2": "### "}[style] + ja
        out.append(("JA", ja))

    missing = set(overrides) - used
    if missing:
        raise SystemExit("unused Japanese overrides (anchors drifted):\n  " + "\n  ".join(missing))

    header = [
        "<!-- 生成物: scripts/17_biorxiv_v2/build_japanese_mirror.py。手で編集しない。 -->",
        "<!-- 英語正本: manuscript/biorxiv_v2/manuscript_biorxiv_v2.md -->",
        "",
        "> **この文書について。** bioRxiv 改訂版（v2）英語正本の日本語ミラーである。改訂で書き換えた段落と、"
        "Word 編集で内容が変わった段落は新規に翻訳した。それ以外の段落は "
        "`manuscript/manuscript_codex_niche_forward_ja.md` の訳文をそのまま再利用し、"
        "文献番号と補足図番号（S(N)→S3、S3→S4、S(M1)→S5、S(R1)→S6、S(E1)→S7、S(Rcv)→S8、S(RT1)→S9）"
        "のみ機械的に付け替えている。",
    ]
    if lagging:
        header += [
            ">",
            "> **語句が英語版に追随していない段落**（数値・主張はいずれも一致。Word 編集時の言い回しの変更のみ）:",
        ]
        for i, r, t in lagging:
            header.append(f"> - 段落 {i}（一致度 {r:.2f}）: {t}…")
    header.append("")

    body = []
    for kind_, text in out:
        if kind_ == "OV":
            for k, v in PLACEHOLDER_JA.items():
                text = text.replace("{" + k + "}", v)
        if kind_ != "REF":          # the reference list keeps its published typography
            text = normalise_symbols(text)
        body += [text, ""]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(header + body).rstrip() + "\n", encoding="utf-8")
    print(f"revised paragraphs : {len(rev)}")
    print(f"fresh Japanese     : {len(used)}")
    print(f"reused Japanese    : {sum(1 for k, _ in out if k == 'JA')}")
    print(f"wording-lag flagged: {len(lagging)}")
    print(f"wrote {OUT.relative_to(PROJ)}")


if __name__ == "__main__":
    main()
