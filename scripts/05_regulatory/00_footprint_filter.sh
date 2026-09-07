#!/bin/bash
# FP step 1 (TRACKED, reproducible) — build the YA/HA niche-cell barcode map and filter the 81 GB
# aggregated ATAC fragment file to per-celltype Tn5 insertions (donor-tagged) for donor-level
# footprinting. Large outputs go to processed/footprint_tmp/ (gitignored); THIS pipeline logic is
# tracked so footprint_test.csv / panel D are regenerable from a fresh clone (given the raw deposit).
# Usage:  bash scripts/05_regulatory/00_footprint_filter.sh
set -euo pipefail
cd "$(dirname "$0")/../.."                      # -> project root
TMP=processed/footprint_tmp; mkdir -p "$TMP"
PY=/home/neurofuture/miniforge3/envs/bio/bin/python
FRAG=raw/GSE268609_lazarov2024/GSE268609_Aggregated_atac_fragments.tsv.gz

# (a) barcode -> ct,group,donor map (Astro/Micro/Endo x YA/HA) from the Seurat metadata.
#     Fragment barcode "SEQ-N" matches the metadata index exactly (N = orig.ident).
$PY - <<'PYEOF'
import pandas as pd
m = pd.read_csv("processed/per_dataset/GSE268609_metadata.csv", index_col=0)
CL = {"Astrocytes": "Astro", "Microglia": "Micro", "Endothelial": "Endo"}
m["ct"] = m["Cluster"].map(CL)
t = m[(m.ct.notna()) & (m.Group.isin(["YA", "HA"]))].reset_index()
t = t.rename(columns={"index": "barcode", "Group": "group", "SampleNumber": "donor"})
t[["barcode", "ct", "group", "donor"]].to_csv(
    "processed/footprint_tmp/barcode_map.tsv", sep="\t", index=False, header=False)
print("barcode map rows:", len(t))
PYEOF

# (b) one-pass filter: keep target-cell fragments, emit BOTH Tn5 cut sites (start,end).
#     NOTE: cut sites are UNSHIFTED (the 10x +4/-5 Tn5 offset is not applied). This is immaterial
#     to the YA-vs-HA *same-site* differential (any constant offset cancels), and likewise the Tn5
#     sequence bias cancels in the differential, so no bias model is used. Absolute footprint shape
#     is therefore only approximate; only the YA-vs-HA contrast is interpreted.
zcat "$FRAG" | awk -v OFS='\t' '
 FNR==NR { ct[$1]=$2; grp[$1]=$3; don[$1]=$4; next }
 /^#/ { next }
 ($4 in ct) {
   c=ct[$4]; cmd="gzip -c > processed/footprint_tmp/ins_" c ".tsv.gz"
   print $1, $2, grp[$4], don[$4] | cmd
   print $1, $3, grp[$4], don[$4] | cmd
 }
' processed/footprint_tmp/barcode_map.tsv -
echo "FILTER_DONE"
for c in Astro Micro Endo; do
  echo "  ins_$c: $(zcat processed/footprint_tmp/ins_$c.tsv.gz | wc -l) insertions"
done
