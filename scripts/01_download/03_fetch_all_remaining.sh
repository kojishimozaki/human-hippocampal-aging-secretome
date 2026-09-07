#!/usr/bin/env bash
# Download ALL remaining datasets + references (idempotent; wget -c). Tier-2.
set -uo pipefail
P=${PROJ:-/home/neurofuture/Bioanalysis/01.SGZ_aging_project}; RAW="$P/raw"; REFS="$P/refs"
mkdir -p "$RAW" "$REFS"
echo "=== disk before ==="; df -h "$P" | tail -1
prefix(){ echo "$1" | sed 's/[0-9]\{3\}$/nnn/'; }
fetch_geo(){
  local gse="$1" dest="$2" base files
  base="https://ftp.ncbi.nlm.nih.gov/geo/series/$(prefix "$gse")/$gse/suppl/"
  mkdir -p "$dest"; echo "[$(date +%H:%M:%S)] $gse"
  files=$(curl -s --max-time 90 "$base" | grep -oE 'href="[^"?]+"' | sed -e 's/href="//' -e 's/"//' | grep -vE '^(/|\.\.|https?:|\?)')
  [ -z "$files" ] && { echo "  !! no listing for $gse"; return; }
  for f in $files; do
    wget -q -c -P "$dest" "${base}${f}" && echo "  ok $f ($(du -h "$dest/$f" 2>/dev/null|cut -f1))" || echo "  FAIL $f"
  done
}
fetch_geo GSE185553 "$RAW/GSE185553_zhou_adult"
fetch_geo GSE160189 "$RAW/GSE160189_HPC_AP_axis"
fetch_geo GSE264624 "$RAW/GSE264624_tran2025_snrna"
fetch_geo GSE264692 "$RAW/GSE264692_tran2025_visium"
fetch_geo GSE248545 "$RAW/GSE248545_spatial_FISH"
fetch_geo GSE95753  "$RAW/GSE95753_mouse_DG"
fetch_geo GSE198323 "$RAW/GSE198323_zhou_AD"
# GSE268609 ATAC fragments (82GB); RDS/matrix handled by the Tier-1 job
B268="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE268nnn/GSE268609/suppl"
echo "[$(date +%H:%M:%S)] GSE268609 ATAC fragments (82GB)"
wget -q -c -P "$RAW/GSE268609_lazarov2024" "$B268/GSE268609_Aggregated_atac_fragments.tsv.gz" && echo "  ok ATAC" || echo "  FAIL ATAC"
# NicheNet v2 priors (Zenodo)
mkdir -p "$REFS/nichenet"
for u in \
  https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds \
  https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds \
  https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds ; do
  wget -q -c -P "$REFS/nichenet" "$u" && echo "  ok $(basename "$u")" || echo "  FAIL $(basename "$u")"
done
echo "[$(date +%H:%M:%S)] FETCH_ALL_REMAINING_DONE"
echo "=== disk after ==="; df -h "$P" | tail -1
du -sh "$RAW"/* 2>/dev/null
