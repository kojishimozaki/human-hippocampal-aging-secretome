#!/usr/bin/env bash
# M1 methylome — fetch GENCODE v44 (hg38) gene annotation for strand-aware TSS / gene-body.
# The repo's CellRanger features.tsv lacks strand; minus-strand promoters (TSS+/-1kb) need it.
# Output: refs/gencode_v44_genes.tsv.gz (ensg, symbol, chrom, start, end, strand, gene_type), gene-level only.
# Build: hg38 / GRCh38, chr-prefixed (matches ALLC, BSgenome.Hsapiens.UCSC.hg38, GSE268609 peaks). seed-free.
set -euo pipefail
PROJ="${PROJ:-/home/neurofuture/Bioanalysis/01.SGZ_aging_project}"
cd "$PROJ/refs"
URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gtf.gz"
echo "fetching $URL"
curl -fsS --max-time 600 -o gencode.v44.basic.annotation.gtf.gz "$URL"
gunzip -t gencode.v44.basic.annotation.gtf.gz && echo "gzip OK"

# extract gene-level records -> compact TSV (version-stripped ENSG)
zcat gencode.v44.basic.annotation.gtf.gz \
| awk -F'\t' 'BEGIN{OFS="\t"; print "ensg","symbol","chrom","start","end","strand","gene_type"}
  $3=="gene"{
    g=""; s=""; t="";
    n=split($9,a,";");
    for(i=1;i<=n;i++){
      if(a[i]~/gene_id/){gsub(/.*gene_id "/,"",a[i]); gsub(/".*/,"",a[i]); sub(/\.[0-9]+$/,"",a[i]); g=a[i]}
      else if(a[i]~/gene_name/){gsub(/.*gene_name "/,"",a[i]); gsub(/".*/,"",a[i]); s=a[i]}
      else if(a[i]~/gene_type/){gsub(/.*gene_type "/,"",a[i]); gsub(/".*/,"",a[i]); t=a[i]}
    }
    print g,s,$1,$4,$5,$7,t
  }' | gzip > gencode_v44_genes.tsv.gz

echo "gene records: $(zcat gencode_v44_genes.tsv.gz | tail -n +2 | wc -l)"
rm -f gencode.v44.basic.annotation.gtf.gz   # keep only the compact table
echo "wrote refs/gencode_v44_genes.tsv.gz"