#!/bin/bash
# pull variants from a region across multiple VCFs

# change these to your gene of interest
CHR="chr2"
START=166845670
END=166984524
OUT="region_variants"

mkdir -p $OUT

# add your VCF paths here
VCFS=(
  "gs://fc-secure-path/to/proband.vcf.gz"
  "gs://fc-secure-path/to/mother.vcf.gz"
  "gs://fc-secure-path/to/father.vcf.gz"
)

# pull the region from each vcf
for vcf in "${VCFS[@]}"; do
    name=$(basename "$vcf" | cut -d'.' -f1)
    echo "extracting $name..."
    
    gsutil -u terra-ababc052 cat "$vcf" 2>/dev/null | \
      zcat 2>/dev/null | \
      awk -v c="$CHR" -v s="$START" -v e="$END" \
        '$1==c && $2>=s && $2<=e && !/^#/' > "$OUT/${name}.txt"
done

# quick summary
echo ""
echo "=== Summary ==="
echo "Region: $CHR:$START-$END"

total=0
for f in $OUT/*.txt; do
    n=$(wc -l < "$f" | tr -d ' ')
    total=$((total + n))
    echo "$(basename $f .txt): $n"
done

echo ""
echo "Total: $total"
echo "Unique: $(cat $OUT/*.txt | awk '{print $1,$2,$4,$5}' | sort -u | wc -l | tr -d ' ')"
