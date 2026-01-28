#!/bin/bash

# Change this region to fit your needs, this is only an example
REGION_CHR="chr2"
REGION_START=166845670
REGION_END=166984524
OUTPUT_DIR="scn1a_variants"
SUMMARY="$OUTPUT_DIR/summary.txt"

mkdir -p $OUTPUT_DIR

# GA4K workspace VCFs
YOUR_VCFS=(
  "gs://fc-secure-path/to/the/vcf"
  "gs://fc-secure-path/to/the/vcf"
)

# DCC workspace VCFs
DCC_VCFS=(
  "gs://fc-secure-path/to/the/vcf"
  "gs://fc-secure-path/to/the/vcf"
)

# Extract region from VCF
extract_region() {
    local vcf_path="$1"
    local output_file="$2"
    
    gsutil -u terra-ababc052 cat "$vcf_path" 2>/dev/null | \
      zcat 2>/dev/null | \
      awk -v chr="$REGION_CHR" -v start="$REGION_START" -v end="$REGION_END" \
        '$1==chr && $2>=start && $2<=end && !/^#/' > "$output_file"
}

# Process all VCFs
for vcf in "${YOUR_VCFS[@]}" "${DCC_VCFS[@]}"; do
    sample=$(basename "$vcf" | cut -d'.' -f1)
    extract_region "$vcf" "$OUTPUT_DIR/${sample}.txt"
done

# Generate summary
{
    echo "SCN1A Variant Extraction Summary"
    echo "Region: $REGION_CHR:$REGION_START-$REGION_END"
    echo "Date: $(date)"
    echo ""
    
    total=0
    for f in $OUTPUT_DIR/*.txt; do
        [ "$f" = "$SUMMARY" ] && continue
        count=$(wc -l < "$f")
        total=$((total + count))
        echo "$(basename $f .txt): $count variants"
    done
    
    echo ""
    echo "Total variants: $total"
    
    unique=$(cat $OUTPUT_DIR/*.txt 2>/dev/null | awk '{print $1,$2,$4,$5}' | sort -u | wc -l)
    echo "Unique positions: $unique"
    
    common=$(cat $OUTPUT_DIR/*.txt 2>/dev/null | awk '{print $1,$2,$4,$5}' | sort | uniq -c | awk '$1==4' | wc -l)
    echo "Common to all 4 samples: $common"
    
} > "$SUMMARY"

cat "$SUMMARY"
