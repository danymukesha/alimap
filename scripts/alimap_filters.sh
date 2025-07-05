#!/bin/bash
# alimap_enhanced.sh
# Enhanced BAM filtering and QC tool with blacklist exclusion
# Author: Dany Mukesha + enhancements 2025-07-05

set -euo pipefail

# Check dependencies
command -v samtools >/dev/null || { echo >&2 "Error: samtools not installed."; exit 1; }
command -v Rscript >/dev/null || { echo >&2 "Error: Rscript not installed."; exit 1; }

usage() {
    cat <<EOF
Usage: $0 <input_bam> [--blacklist <bed_file>] [--threads <num>]

Arguments:
  <input_bam>           Input BAM file (indexed)
  --blacklist <bed>     Optional: BED file of blacklist regions to exclude
  --threads <num>       Number of threads for samtools (default: 1)

Example:
  $0 input.bam --blacklist hg38_blacklist.bed --threads 4

Output:
  filtered_reads/          Filtered BAM files
  qc/                      QC reports and stats
  filter_counts.txt        Read counts after each filter
  alimap_filters_plot.pdf  Bar plot of read counts per filter
EOF
    exit 1
}

if [ $# -lt 1 ]; then
    usage
fi

input_bam=""
blacklist_bed=""
threads=1

# Parse args
while (( "$#" )); do
    case "$1" in
        --blacklist)
            blacklist_bed=$2
            shift 2
            ;;
        --threads)
            threads=$2
            shift 2
            ;;
        -*)
            echo "Unknown option: $1"
            usage
            ;;
        *)
            if [ -z "$input_bam" ]; then
                input_bam=$1
                shift
            else
                echo "Unexpected argument: $1"
                usage
            fi
            ;;
    esac
done

if [ ! -f "$input_bam" ]; then
    echo "Error: Input BAM file '$input_bam' does not exist."
    exit 1
fi

if [ ! -f "${input_bam}.bai" ]; then
    echo "Error: BAM index file '${input_bam}.bai' missing. Please index BAM."
    exit 1
fi

if [ -n "$blacklist_bed" ] && [ ! -f "$blacklist_bed" ]; then
    echo "Error: Blacklist BED file '$blacklist_bed' does not exist."
    exit 1
fi

output_dir="filtered_reads"
qc_dir="qc"
mkdir -p "$output_dir" "$qc_dir"

filter_counts_file="filter_counts.txt"
[ -e "$filter_counts_file" ] && rm "$filter_counts_file"

declare -A filters
filters["no_filter"]=""
filters["exclude_unmapped"]="-F 4"
filters["exclude_secondary"]="-F 256"
filters["mapq_30"]="-q 30"
filters["properly_paired"]="-f 2"

command -v bedtools >/dev/null || {
    echo "Warning: bedtools not found; blacklist exclusion will be skipped."
    blacklist_bed=""
}

filter_and_count_reads() {
    local filter_name=$1
    local filter_flag=$2
    local base_out="$output_dir/filtered_${filter_name}.bam"

    echo "Applying filter: $filter_name with flag '$filter_flag'..."

    if [ -n "$filter_flag" ]; then
        samtools view -@ "$threads" -b $filter_flag "$input_bam" -o "$base_out"
    else
        cp "$input_bam" "$base_out"
    fi

    if [ -n "$blacklist_bed" ]; then
        echo "Excluding blacklist regions from $base_out..."
        # bedtools subtract requires BAM -> BED conversion and back
        # To keep it efficient, use bedtools intersect with -v to exclude blacklist

        bedtools bamtobed -bedpe -i "$base_out" > "$output_dir/${filter_name}.bedpe"

        bedtools intersect -v -a "$output_dir/${filter_name}.bedpe" -b "$blacklist_bed" > "$output_dir/${filter_name}.filtered.bedpe"

        awk '{print $4}' "$output_dir/${filter_name}.filtered.bedpe" | sort | uniq > "$output_dir/${filter_name}.keep_reads.txt"

        samtools view -@ "$threads" -b -N "$output_dir/${filter_name}.keep_reads.txt" "$base_out" -o "${base_out%.bam}.no_blacklist.bam"

        mv "${base_out%.bam}.no_blacklist.bam" "$base_out"

        rm "$output_dir/${filter_name}.bedpe" "$output_dir/${filter_name}.filtered.bedpe" "$output_dir/${filter_name}.keep_reads.txt"
    fi

    samtools index "$base_out"

    count=$(samtools view -c "$base_out")
    echo "$filter_name: $count" >> "$filter_counts_file"
    echo "Filter '$filter_name' remaining reads: $count"
}

for name in "${!filters[@]}"; do
    filter_and_count_reads "$name" "${filters[$name]}"
done
echo "Generating QC reports..."

samtools flagstat "$input_bam" > "$qc_dir/input_flagstat.txt"
samtools idxstats "$input_bam" > "$qc_dir/input_idxstats.txt"

for bam in "$output_dir"/*.bam; do
    base=$(basename "$bam" .bam)
    samtools flagstat "$bam" > "$qc_dir/${base}_flagstat.txt"
    samtools idxstats "$bam" > "$qc_dir/${base}_idxstats.txt"
done

cat << 'EOF' > visualize_filters.R
library(ggplot2)
library(scales)

filter_data <- read.table("filter_counts.txt", sep=":", header=FALSE, col.names=c("Filter", "Reads"))
filter_data$Filter <- factor(filter_data$Filter, levels=filter_data$Filter)

p <- ggplot(filter_data, aes(x=Filter, y=Reads, fill=Filter)) +
  geom_bar(stat="identity", color="black") +
  geom_text(aes(label=scales::comma(Reads)), vjust=-0.3, size=3) +
  scale_y_continuous(labels = scales::comma) +
  labs(title="Reads Remaining After Filters",
       x="Filter Applied",
       y="Number of Reads",
       caption="alimap_enhanced.sh output") +
  theme_minimal() +
  theme(legend.position="none") +
  scale_fill_brewer(palette="Set2")

pdf("alimap_filters_plot.pdf", width=8, height=5)
print(p)
dev.off()
EOF

echo "Running R script to generate plot..."
Rscript visualize_filters.R

echo "Done! Outputs:"
echo "- Filtered BAMs in $output_dir/"
echo "- QC reports in $qc_dir/"
echo "- Filter counts in $filter_counts_file"
echo "- Plot: alimap_filters_plot.pdf"
