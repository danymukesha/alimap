#!/usr/bin/env bash
# alimap_enhanced.sh
# Refactored, robust BAM filtering + QC tool with optional blacklist exclusion
# Author: Dany Mukesha (original) + refactor and enhancements (2025-08-11)
# Purpose: produce cumulative filters, robust logging, optional blacklist exclusion,
#          samtools + bedtools + R-based plotting, samtools-stats QC, and reproducible outputs.

set -euo pipefail
IFS=$'\n\t'

VERSION="2025-08-11"

print_help() {
    cat <<EOF
alimap_enhanced.sh v$VERSION

Usage: $(basename "$0") -i <input.bam> [-b <blacklist.bed>] [-t <threads>] [-o <outdir>] [--mode cumulative|independent] [--keep-temp]

Options:
  -i|--input        Input BAM file (required). Must be indexed (.bai or .csi). If missing index, script can create it.
  -b|--blacklist    Optional BED file of blacklist regions to exclude (requires bedtools >=2.25)
  -t|--threads      Number of threads for samtools (default: 1)
  -o|--outdir       Output directory (default: filtered_reads)
  -m|--mode         "cumulative" (apply filters sequentially) or "independent" (each filter applied to input). Default: cumulative
  --keep-temp       Keep temporary files (useful for debugging)
  -h|--help         Show this help and exit

Example:
  $(basename "$0") -i sample.bam -b hg38_blacklist.bed -t 4 -o results --mode cumulative

Outputs written to <outdir> and <outdir>/qc. Files:
  filter_counts.txt        Tab-delimited: filter\treads_remaining
  <outdir>/filtered_*.bam Filtered BAMs per filter step
  <outdir>/qc/*            QC outputs (samtools flagstat, idxstats, stats)
  alimap_filters_plot.pdf  Barplot of read counts per filter (requires R + ggplot2)

Notes:
  - If bedtools is not present, blacklist exclusion will be skipped and a warning will be shown.
  - This script aims for reproducible, audit-friendly filtering of BAM files.
EOF
}

# default params
threads=1
outdir="filtered_reads"
mode="cumulative"
keep_temp=0
blacklist_bed=""

# parsing args (getopt for long options)
if ! getopt --test >/dev/null; then
    echo "getopt not available; falling back to basic parsing." >&2
fi

ARGS=$(getopt -o i:b:t:o:m:kh --long input:,blacklist:,threads:,outdir:,mode:,keep-temp,help -n "$(basename "$0")" -- "$@") || { print_help; exit 1; }
eval set -- "$ARGS"
while true; do
    case "$1" in
        -i|--input) input_bam="$2"; shift 2;;
        -b|--blacklist) blacklist_bed="$2"; shift 2;;
        -t|--threads) threads="$2"; shift 2;;
        -o|--outdir) outdir="$2"; shift 2;;
        -m|--mode) mode="$2"; shift 2;;
        --keep-temp) keep_temp=1; shift;;
        -h|--help) print_help; exit 0;;
        --) shift; break;;
        *) echo "Internal error while parsing args"; exit 1;;
    esac
done

# Basic validation
if [ -z "${input_bam-}" ]; then
    echo "Error: input BAM required." >&2
    print_help
    exit 1
fi

if [ ! -f "$input_bam" ]; then
    echo "Error: Input BAM '$input_bam' not found." >&2
    exit 1
fi

# output dirs
mkdir -p "$outdir"
qc_dir="$outdir/qc"
mkdir -p "$qc_dir"

# temporary dir
tmpdir=$(mktemp -d -t alimap.XXXXXX)
trap '[[ $keep_temp -eq 1 ]] || rm -rf "$tmpdir"' EXIT

# checking the dependencies
missing_deps=()
for cmd in samtools awk sort uniq Rscript; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        missing_deps+=("$cmd")
    fi
done

# bedtools optional
bedtools_ok=1
if ! command -v bedtools >/dev/null 2>&1; then
    bedtools_ok=0
    echo "Warning: bedtools not found. Blacklist exclusion will be skipped if requested." >&2
fi

if [ ${#missing_deps[@]} -gt 0 ]; then
    echo "Error: missing required dependencies: ${missing_deps[*]}" >&2
    exit 1
fi

# controlling threads is integer
if ! [[ "$threads" =~ ^[0-9]+$ ]]; then
    echo "Error: threads must be an integer." >&2
    exit 1
fi

# controlling index; create if needed
index_file="${input_bam}.bai"
if [ ! -f "$index_file" ]; then
    echo "Index for $input_bam not found. Creating index with samtools index..."
    samtools index -@ "$threads" "$input_bam"
fi

# validating blacklist file if provided
if [ -n "$blacklist_bed" ]; then
    if [ ! -f "$blacklist_bed" ]; then
        echo "Error: blacklist BED '$blacklist_bed' not found." >&2
        exit 1
    fi
    if [ $bedtools_ok -eq 0 ]; then
        echo "Warning: bedtools not available; blacklist exclusion will be skipped." >&2
        blacklist_bed=""
    fi
fi

# Here, we define filters (ordered). Each item is: key|samtools_flags|description
# samtools_flags are applied using samtools view -b <flags>
# Use cumulative mode: each step starts from previous output. Independent mode: always from input_bam.
filters=(
  "no_filter||No filtering (copy input)"
  "exclude_unmapped|-F 4|Exclude unmapped reads"
  "exclude_secondary|-F 256|Exclude secondary alignments"
  "mapq_30|-q 30|Keep reads with MAPQ >= 30"
  "properly_paired|-f 2|Keep properly paired reads"
)

filter_counts_file="$outdir/filter_counts.txt"
: > "$filter_counts_file"

# a helper to run a filter step
run_filter() {
    local name="$1" flags="$2" desc="$3"
    local in_bam="$4"
    local out_bam="$5"

    echo "[INFO] Step: $name — $desc"

    if [ -z "$flags" ]; then
        # make a copy
        samtools view -@ "$threads" -b "$in_bam" -o "$out_bam"
    else
        # we use eval to preserve flags spacing
        eval samtools view -@ "$threads" -b $flags "\"$in_bam\"" -o "\"$out_bam\""
    fi

    # optional: exclude blacklist using bedtools intersect -abam <in> -b <blacklist> -v -bedpe? -ubam
    if [ -n "$blacklist_bed" ]; then
        echo "[INFO] Excluding blacklist regions using bedtools (this will replace $out_bam)"
        # bedtools intersect -abam IN -b blacklist -v > OUT
        tmp_no_black="$tmpdir/${name}.no_blacklist.bam"
        bedtools intersect -abam "$out_bam" -b "$blacklist_bed" -v > "$tmp_no_black"
        mv "$tmp_no_black" "$out_bam"
    fi

    # indexing the output
    samtools index -@ "$threads" "$out_bam"

    local cnt
    cnt=$(samtools view -c "$out_bam")
    printf "%s\t%s\n" "$name" "$cnt" >> "$filter_counts_file"
    echo "[INFO] $name -> $out_bam (reads: $cnt)"
}

# now we apply the filters
prev_bam="$input_bam"
for entry in "${filters[@]}"; do
    IFS='|' read -r fname fflags fdesc <<< "$entry"
    out_bam="$outdir/filtered_${fname}.bam"

    if [ "$mode" = "cumulative" ]; then
        run_filter "$fname" "$fflags" "$fdesc" "$prev_bam" "$out_bam"
        prev_bam="$out_bam"
    else
        # independent: always start from original input
        run_filter "$fname" "$fflags" "$fdesc" "$input_bam" "$out_bam"
    fi
done

# in the end we produce QC reports
echo "[INFO] All filters applied. Generating QC reports..."
echo "[INFO] Generating QC reports (flagstat, idxstats, stats)..."

samtools flagstat "$input_bam" > "$qc_dir/input_flagstat.txt"
samtools idxstats "$input_bam" > "$qc_dir/input_idxstats.txt"
samtools stats "$input_bam" > "$qc_dir/input_stats.txt"

for bam in "$outdir"/filtered_*.bam; do
    base=$(basename "$bam" .bam)
    samtools flagstat "$bam" > "$qc_dir/${base}_flagstat.txt"
    samtools idxstats "$bam" > "$qc_dir/${base}_idxstats.txt"
    samtools stats "$bam" > "$qc_dir/${base}_stats.txt"
done

# plotting script, with ggplot2, to visualize filter counts, for future implementation, 
# we could switch to a more advanced plotting library if needed
cat > "$outdir/visualize_filters.R" <<'RSCRIPT'
library(ggplot2)
library(scales)

filter_data <- read.table("filter_counts.txt", header=FALSE, sep="\t", col.names=c("Filter","Reads"))
# Keep order as in file
filter_data$Filter <- factor(filter_data$Filter, levels=filter_data$Filter)

p <- ggplot(filter_data, aes(x=Filter, y=Reads)) +
  geom_col() +
  geom_text(aes(label=scales::comma(Reads)), vjust=-0.5, size=3) +
  scale_y_continuous(labels = scales::comma) +
  labs(title="Reads Remaining After Filters",
       x="Filter Applied",
       y="Number of Reads",
       caption="alimap_enhanced.sh output") +
  theme_minimal()

pdf("alimap_filters_plot.pdf", width=8, height=5)
print(p)
dev.off()
RSCRIPT

# (if available and ggplot2 installed)
if command -v Rscript >/dev/null 2>&1; then
    echo "[INFO] Running Rscript to generate barplot (alimap_filters_plot.pdf)"
    (cd "$outdir" && Rscript visualize_filters.R) || echo "[WARN] Rscript failed to create plot. Ensure ggplot2 is installed."
fi

# cleanup temporary files if not keeping them
if [ $keep_temp -eq 0 ]; then
    echo "[INFO] Cleaning up temporary files..."
    rm -rf "$tmpdir"
else
    echo "[INFO] Temporary files kept in $tmpdir for debugging."
fi
echo "[DONE] Pipeline finished. Summary written to $filter_counts_file"
echo "QC files in: $qc_dir"

cat "$filter_counts_file"

exit 0
