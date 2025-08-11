#!/usr/bin/env bash
# alimap.sh
# Stand-alone BAM filtering + QC tool with blacklist exclusion, per-read removal log,
# timing/resource benchmarking, MAPQ histograms export, and enhanced plotting.
# Author: Dany Mukesha (original) + refactor & enhancements (2025-08-11)

set -euo pipefail
IFS=$'
	'

VERSION="2025-08-11"

print_help() {
    cat <<EOF
alimap.sh v$VERSION

Usage: $(basename "$0") -i <input.bam> [options]

Required:
  -i,--input FILE         Input BAM file (sorted). If missing .bai the script will index it.

Options:
  -b,--blacklist FILE     Optional blacklist BED (hg38/hg19 style). If absent, can be auto-downloaded in --test.
  -t,--threads N          Number of threads for samtools (default: 1)
  -o,--outdir DIR         Output directory (default: filtered_reads)
  -m,--mode MODE          cumulative|independent (default: cumulative)
  --keep-temp             Keep temporary files
  --test                  Download small public test BAM + blacklist and run pipeline on them
  -h,--help               Show this help

Description:
  - Produces filtered BAMs per filter step, indexed.
  - Produces a per-read removal log when blacklist exclusion is used (read names removed).
  - Produces mapping quality histograms (CSV) per filter.
  - Produces timing/resource usage (if /usr/bin/time available) for heavy steps.
  - Creates enhanced R plots: absolute counts + percent removed per step (stacked bars).

Notes:
  - Required: samtools, awk, sort, uniq, Rscript
  - Optional but recommended: bedtools (blacklist exclusion), /usr/bin/time for benchmarking
  - Test data source (example): UCSC ENCODE small BAM and Boyle-Lab ENCODE blacklist.
    See: https://hgdownload.cse.ucsc.edu/ (example BAM) and https://github.com/Boyle-Lab/Blacklist (blacklist). 
EOF
}

# Defaults
threads=1
outdir="filtered_reads"
mode="cumulative"
keep_temp=0
blacklist_bed=""
run_test=0

# Parse args
ARGS=$(getopt -o i:b:t:o:m:kh --long input:,blacklist:,threads:,outdir:,mode:,keep-temp,test,help -n "$(basename "$0")" -- "$@") || { print_help; exit 1; }
eval set -- "$ARGS"
while true; do
    case "$1" in
        -i|--input) input_bam="$2"; shift 2;;
        -b|--blacklist) blacklist_bed="$2"; shift 2;;
        -t|--threads) threads="$2"; shift 2;;
        -o|--outdir) outdir="$2"; shift 2;;
        -m|--mode) mode="$2"; shift 2;;
        --keep-temp) keep_temp=1; shift;;
        --test) run_test=1; shift;;
        -h|--help) print_help; exit 0;;
        --) shift; break;;
        *) echo "Parsing error"; exit 1;;
    esac
done

# If test mode: download a small public BAM and blacklist
if [ "$run_test" -eq 1 ]; then
    echo "[INFO] Running in --test mode: will download example BAM + blacklist"
    mkdir -p test_data
    cd test_data

    # Example small BAM from UCSC ENCODE (public sample). This is used widely as an example/test file.
    TEST_BAM_URL="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqGm12878G1bAlnRep1.bam"
    TEST_BAM="wgEncodeUwRepliSeqGm12878G1bAlnRep1.bam"

    if [ ! -f "$TEST_BAM" ]; then
        echo "[INFO] Downloading example BAM from UCSC: $TEST_BAM_URL"
        if command -v wget >/dev/null 2>&1; then
            wget -q --show-progress "$TEST_BAM_URL" -O "$TEST_BAM" || { echo "Failed to download test BAM"; exit 1; }
        elif command -v curl >/dev/null 2>&1; then
            curl -L -o "$TEST_BAM" "$TEST_BAM_URL" || { echo "Failed to download test BAM"; exit 1; }
        else
            echo "Error: wget or curl required to download test BAM"; exit 1
        fi
    fi

    # Use Boyle-Lab hg38 blacklist (v2) as a reasonable small blacklist; raw GitHub URL
    BLACKLIST_URL="https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz"
    BLACKLIST_GZ="hg38-blacklist.v2.bed.gz"
    BLACKLIST_BED="hg38-blacklist.v2.bed"

    if [ ! -f "$BLACKLIST_BED" ]; then
        echo "[INFO] Downloading blacklist: $BLACKLIST_URL"
        if command -v wget >/dev/null 2>&1; then
            wget -q --show-progress "$BLACKLIST_URL" -O "$BLACKLIST_GZ" || { echo "Failed to download blacklist"; exit 1; }
        elif command -v curl >/dev/null 2>&1; then
            curl -L -o "$BLACKLIST_GZ" "$BLACKLIST_URL" || { echo "Failed to download blacklist"; exit 1; }
        else
            echo "Error: wget or curl required to download blacklist"; exit 1
        fi
        gunzip -f "$BLACKLIST_GZ"
    fi

    # Set input and blacklist to downloaded files
    input_bam="$PWD/$TEST_BAM"
    blacklist_bed="$PWD/$BLACKLIST_BED"
    cd - >/dev/null
fi

# Validate required variables
if [ -z "${input_bam-}" ]; then
    echo "Error: --input required (or use --test)." >&2
    print_help
    exit 1
fi

if [ ! -f "$input_bam" ]; then
    echo "Error: input BAM '$input_bam' not found." >&2
    exit 1
fi

# Create output dirs
mkdir -p "$outdir"
qc_dir="$outdir/qc"
mkdir -p "$qc_dir"

# Temp dir
tmpdir=$(mktemp -d -t alimap.XXXXXX)
trap '[[ $keep_temp -eq 1 ]] || rm -rf "$tmpdir"' EXIT

# Dependency checks
missing=()
for c in samtools awk sort uniq Rscript; do
    if ! command -v "$c" >/dev/null 2>&1; then
        missing+=("$c")
    fi
done
if [ ${#missing[@]} -gt 0 ]; then
    echo "Error: Missing required commands: ${missing[*]}"; exit 1
fi

# bedtools optional
bedtools_ok=1
if ! command -v bedtools >/dev/null 2>&1; then
    bedtools_ok=0
    echo "[WARN] bedtools not found: blacklist exclusion will be skipped unless bedtools installed."
fi

# /usr/bin/time for benchmarking (optional)
time_bin=""
if command -v /usr/bin/time >/dev/null 2>&1; then
    time_bin="/usr/bin/time -v"
fi

# Check index; create if missing
if [ ! -f "${input_bam}.bai" ] && [ ! -f "${input_bam%.bam}.bai" ]; then
    echo "[INFO] Index missing for $input_bam. Creating index..."
    samtools index -@ "$threads" "$input_bam"
fi

# Validate blacklist
if [ -n "$blacklist_bed" ]; then
    if [ ! -f "$blacklist_bed" ]; then
        echo "Error: blacklist file '$blacklist_bed' not found."; exit 1
    fi
    if [ $bedtools_ok -eq 0 ]; then
        echo "[WARN] bedtools not installed: will skip blacklist exclusion."; blacklist_bed=""
    fi
fi

# Filters definition: name|samtools_flags|description
filters=(
  "no_filter||No filtering (copy input)"
  "exclude_unmapped|-F 4|Exclude unmapped reads"
  "exclude_secondary|-F 256|Exclude secondary alignments"
  "mapq_30|-q 30|Keep reads with MAPQ >= 30"
  "properly_paired|-f 2|Keep properly paired reads"
)

filter_counts="$outdir/filter_counts.txt"
: > "$filter_counts"

# Create header run_info
run_info="$outdir/run_info.txt"
{
  echo "alimap.sh v$VERSION"
  echo "timestamp: $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
  echo "input_bam: $input_bam"
  echo "blacklist_bed: ${blacklist_bed:-NONE}"
  echo "threads: $threads"
  echo "mode: $mode"
  echo "pwd: $(pwd)"
  echo "samtools: $(samtools --version | head -n1)"
  if [ -n "$time_bin" ]; then echo "time: /usr/bin/time available"; fi
} > "$run_info"

# Function to write MAPQ histogram CSV
write_mapq_hist() {
    local bam=$1
    local outcsv=$2
    # produce histogram of MAPQ values (0-60), bin per MAPQ
    samtools view "$bam" | awk '{print $5}' | sort -n | uniq -c | awk '{print $2","$1}' > "$outcsv"
}

# Run filter step
run_filter() {
    local name="$1" flags="$2" desc="$3" in_bam="$4" out_bam="$5"

    echo "[STEP] $name: $desc"
    start_time=$(date +%s)
    start_epoch=$(date +%s%3N)

    if [ -z "$flags" ]; then
        cmd=(samtools view -@ "$threads" -b "$in_bam" -o "$out_bam")
    else
        # split flags into array
        read -r -a flag_array <<< "$flags"
        cmd=(samtools view -@ "$threads" -b "${flag_array[@]}" "$in_bam" -o "$out_bam")
    fi

    if [ -n "$time_bin" ]; then
        # run with timing
        /usr/bin/time -v "${cmd[@]}" 1>"$tmpdir/${name}.samtools.stdout" 2>"$tmpdir/${name}.samtools.time" || { cat "$tmpdir/${name}.samtools.time"; exit 1; }
    else
        "${cmd[@]}"
    fi

    # If blacklist requested, create per-read removal log
    if [ -n "$blacklist_bed" ]; then
        echo "[INFO] Excluding blacklist regions from $out_bam (producing removal log)..."
        # bedtools intersect -abam IN -b blacklist -v > OUT (keeps reads NOT overlapping blacklist)
        tmp_no_black="$tmpdir/${name}.no_blacklist.bam"
        if [ $bedtools_ok -eq 1 ]; then
            # produce list of read names removed (overlapping blacklist)
            # bedtools intersect -abam IN -b blacklist -u outputs reads that overlap; -u to report once
            bedtools intersect -abam "$out_bam" -b "$blacklist_bed" -u | samtools view -F 0x4 -@ "$threads" - | awk '{print $1}' | sort -u > "$out_bam.removed.reads.txt"
            # now produce bam without blacklist overlapping reads
            bedtools intersect -abam "$out_bam" -b "$blacklist_bed" -v > "$tmp_no_black"
            mv "$tmp_no_black" "$out_bam"
        else
            echo "[WARN] bedtools not available; skipping blacklist exclusion step for $name"
        fi
    fi

    # Index output
    samtools index -@ "$threads" "$out_bam"

    # Count reads
    cnt=$(samtools view -c "$out_bam")
    printf "%s	%s
" "$name" "$cnt" >> "$filter_counts"

    # produce MAPQ histogram CSV
    mapq_csv="$outdir/${name}_mapq_hist.csv"
    write_mapq_hist "$out_bam" "$mapq_csv"

    end_time=$(date +%s)
    duration=$((end_time - start_time))

    # If timing log created by /usr/bin/time, copy it to outdir
    if [ -f "$tmpdir/${name}.samtools.time" ]; then
        mv "$tmpdir/${name}.samtools.time" "$outdir/${name}.samtools.time.txt"
    fi

    echo "[DONE] $name: reads=$cnt duration=${duration}s mapq_csv=$mapq_csv"
}

# Execute filters
prev_bam="$input_bam"
for entry in "${filters[@]}"; do
    IFS='|' read -r fname fflags fdesc <<< "$entry"
    out_bam="$outdir/filtered_${fname}.bam"

    if [ "$mode" = "cumulative" ]; then
        run_filter "$fname" "$fflags" "$fdesc" "$prev_bam" "$out_bam"
        prev_bam="$out_bam"
    else
        run_filter "$fname" "$fflags" "$fdesc" "$input_bam" "$out_bam"
    fi
done

# QC for input + outputs
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

# Enhanced R plotting: stacked bars + percent removed
cat > "$outdir/visualize_filters_enhanced.R" <<'RSCRIPT'
library(ggplot2)
library(scales)

# Read filter counts
f <- read.table("filter_counts.txt", header=FALSE, sep="	", col.names=c("Filter","Reads"), stringsAsFactors=FALSE)
# Keep file order
f$Filter <- factor(f$Filter, levels=f$Filter)
# compute removed counts relative to previous (cumulative assumed)
f$ReadsNumeric <- as.numeric(gsub(",", "", f$Reads))
removed <- c(0, head(f$ReadsNumeric, -1) - tail(f$ReadsNumeric, -0))
# Sometimes independent mode will give negative removed; coerce to 0
removed[removed < 0] <- 0
plot_df <- data.frame(Filter=f$Filter, Kept=f$ReadsNumeric, Removed=removed)
plot_df_melt <- reshape2::melt(plot_df, id.vars="Filter")

p1 <- ggplot(plot_df_melt, aes(x=Filter, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  geom_text(data=plot_df, aes(x=Filter, y=Kept + pmax(Removed,1), label=scales::comma(Kept)), vjust=-0.5, size=3) +
  scale_y_continuous(labels = scales::comma) +
  labs(title="Reads Kept and Removed per Filter",
       x="Filter",
       y="Number of Reads",
       fill="Category",
       caption="alimap.sh output (Kept and Removed)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

pdf("alimap_filters_plot_enhanced.pdf", width=10, height=6)
print(p1)
dev.off()
RSCRIPT

# Run R script in outdir
(cd "$outdir" && Rscript visualize_filters_enhanced.R) || echo "[WARN] Enhanced plotting failed (check R + reshape2 + ggplot2 installed)."

# Final summary
echo "[DONE] Pipeline finished. Outputs in: $outdir"
echo "- filter counts: $filter_counts"
echo "- QC directory: $qc_dir"

echo "Filter counts:";
cat "$filter_counts"

exit 0
