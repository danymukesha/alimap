#!/usr/bin/env bash
# alimap.sh
# Stand-alone BAM filtering + QC + benchmarking tool
# Version: 2025-08-11 Dany MUKESHA
#
# Features:
#  - cumulative/independent modes
#  - optional --test to download small public BAM + blacklists
#  - per-read removal log for blacklist exclusion
#  - /usr/bin/time resource usage capture (if available)
#  - per-filter MAPQ histograms (CSV)
#  - enhanced stacked-bar PDF (kept vs removed + percent)
#  - auto blacklist selection (hg19 vs hg38) by contig naming
#  - optional checksum verification support (pass --checksum checksum_file)
#
# NOTE: If you see "/usr/bin/env: ‘bash\r’: No such file or directory" convert line endings:
#       sed -i 's/\r$//' alimap.sh

set -euo pipefail
IFS=$'\n\t'

VERSION="2025-08-11"

print_help() {
  cat <<EOF
alimap.sh v$VERSION (Author: Dany MUKESHA)

Usage: $(basename "$0") -i <input.bam> [options]

Required:
  -i,--input FILE         Input BAM file (sorted). If missing .bai the script will index it.

Options:
  -b,--blacklist FILE     Optional blacklist BED (if absent and --test used, script downloads hg19/hg38 blacklists)
  -t,--threads N          Number of threads for samtools (default: 1)
  -o,--outdir DIR         Output directory (default: filtered_reads)
  -m,--mode MODE          cumulative|independent (default: cumulative)
  --keep-temp             Keep temporary files
  --test                  Download small public example BAM and blacklists and run pipeline using them
  --checksum FILE         Optional: path to a checksum file (tab-delimited: filename<TAB>md5) to validate downloads
  -h,--help               Show this help

Outputs (written to --outdir):
  - filtered_*.bam (.bai generated)
  - filter_counts.txt    (tab: filter<TAB>reads)
  - *_mapq_hist.csv      (per-filter mapping quality histogram)
  - *.removed.reads.txt  (per-filter list of read names removed by blacklist exclusion)
  - *.samtools.time.txt  (timing/resource usage from /usr/bin/time, if available)
  - qc/*                 (samtools flagstat/idxstats/stats)
  - alimap_filters_plot.pdf (stacked bars: kept vs removed + percent)
  - run_info.txt         (provenance)
EOF
}

threads=1
outdir="filtered_reads"
mode="cumulative"
keep_temp=0
blacklist_bed=""
run_test=0
checksum_file=""

ARGS=$(getopt -o i:b:t:o:m:kh --long input:,blacklist:,threads:,outdir:,mode:,keep-temp,test,checksum:,help -n "$(basename "$0")" -- "$@") || { print_help; exit 1; }
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
    --checksum) checksum_file="$2"; shift 2;;
    -h|--help) print_help; exit 0;;
    --) shift; break;;
    *) echo "Parsing error"; exit 1;;
  esac
done

# --- Test-mode downloads (small real files) ---
if [ "${run_test:-0}" -eq 1 ]; then
  echo "[INFO] --test: will download example BAM + blacklists into ./test_data"
  mkdir -p test_data
  cd test_data

  # small BAM (public UCSC ENCODE example). If inaccessible, replace with another public small BAM.
  TEST_BAM_URL="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqGm12878G1bAlnRep1.bam"
  TEST_BAM="wgEncodeUwRepliSeqGm12878G1bAlnRep1.bam"

  if [ ! -f "$TEST_BAM" ]; then
    echo "[INFO] Downloading example BAM..."
    if command -v wget >/dev/null 2>&1; then
      wget -c "$TEST_BAM_URL" -O "$TEST_BAM" || { echo "Download failed"; exit 1; }
    elif command -v curl >/dev/null 2>&1; then
      curl -L "$TEST_BAM_URL" -o "$TEST_BAM" || { echo "Download failed"; exit 1; }
    else
      echo "Error: wget or curl required to download test BAM"; exit 1
    fi
  fi

  # Boyle-Lab blacklists (hg19/hg38)
  BL_HG19_URL="https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz"
  BL_HG38_URL="https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz"
  BL_HG19_GZ="hg19-blacklist.v2.bed.gz"
  BL_HG38_GZ="hg38-blacklist.v2.bed.gz"

  if [ ! -f "hg19-blacklist.v2.bed" ]; then
    echo "[INFO] Downloading hg19 blacklist..."
    wget -c "$BL_HG19_URL" -O "$BL_HG19_GZ" && gunzip -f "$BL_HG19_GZ" || { echo "Failed to get hg19 blacklist"; exit 1; }
  fi
  if [ ! -f "hg38-blacklist.v2.bed" ]; then
    echo "[INFO] Downloading hg38 blacklist..."
    wget -c "$BL_HG38_URL" -O "$BL_HG38_GZ" && gunzip -f "$BL_HG38_GZ" || { echo "Failed to get hg38 blacklist"; exit 1; }
  fi

  input_bam="$PWD/$TEST_BAM"
  # user didn't specify blacklist; will auto-detect later
  cd - >/dev/null
fi

# basic validations
if [ -z "${input_bam-}" ]; then
  echo "Error: --input is required (or use --test)."; print_help; exit 1
fi
if [ ! -f "$input_bam" ]; then echo "Error: input BAM not found: $input_bam"; exit 1; fi

mkdir -p "$outdir"
qc_dir="$outdir/qc"
mkdir -p "$qc_dir"

tmpdir=$(mktemp -d -t alimap.XXXXXX)
trap '[[ $keep_temp -eq 1 ]] || rm -rf "$tmpdir"' EXIT

# dependencies
req=(samtools awk sort uniq Rscript)
missing=()
for c in "${req[@]}"; do command -v "$c" >/dev/null 2>&1 || missing+=("$c"); done
if [ ${#missing[@]} -gt 0 ]; then echo "Error: missing required commands: ${missing[*]}"; exit 1; fi

# optional tools
bedtools_ok=1
if ! command -v bedtools >/dev/null 2>&1; then bedtools_ok=0; echo "[WARN] bedtools not found; blacklist features will be skipped unless bedtools installed."; fi

time_bin=""
if command -v /usr/bin/time >/dev/null 2>&1; then time_bin="/usr/bin/time -v"; fi

# create index if missing
if [ ! -f "${input_bam}.bai" ] && [ ! -f "${input_bam%.bam}.bai" ]; then
  echo "[INFO] BAM index missing — creating index with samtools index"
  samtools index -@ "$threads" "$input_bam"
fi

# checksum verification function
verify_checksums_if_provided() {
  local chfile="$1"
  if [ -z "$chfile" ]; then return 0; fi
  if [ ! -f "$chfile" ]; then echo "[WARN] checksum file not found: $chfile"; return 0; fi
  echo "[INFO] Verifying checksums from: $chfile"
  # expected format: filename<TAB>md5
  while IFS=$'\t' read -r fname md5; do
    if [ -z "$fname" ] || [ -z "$md5" ]; then continue; fi
    if [ ! -f "$fname" ]; then echo "[ERROR] checksum file references missing file: $fname"; continue; fi
    got=$(md5sum "$fname" | awk '{print $1}')
    if [ "$got" != "$md5" ]; then
      echo "[ERROR] MD5 mismatch for $fname: expected $md5 got $got"
    else
      echo "[OK] $fname matches MD5"
    fi
  done < "$chfile"
}

# run provenance header
run_info="$outdir/run_info.txt"
{
  echo "alimap.sh v$VERSION"
  echo "timestamp: $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
  echo "input_bam: $input_bam"
  echo "blacklist_bed (user): ${blacklist_bed:-NONE}"
  echo "threads: $threads"
  echo "mode: $mode"
  echo "pwd: $(pwd)"
  echo "samtools: $(samtools --version | head -n1)"
  [ -n "$time_bin" ] && echo "time: /usr/bin/time available"
} > "$run_info"

# If checksum file provided, run verify (useful if user passed checksums for downloaded files)
if [ -n "${checksum_file}" ]; then
  verify_checksums_if_provided "$checksum_file" || true
fi

# If --test downloaded files, record sizes & md5s
if [ "${run_test:-0}" -eq 1 ] && [ -d "test_data" ]; then
  echo "[INFO] Recording MD5 and sizes for downloaded files (test mode)"
  (cd test_data && ls -lh && md5sum * | tee "$outdir/test_downloads.md5") || true
fi

# detect contig naming in BAM
detect_genome_from_bam() {
  local bam=$1
  sn=$(samtools view -H "$bam" | awk '/^@SQ/ {print $2; exit}' | sed 's/SN://') || true
  if [ -z "$sn" ]; then echo "unknown"; return; fi
  if [[ "$sn" =~ ^chr ]]; then echo "chr_prefixed"; else echo "no_chr"; fi
}

# auto select blacklist if not provided
if [ -z "$blacklist_bed" ]; then
  detect=$(detect_genome_from_bam "$input_bam")
  echo "[INFO] Detected contig style: $detect"
  if [ "$detect" = "chr_prefixed" ]; then
    if [ -f "./test_data/hg38-blacklist.v2.bed" ]; then
      blacklist_bed="./test_data/hg38-blacklist.v2.bed"
      echo "[INFO] Auto-selected blacklist: $blacklist_bed"
    elif [ -f "./test_data/hg19-blacklist.v2.bed" ]; then
      blacklist_bed="./test_data/hg19-blacklist.v2.bed"
      echo "[INFO] Fallback selected blacklist: $blacklist_bed"
    fi
  elif [ "$detect" = "no_chr" ]; then
    if [ -f "./test_data/hg19-blacklist.v2.bed" ]; then
      blacklist_bed="./test_data/hg19-blacklist.v2.bed"
      echo "[INFO] Auto-selected blacklist: $blacklist_bed"
    elif [ -f "./test_data/hg38-blacklist.v2.bed" ]; then
      blacklist_bed="./test_data/hg38-blacklist.v2.bed"
      echo "[INFO] Fallback selected blacklist: $blacklist_bed"
    fi
  else
    echo "[WARN] Could not detect contig naming. If you want blacklist exclusion, provide --blacklist PATH"
  fi
fi

if [ -n "$blacklist_bed" ] && [ ! -f "$blacklist_bed" ]; then
  echo "Error: specified blacklist not found: $blacklist_bed"; exit 1
fi
if [ -n "$blacklist_bed" ] && [ $bedtools_ok -eq 0 ]; then
  echo "[WARN] blacklist requested but bedtools not found; blacklist will be skipped"; blacklist_bed=""
fi

# Filters: name|samtools_flags|description
filters=(
  "no_filter||No filtering (copy input)"
  "exclude_unmapped|-F 4|Exclude unmapped reads"
  "exclude_secondary|-F 256|Exclude secondary alignments"
  "mapq_30|-q 30|Keep reads with MAPQ >= 30"
  "properly_paired|-f 2|Keep properly paired reads"
)

filter_counts="$outdir/filter_counts.txt"
: > "$filter_counts"

# mapq histogram writer
write_mapq_hist() {
  local bam=$1
  local outcsv=$2
  samtools view "$bam" | awk '{print $5}' | sort -n | uniq -c | awk '{print $2","$1}' > "$outcsv" || echo "0,0" > "$outcsv"
}

# run filter step
run_filter() {
  local name="$1" flags="$2" desc="$3" in_bam="$4" out_bam="$5"
  echo "[STEP] $name — $desc"
  start_ts=$(date +%s)

  # build samtools command
  if [ -z "$flags" ]; then
    cmd=(samtools view -@ "$threads" -b "$in_bam" -o "$out_bam")
  else
    read -r -a flag_array <<< "$flags"
    cmd=(samtools view -@ "$threads" -b "${flag_array[@]}" "$in_bam" -o "$out_bam")
  fi

  # run with /usr/bin/time if available
  if [ -n "$time_bin" ]; then
    /usr/bin/time -v "${cmd[@]}" 1>"$tmpdir/${name}.samtools.stdout" 2>"$tmpdir/${name}.samtools.time" || { cat "$tmpdir/${name}.samtools.time"; exit 1; }
    mv "$tmpdir/${name}.samtools.time" "$outdir/${name}.samtools.time.txt"
  else
    "${cmd[@]}"
  fi

  # Blacklist exclusion (per-read removal log)
  if [ -n "$blacklist_bed" ] && [ $bedtools_ok -eq 1 ]; then
    echo "[INFO] Removing reads overlapping blacklist for step $name"
    bedtools intersect -abam "$out_bam" -b "$blacklist_bed" -u | samtools view -@ "$threads" -F 0x4 - | awk '{print $1}' | sort -u > "$outdir/${name}.removed.reads.txt"
    bedtools intersect -abam "$out_bam" -b "$blacklist_bed" -v > "$tmpdir/${name}.no_blacklist.bam"
    mv "$tmpdir/${name}.no_blacklist.bam" "$out_bam"
  fi

  # Index
  samtools index -@ "$threads" "$out_bam"

  # Count and write
  cnt=$(samtools view -c "$out_bam")
  printf "%s\t%s\n" "$name" "$cnt" >> "$filter_counts"

  # Mapq hist
  write_mapq_hist "$out_bam" "$outdir/${name}_mapq_hist.csv"

  end_ts=$(date +%s)
  dur=$((end_ts - start_ts))
  echo "[DONE] $name : reads=$cnt duration=${dur}s mapq_csv=${name}_mapq_hist.csv"
}

# iterate filters
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

echo "[INFO] Running samtools flagstat/idxstats/stats"
samtools flagstat "$input_bam" > "$qc_dir/input_flagstat.txt"
samtools idxstats "$input_bam" > "$qc_dir/input_idxstats.txt"
samtools stats "$input_bam" > "$qc_dir/input_stats.txt"
for bam in "$outdir"/filtered_*.bam; do
  base=$(basename "$bam" .bam)
  samtools flagstat "$bam" > "$qc_dir/${base}_flagstat.txt"
  samtools idxstats "$bam" > "$qc_dir/${base}_idxstats.txt"
  samtools stats "$bam" > "$qc_dir/${base}_stats.txt"
done

cat > "$outdir/isualize_filters.R" <<'RSCRIPT'
library(ggplot2)
library(scales)
library(reshape2)

f <- read.table("filter_counts.txt", header=FALSE, sep="\t", col.names=c("Filter","Reads"), stringsAsFactors=FALSE)
f$Filter <- factor(f$Filter, levels=f$Filter)
f$ReadsNumeric <- as.numeric(gsub(",","", f$Reads))
removed <- c(0, head(f$ReadsNumeric, -1) - tail(f$ReadsNumeric, -0))
removed[removed < 0] <- 0
plot_df <- data.frame(Filter=f$Filter, Kept=f$ReadsNumeric, Removed=removed)
plot_df_melt <- melt(plot_df, id.vars="Filter")
plot_df$Total <- plot_df$Kept + plot_df$Removed
plot_df$PercentRemoved <- ifelse(plot_df$Total>0, round(100 * plot_df$Removed / plot_df$Total, 2), 0)

p <- ggplot(plot_df_melt, aes(x = Filter, y = value, fill = variable)) +
    geom_bar(stat = 'identity') +
    geom_text(data = plot_df, inherit.aes = FALSE, 
        aes(x = Filter, y = Kept + Removed + 1, 
            label = paste0(PercentRemoved, "% removed")), 
        size = 3, angle = 45, hjust = 0) +
    labs(title='Reads Kept and Removed per Filter', 
        x = 'Filter', y = 'Reads', fill = 'Category') +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(labels = scales::comma)

pdf('alimap_filters_plot.pdf', width=10, height = 6)
print(p)
dev.off()
RSCRIPT

missing_rpkgs=()
for pkg in ggplot2 reshape2 scales; do
  if ! Rscript -e "if (!requireNamespace('$pkg', quietly=TRUE)) quit(status=1)" >/dev/null 2>&1; then
    missing_rpkgs+=("$pkg")
  fi
done
if [ ${#missing_rpkgs[@]} -gt 0 ]; then
  echo "[WARN] R packages missing: ${missing_rpkgs[*]}. The plot step may fail. Install in R: install.packages(c(${missing_rpkgs[@]}))"
fi

(cd "$outdir" && Rscript isualize_filters.R) || echo "[WARN] Plotting failed - ensure R + ggplot2/reshape2/scales installed"

echo "[DONE] Pipeline finished. Outputs in $outdir"
cat "$filter_counts"
exit 0
