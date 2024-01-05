#!/bin/bash
# Created on Wed Nov 15 11:45:07 2023
# @author: Dany Mukesha

# Check if input BAM file is provided
if [ -z "$1" ]; then
    echo "Usage: ./alimap_filters.sh <input_bam>"
    exit 1
fi

# Input BAM file
input_bam="$1"

# Output directory
output_dir="filtered_reads"

# Make sure the output directory exists
mkdir -p "$output_dir"

# Define a function to filter and count reads
filter_and_count_reads() {
    filter="$1"
    output_bam="$output_dir/filtered_reads_${filter// /}.bam"  # Remove spaces from filter for file naming

    # Filter BAM file
    samtools view -b $filter -o "$output_bam" "$input_bam"

    # Count the number of reads in the filtered BAM
    count=$(samtools view -c "$output_bam")

    # Print filter and count
    echo "Filter: $filter, Remaining Reads: $count"
    
    # Print the filter name and count to a file
    echo "$filter: $count" >> filter_counts.txt 
}

# Apply filters and count reads
filters=(
    ""           # No additional filter
    "-F 4"       # Exclude unmapped reads
    "-F 256"     # Exclude secondary alignments
    "-q 30"      # Minimum mapping quality of 30
    "-f 2"       # Properly paired reads
)

# Remove the file if exists
[ -e filter_counts.txt ] && rm filter_counts.txt
    
# Apply filters and count reads
for filter in "${filters[@]}"; do
    filter_and_count_reads "$filter"
done

# R script for plotting
cat << 'RSCRIPT' > visualize_filters.R
library(ggplot2)

# Read filter names and counts from the file
filter_data <- read.table("filter_counts.txt", sep=":", col.names=c("Filter", "Reads"))

# Define custom colors
my_colors <- c("#4e79a7", "#e15759", "#76b7b2", "#ff9da7", "#9c755f")

# Plot
p <- ggplot(filter_data, aes(x=Filter, y=Reads, fill=Filter)) +
  geom_bar(stat="identity", color="black") +
  labs(title="Number of Reads After Filtering",
       x="Filter",
       y="Number of Reads",
       caption="Filters applied to aligned_reads.bam") +
  theme_minimal() +
  theme(legend.position="none") +  # Remove legend
  scale_fill_manual(values=my_colors)  # Use custom colors

pdf("alimap_filters_plot.pdf")
print(p)
dev.off()

RSCRIPT

# Run the R script for visualization
Rscript visualize_filters.R
