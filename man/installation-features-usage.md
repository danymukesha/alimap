## Installation

1. **Dependencies:**
   - [samtools](http://www.htslib.org/)
   - [R](https://www.r-project.org/)
   - [ggplot2](https://ggplot2.tidyverse.org/)

2. **Clone the Repository:**
   ```bash
   git clone https://github.com/danymukesha/alimap.git
   cd alimap

## Features

- **Filtering:** The tool applies various filters to the input BAM file, such as excluding unmapped reads, secondary alignments, setting a minimum mapping quality, and selecting properly paired reads.
  
- **Visualization:** It generates a bar plot visualizing the number of reads remaining after each filtering step using the ggplot2 library in R.

- **Flexibility:** It permits to easily customize the filters applied and adjust the R visualization script to suit your analysis needs.

## Usage

```bash
chmod +x alimap_filters.sh
./alimap_filters.sh <input_bam>
```

- `<input_bam>`: Path to the input BAM file.

### Example 1:

```bash
./alimap_filters.sh input_aligned.bam
```

### Example 2:

```bash
./alimap_filters.sh <path>/<to>/<the>/<name>.bam
```

## Output

The tool generates the following output files and directories:

- **Filtered Reads Directory:**
  - Directory: `filtered_reads`
  - Contains filtered BAM files after each filtering step.

- **Filter Counts File:**
  - File: `filter_counts.txt`
  - Records the read counts after each filtering step.

- **Visualization Plot:**
  - File: `alimap_filters.pdf`
  - A bar plot displaying the number of reads after each filtering step.

- **R script**
  - file: `visualize_filters.R`
  - he R script uses ggplot2 to create a bar plot displaying the number of reads after each filtering step.
  
### Example Output:

```bash
Filter: , Remaining Reads: 2599298  # No additional filter
Filter: -F 4, Remaining Reads: 1816425  # Exclude unmapped reads
Filter: -F 256, Remaining Reads: 2599298  # Exclude secondary alignments
Filter: -q 30, Remaining Reads: 1459959  # Minimum mapping quality of 30
Filter: -f 2, Remaining Reads: 1688814  # Properly paired reads
```

## Summary

The tool generates a directory named `filtered_reads` containing filtered BAM files, and a file named `filter_counts.txt` recording the read counts after each filtering step. The final visualization is stored as `alimap_filters.pdf`.
