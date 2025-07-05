We eh new enhancements, this tools underwent this updates:

* Improved filtering with blacklist exclusion
* Multithreading support
* Expanded QC stats
* Clean directory structure
* Complete README with usage
* Example `hg38_blacklist.bed` download command
* R visualization enhanced
* Sample commands for preparing a toy BAM


**alimap** as a command-line tool will now be used to apply common BAM filtering steps
 and visualize read counts after filtering. It includes blacklist region exclusion, 
 multi-threading support, and expanded QC reports.

## all features

- Filter BAM by:
  - Excluding unmapped reads
  - Excluding secondary alignments
  - Minimum mapping quality (default ≥30)
  - Selecting properly paired reads
- Optional blacklist region exclusion using BED file
- Multi-threading for faster BAM processing
- QC reports with `samtools flagstat` and `idxstats`
- Publication-quality bar plot with R/ggplot2 showing read counts
- Easy customization of filters and blacklist input

## installations

### depending one

- samtools (≥1.9 recommended)
- bedtools (optional, for blacklist exclusion)
- R with `ggplot2`
- bash shell

## new usage

```bash
./alimap_enhanced.sh <input_bam> [--blacklist <blacklist.bed>] [--threads <num>]
```

example with blacklist and 4 threads:

```bash
./alimap_enhanced.sh input.bam --blacklist hg38_blacklist.bed --threads 4
```

## outputs

* `filtered_reads/`: filtered BAM files and their indexes
* `qc/`: QC reports (`flagstat` and `idxstats`)
* `filter_counts.txt`: counts of reads after each filter
* `alimap_filters_plot.pdf`: barplot of filtered read counts

## obtain Blacklist BED for hg38

The ENCODE blacklist BED file for hg38 can be downloaded:

```bash
wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz
gunzip hg38.blacklist.bed.gz
mv hg38.blacklist.bed hg38_blacklist.bed
```

## Preparing a Toy BAM for Testing

Use Gencode primary assembly:

```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
mv GRCh38.primary_assembly.genome.fa reference.fa
samtools faidx reference.fa

# Simulate reads (requires wgsim)
wgsim -N 1000 -1 reads_1.fq -2 reads_2.fq reference.fa /dev/null

# Align with bwa mem
bwa index reference.fa
bwa mem -t 4 reference.fa reads_1.fq reads_2.fq > aln.sam

samtools view -bS aln.sam | samtools sort -o input.bam
samtools index input.bam
```

### 1. **Illumina Platinum Genomes - NA12878 (CEPH/Utah pedigree)**

* Sample: NA12878 (well-studied human genome sample)
* Sequencing: Illumina HiSeq paired-end reads
* Size: \~30X coverage (you can download subsets for testing)

#### Download Example Subset (2x100bp paired-end reads)

You can download smaller fastq subsets from **EBI ENA** or **SRA**.

**From ENA (European Nucleotide Archive):**

```bash
# Read 1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194146/ERR194146_1.fastq.gz

# Read 2
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194146/ERR194146_2.fastq.gz
```

These are part of the NA12878 dataset (ERR194146).

---

### 2. **1000 Genomes Project - Sample NA12878**

From SRA with fastq dump or from ENA:

```bash
# Read 1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz

# Read 2
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz
```

---

### 3. **Example from ENCODE Project**

ENCODE offers many publicly accessible RNA-seq and DNA-seq datasets with raw FASTQ files.

Example DNA-seq sample:

```bash
wget https://www.encodeproject.org/files/ENCFF001UUS/@@download/ENCFF001UUS.fastq.gz
wget https://www.encodeproject.org/files/ENCFF001UUV/@@download/ENCFF001UUV.fastq.gz
```

(Replace with your sample of interest from [https://www.encodeproject.org/](https://www.encodeproject.org))

---

### Quick Setup to Download NA12878 Reads for Testing

```bash
mkdir -p test_reads && cd test_reads

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194146/ERR194146_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194146/ERR194146_2.fastq.gz

# Check files
gunzip -c ERR194146_1.fastq.gz | head -n 8
gunzip -c ERR194146_2.fastq.gz | head -n 8
```

---

### After Download: Align Reads (Example with bwa mem)

```bash
# Index your reference (if not done)
bwa index reference.fa

# Align
bwa mem -t 4 reference.fa ERR194146_1.fastq.gz ERR194146_2.fastq.gz > aln.sam

# Convert and sort
samtools view -bS aln.sam | samtools sort -@4 -o input.bam
samtools index input.bam
```

Now `input.bam` is ready to run with your enhanced `alimap_enhanced.sh`!
---

# 3. Example Blacklist BED Download (already in README, repeat here)

```sh
wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz
gunzip hg38.blacklist.bed.gz
mv hg38.blacklist.bed hg38_blacklist.bed
````



# with small sized samples:
---

### 1. **ENCODE Pilot Datasets - Small Test FASTQ**

ENCODE provides many test FASTQ files, some are very small:

* Example: **ENCFF001TDO** (Read 1, \~2.7 MB)
* Example: **ENCFF001TDP** (Read 2, \~2.7 MB)

Download:

```bash
wget https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.fastq.gz
wget https://www.encodeproject.org/files/ENCFF001TDP/@@download/ENCFF001TDP.fastq.gz
```

Perfect for quick testing!

---

### 2. **SRA - Small Subsets (1000 Genomes Project)**

A few hundred thousand reads (\~10 MB):

* ERR031844 (low coverage subset of NA12878)

Download:

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR031/ERR031844/ERR031844_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR031/ERR031844/ERR031844_2.fastq.gz
```

---

### 3. **Simulated Small Dataset**

If you prefer synthetic data, `wgsim` (from samtools) can generate small paired-end datasets in seconds:

```bash
# Generate 1000 paired-end reads (small dataset)
wgsim -N 1000 -1 100 -2 100 reference.fa reads_1.fq reads_2.fq
```


---

### Summary:

| Dataset                          | Size          | Use Case               | Download Links                                                                                                                                                                                                                                                                                                                              |
| -------------------------------- | ------------- | ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| ENCODE ENCFF001TDO + ENCFF001TDP | \~2.7 MB each | Small test dataset     | [https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.fastq.gz](https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.fastq.gz)  <br> [https://www.encodeproject.org/files/ENCFF001TDP/@@download/ENCFF001TDP.fastq.gz](https://www.encodeproject.org/files/ENCFF001TDP/@@download/ENCFF001TDP.fastq.gz) |
| 1000 Genomes ERR031844           | \~10 MB       | Low coverage real data | ftp\://ftp.sra.ebi.ac.uk/vol1/fastq/ERR031/ERR031844/ERR031844\_1.fastq.gz  <br> ftp\://ftp.sra.ebi.ac.uk/vol1/fastq/ERR031/ERR031844/ERR031844\_2.fastq.gz                                                                                                                                                                                 |
| Simulated wgsim                  | User-defined  | Synthetic data         | `wgsim` command-line tool (part of samtools)                                                                                                                                                                                                                                                                                                |
