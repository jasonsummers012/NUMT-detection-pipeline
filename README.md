NUMT Detection Pipeline
A bioinformatics pipeline for detecting Nuclear Mitochondrial DNA segments (NUMTs) from paired-end sequencing data. This project identifies mitochondrial DNA sequences that have been integrated into the nuclear genome by analyzing discordant read pairs and split reads.

Overview

Nuclear Mitochondrial DNA segments (NUMTs) are sequences of mitochondrial origin that have been inserted into the nuclear genome during evolution. This pipeline uses a multi-step approach to identify NUMT candidates from BAM alignment files by:

1. Initial filtering - Identifying high-quality read pairs where one read maps to nuclear DNA and its mate maps to mitochondrial DNA
2. Clustering - Grouping nearby NUMT candidates into clusters to reduce noise
3. Split-read analysis - Extracting and analyzing soft-clipped sequences to confirm NUMT boundaries
4. Validation - Confirming NUMT candidates by aligning split reads back to the mitochondrial reference

Requirements

Software Dependencies

Python 3.9+ with the following packages:
pysam - BAM/SAM file manipulation
pandas - Data analysis and manipulation
matplotlib - Plotting and visualization
seaborn - Statistical data visualization
bwa - Burrows-Wheeler Aligner for sequence alignment
samtools - SAM/BAM file utilities

Data Requirements

Paired-end FASTQ files (e.g., SRR13269374_1.fastq, SRR13269374_2.fastq)
Combined reference genome (nuclear + mitochondrial)
Human mitochondrial DNA reference (human_mtDNA.fa)

Installation

1. Clone this repository:
bash
git clone <github.com/jasonsummers012/NUMT-detection-pipeline>
cd numt-detection-pipeline
2. Create and activate a conda environment:
bash
conda create -n bioinfo_env python=3.9
conda activate bioinfo_env
3. Install Python dependencies:
bash
pip install pysam pandas matplotlib seaborn
4. Install bioinformatics tools:
bash
conda install -c bioconda bwa
conda install -c bioconda samtools

Directory Structure

project/
├── data/
│   ├── SRR13269374_1.fastq          # Forward reads
│   ├── SRR13269374_2.fastq          # Reverse reads
│   ├── combined_reference.fa         # Nuclear + mitochondrial reference
│   └── human_mtDNA.fa               # Mitochondrial reference only
├── results/                         # Output directory (created automatically)
├── notebooks/
│   ├── 01_bam_parsing_NUMT.ipynb
│   ├── 02_candidate_clustering.ipynb
│   ├── 03_split_read_identification.ipynb
│   └── 04_split_read_counting.ipynb
└── run_alignment.sh                 # Initial alignment script

Usage

Step 1: Initial Alignment
Run the alignment script to generate sorted BAM files:

bash
chmod +x run_alignment.sh
./run_alignment.sh
This script will:

Index the reference genome (if not already done)
Align paired-end reads using BWA MEM
Convert SAM to BAM format
Sort and index the BAM file

Step 2: NUMT Detection Pipeline
Execute the Jupyter notebooks in order:

1. BAM Parsing (01_bam_parsing_NUMT.ipynb)
Filters reads based on mapping quality (MAPQ ≥ 30)
Identifies discordant read pairs (nuclear-mitochondrial pairs)
Exports NUMT candidates to CSV

2. Candidate Clustering (02_candidate_clustering.ipynb)
Groups nearby NUMT candidates into 500bp bins
Filters clusters requiring minimum 2 reads per cluster
Visualizes cluster distribution across chromosomes
Exports filtered clusters to CSV

3. Split Read Identification (03_split_read_identification.ipynb)
Extracts soft-clipped sequences from reads in NUMT regions
Generates FASTQ file of clipped sequences
Aligns clipped sequences to mitochondrial reference
Creates sorted BAM of clipped vs mtDNA alignments

4. Split Read Counting (04_split_read_counting.ipynb)
Counts split reads that align to mitochondrial DNA for each cluster
Applies final filtering (minimum 5 split reads per cluster)
Exports confirmed NUMT clusters

Parameters

Key parameters that can be adjusted:

MAPQ_THRESHOLD = 30 - Minimum mapping quality
BIN_SIZE = 500 - Cluster bin size (bp)
MIN_NUM_READS = 2 - Minimum reads per cluster
MIN_CLIP_LEN = 10 - Minimum clipped sequence length
MIN_SPLIT_READS = 5 - Minimum split reads for confirmation

Output Files

The pipeline generates several intermediate and final output files:

NUMT_candidates.csv - Initial NUMT candidates
filtered_NUMT_candidates.csv - Clustered and filtered candidates
clipped_reads.fq - FASTQ of soft-clipped sequences
clipped_reads_with_cluster.csv - Clipped reads with cluster assignments
clipped_vs_mtDNA_sorted.bam - Alignment of clipped reads to mtDNA
confirmed_NUMT_clusters.csv - Final confirmed NUMT locations

Methodology

This pipeline implements the approach described in "The Mighty NUMT: Mitochondrial DNA Flexing its Code in the Nuclear Genome" and uses:

1. Discordant read pair analysis - Identifies read pairs where one mate maps to nuclear DNA and the other to mitochondrial DNA
2. Spatial clustering - Groups nearby events to distinguish real NUMTs from random noise
3. Split-read validation - Analyzes soft-clipped sequences to confirm NUMT boundaries and validate insertions

Citation

If you use this pipeline in your research, please cite the original methodology paper:
"The Mighty NUMT: Mitochondrial DNA Flexing its Code in the Nuclear Genome"

Contact

Jason Summers
jasonsummers012@gmail.com

Troubleshooting

Common Issues
Memory errors: Large BAM files may require more memory. Consider processing chromosomes individually.
Missing index files: Ensure all BAM files are properly indexed using samtools index.
Reference format: Verify that chromosome naming is consistent between reference and BAM files.
Path issues: Update file paths in notebooks if using different directory structure.

Performance Tips

Use multiple threads for BWA alignment (-t parameter)
Index BAM files for faster random access
Consider filtering by chromosome to reduce memory usage

Version History

v1.0 - Initial release with complete NUMT detection pipeline
