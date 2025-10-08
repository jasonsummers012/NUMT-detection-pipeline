NUMT Detection Pipeline
A comprehensive bioinformatics pipeline for identifying Nuclear Mitochondrial DNA segments (NUMTs) in human genomic data using paired-end sequencing reads.

Overview
NUMTs (Nuclear Mitochondrial DNA segments) are fragments of mitochondrial DNA that have been inserted into the nuclear genome. This pipeline detects NUMT insertions by identifying:

Discordant read pairs - where one read maps to nuclear DNA and its mate maps to mitochondrial DNA
Split reads - reads that span the nuclear-mitochondrial junction
Pipeline Workflow
BAM File → Filter Reads → Cluster Candidates → Extract Split Reads → Validate → Final NUMTs
Step-by-Step Process
BAM Parsing (01_bam_parsing_NUMT.ipynb)
Filters high-quality reads (MAPQ ≥ 30)
Identifies discordant nuclear-mitochondrial read pairs
Output: 1,917 candidate reads
Candidate Clustering (02_candidate_clustering.ipynb)
Groups reads into 500bp genomic bins
Filters clusters with ≥2 supporting reads
Output: 184 candidate clusters
Split-Read Identification (03_split_read_identification.ipynb)
Extracts soft-clipped sequences from cluster regions
Realigns clipped sequences to mitochondrial reference genome
Identifies reads spanning nuclear-MT junctions
Split-Read Counting (04_split_read_counting.ipynb)
Counts validated split reads per cluster
Filters clusters with ≥5 confirmed split reads
Output: 88 high-confidence NUMT candidates
Data Visualization (05_data_visualization.ipynb)
Generates comprehensive multi-panel figure
Shows genome-wide distribution, confidence metrics, and top candidates
Requirements
Software
Python 3.9+
BWA (Burrows-Wheeler Aligner)
SAMtools
Jupyter Notebook
Python Packages
bash
conda install pandas numpy matplotlib seaborn
conda install -c bioconda pysam
Installation
bash
# Clone the repository
git clone https://github.com/jasonsummers012/NUMT-detection-pipeline.git
cd NUMT-detection-pipeline

# Create conda environment
conda create -n bioinfo_env python=3.9
conda activate bioinfo_env

# Install dependencies
conda install pandas numpy matplotlib seaborn
conda install -c bioconda pysam bwa samtools
Usage
1. Prepare Reference Genome
Create a combined reference with both nuclear and mitochondrial genomes:

bash
cat hg38.fa human_mtDNA.fa > data/combined_reference.fa
2. Run Alignment
bash
bash run_alignment.sh
This script:

Indexes the reference genome
Aligns paired-end reads with BWA MEM
Converts to sorted BAM with index
3. Run NUMT Detection Pipeline
Execute notebooks in order:

bash
jupyter notebook
Open and run:

01_bam_parsing_NUMT.ipynb
02_candidate_clustering.ipynb
03_split_read_identification.ipynb
04_split_read_counting.ipynb
05_data_visualization.ipynb
Results
Summary Statistics
88 confirmed NUMT candidates identified
44 high-confidence NUMTs (≥10 split reads)
Top candidate: chr12:51,346,632 with 149 split reads
Output Files
File	Description
NUMT_candidates.csv	Initial discordant read pairs (1,917 reads)
filtered_NUMT_candidates.csv	Clustered candidates (184 clusters)
confirmed_NUMT_clusters.csv	Final validated NUMTs (88 clusters)
NUMT_comprehensive_analysis.png/pdf	Multi-panel visualization figure
Visualization
Show Image

The comprehensive figure includes:

Panel A: Genome-wide NUMT distribution across chromosomes
Panel B: Confidence scoring (split reads vs total reads)
Panel C: Mitochondrial origin map (which MT regions are inserted)
Panel D: Read support distribution
Panel E: Top 10 high-confidence NUMTs
Manual Validation
Top NUMT candidates were manually validated using IGV (Integrative Genomics Viewer):

Show Image

Visual inspection confirms:

Clear junction breakpoint with clustered soft-clipped reads
Discordant read pairs mapping to nuclear and mitochondrial genomes
High mapping quality and consistent read support
Directory Structure
NUMT-detection-pipeline/
├── data/                           # Reference genomes and input data
│   ├── combined_reference.fa
│   ├── human_mtDNA.fa
│   └── SRR13269374_1/2.fastq
├── results/                        # Pipeline outputs
│   ├── sample_alignment_sorted.bam
│   ├── NUMT_candidates.csv
│   ├── confirmed_NUMT_clusters.csv
│   └── NUMT_comprehensive_analysis.pdf
├── notebooks/                      # Jupyter notebooks
│   ├── 01_bam_parsing_NUMT.ipynb
│   ├── 02_candidate_clustering.ipynb
│   ├── 03_split_read_identification.ipynb
│   ├── 04_split_read_counting.ipynb
│   └── 05_data_visualization.ipynb
├── run_alignment.sh               # Alignment script
└── README.md
Key Parameters
Parameter	Value	Description
MAPQ_THRESHOLD	30	Minimum mapping quality
BIN_SIZE	500	Genomic bin size (bp)
MIN_NUM_READS	2	Minimum reads per cluster
MIN_CLIP_LEN	10	Minimum soft-clip length
MIN_SPLIT_READS	5	Minimum split reads for confirmation
Methods
Quality Filters
Mapping quality ≥ 30
No unmapped reads
No duplicates, QC failures, or secondary alignments
Both reads in pair must be mapped
NUMT Detection Strategy
Identify discordant pairs where one read maps to nuclear genome (chr1-22, X, Y) and mate maps to mitochondrial genome
Cluster discordant pairs into genomic bins (500bp windows)
Extract soft-clipped sequences from cluster regions
Realign clips to mitochondrial reference genome
Count validated split reads per cluster
Filter for high-confidence candidates (≥5 split reads)
Validation
NUMTs can be validated against known databases:

UCSC Genome Browser NumtS track
Literature-based NUMT annotations
Manual inspection with IGV
Our top candidate (chr12:51,346,632) shows no overlap with known NUMT databases, suggesting it may be a novel insertion or population-specific variant.

Future Directions
Compare against published NUMT databases for validation
Estimate insertion sizes and breakpoint coordinates
Identify gene disruptions or functional impacts
Population-level analysis across multiple samples
Citation
If you use this pipeline, please cite:

[Jason Summers]. (2025). NUMT Detection Pipeline. 
GitHub: https://github.com/jasonsummers012/NUMT-detection-pipeline
License
MIT License - feel free to use and modify for your research.

Contact
For questions or issues, please open an issue on GitHub or contact [your email].

Acknowledgments
Pipeline developed for genomic analysis of nuclear mitochondrial insertions
Uses standard bioinformatics tools: BWA, SAMtools, pysam
Visualization created with matplotlib and seaborn
