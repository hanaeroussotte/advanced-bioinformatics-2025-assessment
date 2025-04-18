#!/bin/bash

# Create and activate conda environment
# Using conda enables isolated environment management for bioinformatics tools
# Miniconda is recommended due to its minimal footprint
conda create -y -n ngs_pipeline
conda activate ngs_pipeline

# Install necessary bioinformatics tools
# FastQC - Quality control tool for high throughput sequence data
# Used for assessing sequence quality before and after trimming
conda install -y -c bioconda fastqc=0.11.9

# Trimmomatic - Flexible trimming tool for Illumina NGS data
# Removes adapters, low quality bases, and filters short reads
conda install -y -c bioconda trimmomatic=0.39

# BWA - Burrows-Wheeler Aligner for mapping sequence reads to reference genomes
# Specifically optimized for Illumina reads, with BWA-MEM algorithm for longer sequences
conda install -y -c bioconda bwa=0.7.17

# Samtools - Suite of utilities for manipulating alignments in SAM/BAM format
# Essential for processing, sorting, indexing, and analyzing alignment files
conda install -y -c bioconda samtools=1.13

# Picard - Tools for manipulating high-throughput sequencing data and formats
# Required for marking duplicate reads resulting from PCR amplification
conda install -y -c bioconda picard=2.26.0

# Freebayes - Bayesian haplotype-based genetic variant detector
# Powerful variant caller able to detect SNPs, indels, and complex variants
conda install -y -c bioconda freebayes=1.3.5

# BEDTools - Tools for genomic arithmetic and interval manipulations
# Necessary for operations on genomic intervals in BED files
conda install -y -c bioconda bedtools=2.30.0

# BCFTools - Utilities for variant calling and manipulating VCFs and BCFs
# Used for filtering and analyzing variant call files
conda install -y -c bioconda bcftools=1.13

# ANNOVAR - Efficient software tool for functional annotation of genetic variants
# Comprehensive annotation of variants from multiple databases
conda install -y -c bioconda annovar

# SnpEff - Genomic variant annotations and functional effect prediction
# Provides detailed variant classification and potential functional impacts
conda install -y -c bioconda snpeff=5.0


# Create project directories with organised structure
# This facilitates workflow management and output organization
mkdir -p raw_data reference_genome results/{fastqc,trimmed,alignment,variants,annotation}

# Download and rename FastQ files
# The .qz extension is non-standard, so we rename to standard .gz
# -O flag specifies the output filename with proper extension
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz -O raw_data/NGS0001.R1.fastq.gz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz -O raw_data/NGS0001.R2.fastq.gz

# Download annotation BED file containing regions of interest
# BED format contains genomic coordinates to focus variant calling 
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed -O raw_data/annotation.bed

# Download and prepare reference genome hg19 (human genome version 19)
# Reference genome is required for alignment and variant calling
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz -O reference_genome/hg19.fa.gz
gunzip reference_genome/hg19.fa.gz

# Index reference genome for BWA alignment
# Indexing creates auxiliary files needed for efficient sequence alignment
# This is a computationally intensive one-time preprocessing step
cd reference_genome
bwa index hg19.fa     # Creates .amb, .ann, .bwt, .pac, and .sa index files
samtools faidx hg19.fa  # Creates .fai index for random access
cd . .

# Run FastQC on raw reads to assess quality
# -o specifies output directory for HTML reports
# -t 4 uses 4 threads for parallel processing to speed up analysis
fastqc raw_data/NGS0001.R1.fastq.gz raw_data/NGS0001.R2.fastq.gz -o results/fastqc -t 4

# Trim reads using Trimmomatic
# PE mode for paired-end data ensures proper handling of both forward and reverse reads
# ILLUMINACLIP removes adapter sequences:
#   - TruSeq3-PE.fa contains Illumina adapter sequences
#   - 2:30:10 parameters mean:
#     - 2: Maximum mismatch count for seed alignment
#     - 30: Palindrome clip threshold (for adapter read-through)
#     - 10: Simple clip threshold
# LEADING/TRAILING: Trim low quality bases (below quality 3) from start/end
# SLIDINGWINDOW: Scan with 4-base window, cut when average quality falls below 15
# MINLEN: Drop reads below 36 bases long after trimming to avoid ultrashort reads
trimmomatic PE -threads 4 \
  raw_data/NGS0001.R1.fastq.gz raw_data/NGS0001.R2.fastq.gz \
  results/trimmed/NGS0001.R1.trimmed.fastq.gz results/trimmed/NGS0001.R1.unpaired.fastq.gz \
  results/trimmed/NGS0001.R2.trimmed.fastq.gz results/trimmed/NGS0001.R2.unpaired.fastq.gz \
  ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


# Run FastQC on trimmed reads to verify quality improvement
# This allows comparison of quality metrics before and after trimming
# Same parameters as before, examining the trimmed files
fastqc results/trimmed/NGS0001.R1.trimmed.fastq.gz results/trimmed/NGS0001.R2.trimmed.fastq.gz -o results/fastqc -t 4

# The key quality metrics to examine in FastQC reports include:
# - Per base sequence quality: Should show improved quality scores after trimming
# - Adapter content: Should be eliminated or greatly reduced after trimming
# - Per base N content: Should be minimal
# - Sequence length distribution: Will show the effect of trimming
# - Overrepresented sequences: Adapter sequences should be removed

# Set variables for read group information
# Read groups are crucial for tracking sample information
SAMPLE="NGS0001"
PLATFORM="ILLUMINA"
LIBRARY="Lib1"

# Align with BWA MEM including read group information
# BWA-MEM algorithm is selected because:
# - It's optimised for 70-100bp+ Illumina reads
# - Provides better performance for longer reads compared to BWA-aln
# - Supports local and chimeric alignment for higher sensitivity
# - Has robust error correction for noisy reads

# Parameters explained:
# -t 4: Use 4 threads for parallel processing
# -M: Mark shorter split hits as secondary (for compatibility with Picard)
# -R: Add read group information which is critical for:
#     - Identifying which reads came from which sample/library/run
#     - Required by downstream tools like GATK
#     - Enables proper handling of samples in multi-sample analyses
#     - Tracking technical bias by sequencing run
bwa mem -t 4 -M \
  -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:${PLATFORM}\tLB:${LIBRARY}" \
  reference_genome/hg19.fa \
  results/trimmed/NGS0001.R1.trimmed.fastq.gz results/trimmed/NGS0001.R2.trimmed.fastq.gz | \
  samtools view -bS - > results/alignment/NGS0001.bam

# Sort BAM file by coordinate position
# Sorting is required for indexing and essential for many downstream tools
# Coordinate sorting enables efficient access to alignments in specific regions
samtools sort -o results/alignment/NGS0001.sorted.bam results/alignment/NGS0001.bam

# Index the sorted BAM file for random access
# Indexing creates a .bai file that enables fast retrieval of alignments
# This is required for viewing in genome browsers and for many tools
samtools index results/alignment/NGS0001.sorted.bam

# Mark duplicates using Picard
# Duplicate reads likely arise from PCR amplification during library preparation
# Marking instead of removing allows downstream tools to account for duplicates
# Duplicates can cause false confidence in variant allele frequency
# Parameters explained:
# I: Input BAM file
# O: Output BAM with marked duplicates
# M: File to write duplication metrics
# CREATE_INDEX: Automatically create index for the output file
picard MarkDuplicates \
  I=results/alignment/NGS0001.sorted.bam \
  O=results/alignment/NGS0001.marked.bam \
  M=results/alignment/NGS0001.metrics.txt \
  CREATE_INDEX=true

# Quality filter BAM file to remove low-quality and problematic alignments
# Filtering improves the reliability of downstream variant calling
# Parameters explained:
# -b: Output in BAM format
# -F 0x704: Filter flag combining:
#    - 0x4: Unmapped reads (not aligned)
#    - 0x100: Secondary alignments (non-primary)
#    - 0x400: PCR/optical duplicates (marked by Picard)
# -q 20: Minimum mapping quality of 20 (99% probability alignment is correct)
samtools view -b -F 0x704 -q 20 results/alignment/NGS0001.marked.bam > results/alignment/NGS0001.filtered.bam

# Index the filtered BAM file
samtools index results/alignment/NGS0001.filtered.bam

Generate standard alignment statistics (i.e. flagstats, idxstats, depth of coverage, insert size) (4pts)
# Generate flagstat - summary statistics of flags in the BAM file
# Shows counts of read categories (mapped, paired, duplicates, etc.)
# Useful for quick assessment of alignment quality
samtools flagstat results/alignment/NGS0001.filtered.bam > results/alignment/NGS0001.flagstat.txt

# Generate idxstats - chromosome-level alignment statistics
# Reports number of mapped/unmapped reads per chromosome
# Helps identify chromosomal biases or mapping issues
samtools idxstats results/alignment/NGS0001.filtered.bam > results/alignment/NGS0001.idxstats.txt

# Calculate depth of coverage
# Important for assessing sequencing depth across the genome
# Average coverage affects confidence in variant calling
samtools depth -a results/alignment/NGS0001.filtered.bam | \
  awk '{sum+=$3} END {print "Average coverage: "sum/NR}' > results/alignment/NGS0001.coverage.txt

# Analyse insert size metrics with Picard
# Insert size is the distance between paired-end reads
# Distribution of insert sizes is important for:
#  - QC of library preparation
#  - Detection of structural variants
#  - Configuration of variant callers
picard CollectInsertSizeMetrics \
  I=results/alignment/NGS0001.filtered.bam \
  O=results/alignment/NGS0001.insert_size_metrics.txt \
  H=results/alignment/NGS0001.insert_size_histogram.pdf

# Call variants with Freebayes
# Freebayes is selected because:
# - It's a Bayesian haplotype-based variant detector
# - Can detect SNPs, indels, MNPs, and complex events
# - Models samples as a mixture of haplotypes
# - Doesn't require preprocessing steps like base quality recalibration

# Parameters explained:
# -f: Reference genome in FASTA format
# --targets: BED file restricting analysis to specific regions
#   - Using a BED file greatly improves efficiency
#   - Focuses analysis on regions of interest only
#   - Reduces computational requirements
freebayes -f reference_genome/hg19.fa \
  --targets raw_data/annotation.bed \
  results/alignment/NGS0001.filtered.bam > results/variants/NGS0001.vcf

# Filter variants based on quality and depth using bcftools
# Filtering reduces false positives while retaining true variants
# Parameters explained:
# -i: Include variants matching the filter expression
# 'QUAL>20': Only keep variants with quality score >20 (Phred scale)
#   - QUAL of 20 means 99% probability the variant is real
# 'DP>10': Only keep variants with read depth >10
#   - Minimum depth ensures sufficient evidence for the variant
bcftools filter -i 'QUAL>20 && DP>10' \
  results/variants/NGS0001.vcf > results/variants/NGS0001.filtered.vcf

# These filtering criteria represent a balance between:
# - Sensitivity (ability to detect true variants)
# - Specificity (ability to reject false positives)
# More stringent filters could be applied for clinical applications


# 1.Annotate variants using ANNOVAR
# ANNOVAR provides comprehensive annotation from multiple databases
# Essential for understanding the functional impact of variants

# Parameters explained:
# -buildver hg19: Specify reference genome build
# -out: Output file prefix
# -remove: Remove temporary files
# -protocol: Specify annotation databases to use
#   - refGene: Gene-based annotation (exonic, intronic, etc.)
#   - exac03: Population frequency from ExAC database
#   - dbsnp138: dbSNP identifiers and known variants
# -operation: Type of operation for each protocol
#   - g: Gene-based annotation
#   - f: Filter-based annotation
# -nastring: String to use for missing values
# -vcfinput: Input is in VCF format
perl $CONDA_PREFIX/bin/table_annovar.pl \
  results/variants/NGS0001.filtered.vcf \
  $CONDA_PREFIX/share/annovar/humandb/ \
  -buildver hg19 \
  -out results/annotation/NGS0001.annovar \
  -remove \
  -protocol refGene,exac03,dbsnp138 \
  -operation g,f,f \
  -nastring . \
  -vcfinput

# ANNOVAR would be chosen because:
# - It provides extensive functional annotation of genetic variants
# - Supports multiple gene definition systems and databases
# - Can identify disease-associated variants
# - Annotates variants with population frequencies
# - Provides effect predictions (synonymous, nonsynonymous, etc.)

# 2. Annotate variants using SnpEff
# SnpEff provides detailed variant effect and impact predictions
# Complements ANNOVAR with additional annotations and predictions

# Parameters explained:
# -Xmx4g: Allow Java to use up to 4GB of memory
# -v: Verbose output for monitoring progress
# GRCh37.75: Genome version compatible with hg19
#   - Includes ENSEMBL gene models version 75
java -Xmx4g -jar $CONDA_PREFIX/share/snpeff/snpEff.jar \
  -v GRCh37.75 \
  results/variants/NGS0001.filtered.vcf > results/annotation/NGS0001.snpeff.vcf

# SnpEff would be chosen because:
# - It provides detailed functional impact predictions
# - Classifies variants by effect (HIGH, MODERATE, LOW, MODIFIER)
# - Annotates using standardised terms (SO - Sequence Ontology)
# - Can identify loss-of-function and gain-of-function variants
# - Adds annotations directly to the VCF file
# - Generates HTML summary reports of variant effects

# Prioritise variants: select exonic variants not seen in dbSNP
# This filtering strategy focuses on novel coding variants
# Such variants are more likely to be functionally significant

# Steps explained:
# 1. grep "exonic": Select only variants in exonic regions
#    - Exonic variants directly affect protein sequence
#    - More likely to have functional impact than non-coding variants
# 2. grep -v "rs[0-9]": Exclude variants with dbSNP IDs
#    - Variants without rs IDs are novel or rare
#    - Novel variants are more likely to be specific to this sample
grep "exonic" results/annotation/NGS0001.annovar.hg19_multianno.txt | \
  grep -v "rs[0-9]" > results/annotation/NGS0001.prioritized.txt

# This prioritisation approach is valuable for:
# - Discovery of potentially novel disease-causing variants
# - Reducing the number of variants to manual review
# - Focusing on variants with higher likelihood of functional impact
# In a clinical setting, additional filters might be applied based on:
# - Minor allele frequency
# - Predicted pathogenicity scores
# - Inheritance patterns

# OPTIONAL : Generate a summary report of the pipeline results
echo "==== NGS Pipeline Report ====" > results/pipeline_summary.txt
echo "Sample: NGS0001" >> results/pipeline_summary.txt
echo "" >> results/pipeline_summary.txt

echo "=== Alignment Statistics ===" >> results/pipeline_summary.txt
grep "mapped (" results/alignment/NGS0001.flagstat.txt >> results/pipeline_summary.txt
cat results/alignment/NGS0001.coverage.txt >> results/pipeline_summary.txt
echo "" >> results/pipeline_summary.txt

echo "=== Variant Statistics ===" >> results/pipeline_summary.txt
echo "Total variants: $(grep -v "#" results/variants/NGS0001.vcf | wc -l)" >> results/pipeline_summary.txt
echo "Filtered variants: $(grep -v "#" results/variants/NGS0001.filtered.vcf | wc -l)" >> results/pipeline_summary.txt
echo "Exonic variants: $(grep "exonic" results/annotation/NGS0001.annovar.hg19_multianno.txt | wc -l)" >> results/pipeline_summary.txt
echo "Prioritized variants: $(wc -l < results/annotation/NGS0001.prioritized.txt)" >> results/pipeline_summary.txt

echo "Pipeline completed at $(date)" >> results/pipeline_summary.txt

# OPTIONAL : Clean up intermediate files to save disk space
rm results/alignment/NGS0001.bam
rm results/alignment/NGS0001.sorted.bam
rm results/alignment/NGS0001.sorted.bam.bai

# Deactivate conda environment + finish comment 
# Helps to know when the pipeline has finished running 
conda deactivate
echo "Pipeline completed successfully!”
