#!/bin/bash
# NGS Pipeline using Bowtie2 instead of BWA for alignment
# All other steps remain the same as the original pipeline

# Create and activate conda environment
conda create -y -n ngs_pipeline_alt
conda activate ngs_pipeline_alt

# Install tools - same as original except Bowtie2 instead of BWA
conda install -y -c bioconda fastqc=0.11.9
conda install -y -c bioconda trimmomatic=0.39
conda install -y -c bioconda bowtie2=2.4.5  # Alternative aligner
conda install -y -c bioconda samtools=1.13
conda install -y -c bioconda picard=2.26.0
conda install -y -c bioconda freebayes=1.3.5
conda install -y -c bioconda bedtools=2.30.0
conda install -y -c bioconda bcftools=1.13
conda install -y -c bioconda annovar
conda install -y -c bioconda snpeff=5.0

# Create directory structure - same as original
mkdir -p raw_data reference_genome results/{fastqc,trimmed,alignment,variants,annotation}

# Download input files - same as original
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz -O raw_data/NGS0001.R1.fastq.gz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz -O raw_data/NGS0001.R2.fastq.gz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed -O raw_data/annotation.bed
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz -O reference_genome/hg19.fa.gz
gunzip reference_genome/hg19.fa.gz

# Index reference genome using Bowtie2 instead of BWA
# Bowtie2 creates different index files (.bt2 extension)
bowtie2-build reference_genome/hg19.fa reference_genome/hg19

# Create samtools index - same as original
samtools faidx reference_genome/hg19.fa

# Pre-alignment QC - same as original
fastqc raw_data/NGS0001.R1.fastq.gz raw_data/NGS0001.R2.fastq.gz -o results/fastqc -t 4

trimmomatic PE -threads 4 \
  raw_data/NGS0001.R1.fastq.gz raw_data/NGS0001.R2.fastq.gz \
  results/trimmed/NGS0001.R1.trimmed.fastq.gz results/trimmed/NGS0001.R1.unpaired.fastq.gz \
  results/trimmed/NGS0001.R2.trimmed.fastq.gz results/trimmed/NGS0001.R2.unpaired.fastq.gz \
  ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

fastqc results/trimmed/NGS0001.R1.trimmed.fastq.gz results/trimmed/NGS0001.R2.trimmed.fastq.gz -o results/fastqc -t 4

# Alignment using Bowtie2 instead of BWA
SAMPLE="NGS0001"
PLATFORM="ILLUMINA"
LIBRARY="Lib1"

# Bowtie2 alignment with read group information
# This is the key difference from the original pipeline
bowtie2 -p 4 --sensitive-local \
  --rg-id "${SAMPLE}" --rg "SM:${SAMPLE}" --rg "PL:${PLATFORM}" --rg "LB:${LIBRARY}" \
  -x reference_genome/hg19 \
  -1 results/trimmed/NGS0001.R1.trimmed.fastq.gz \
  -2 results/trimmed/NGS0001.R2.trimmed.fastq.gz | \
  samtools view -bS - > results/alignment/NGS0001.bam

# Sort BAM file by coordinate position
samtools sort -o results/alignment/NGS0001.sorted.bam results/alignment/NGS0001.bam

# Index the sorted BAM file for random access
samtools index results/alignment/NGS0001.sorted.bam

# Mark duplicates using Picard
picard MarkDuplicates \
  I=results/alignment/NGS0001.sorted.bam \
  O=results/alignment/NGS0001.marked.bam \
  M=results/alignment/NGS0001.metrics.txt \
  CREATE_INDEX=true

# Quality filter BAM file to remove low-quality and problematic alignments
samtools view -b -F 0x704 -q 20 results/alignment/NGS0001.marked.bam > results/alignment/NGS0001.filtered.bam

# Index the filtered BAM file
samtools index results/alignment/NGS0001.filtered.bam

# Generate flagstat - summary statistics of flags in the BAM file
samtools flagstat results/alignment/NGS0001.filtered.bam > results/alignment/NGS0001.flagstat.txt

# Generate idxstats - chromosome-level alignment statistics
samtools idxstats results/alignment/NGS0001.filtered.bam > results/alignment/NGS0001.idxstats.txt

# Calculate depth of coverage
samtools depth -a results/alignment/NGS0001.filtered.bam | \
  awk '{sum+=$3} END {print "Average coverage: "sum/NR}' > results/alignment/NGS0001.coverage.txt

# Analyse insert size metrics with Picard
picard CollectInsertSizeMetrics \
  I=results/alignment/NGS0001.filtered.bam \
  O=results/alignment/NGS0001.insert_size_metrics.txt \
  H=results/alignment/NGS0001.insert_size_histogram.pdf


# Call variants with Freebayes
freebayes -f reference_genome/hg19.fa \
  --targets raw_data/annotation.bed \
  results/alignment/NGS0001.filtered.bam > results/variants/NGS0001.vcf

# Filter variants based on quality and depth using bcftools
bcftools filter -i 'QUAL>20 && DP>10' \
  results/variants/NGS0001.vcf > results/variants/NGS0001.filtered.vcf


# 1.Annotate variants using ANNOVAR
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


# 2. Annotate variants using SnpEff
java -Xmx4g -jar $CONDA_PREFIX/share/snpeff/snpEff.jar \
  -v GRCh37.75 \
  results/variants/NGS0001.filtered.vcf > results/annotation/NGS0001.snpeff.vcf

# Prioritise variants: select exonic variants not seen in dbSNP
grep "exonic" results/annotation/NGS0001.annovar.hg19_multianno.txt | \
  grep -v "rs[0-9]" > results/annotation/NGS0001.prioritized.txt

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
