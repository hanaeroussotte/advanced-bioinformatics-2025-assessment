# advanced-bioinformatics-2023-assessment

# NGS Variant Calling Pipeline – Assignment Submission

This repository contains two fully automated Bash scripts designed to execute a standard NGS data analysis pipeline from raw Illumina paired-end reads to variant calling and annotation. These scripts were developed as part of an academic assignment for evaluating bioinformatics skills in a Linux environment.

---

## Assignment Context

This project was completed in response to the following assessment brief:

> "You are required to share a Bash script that runs the workflow and takes the provided sequencing data as input. The script will be evaluated for its ability to install tools, perform read alignment, variant discovery, and annotation using a standard command-line NGS pipeline on a Linux OpenStack instance."

- Single-sample analysis using paired-end reads from a HiSeq 2500 run  
- All steps automated via Bash  
- Scripts tested on an Ubuntu OpenStack virtual machine (note: encountered storage limitations)

---

## Contents

This repository includes **two Bash scripts**:

- `detailed-bashscript-pipeline-bwa.sh`: The main pipeline using **BWA-MEM** for alignment. This script includes **detailed explanations and comments** for each command and tool used, justifying parameter choices and highlighting good practices.
- `bashscript-pipeline-bowtie.sh`: An alternative version of the pipeline using **Bowtie** as the aligner. This script is streamlined and does **not include extensive command comments**, but demonstrates the same logic and structure using a different alignment strategy.

---

## Installation and Dependencies

Both scripts use **Miniconda** to create an isolated environment and install required tools from **Bioconda**. Example from the main pipeline:

```bash
# Create and activate conda environment
conda create -y -n ngs_pipeline
conda activate ngs_pipeline

# Install required tools
conda install -y -c bioconda fastqc=0.11.9 trimmomatic=0.39 bwa=0.7.17 samtools=1.13 \
    picard=2.26.0 freebayes=1.3.5 bedtools=2.30.0 bcftools=1.13 annovar snpeff=5.0
