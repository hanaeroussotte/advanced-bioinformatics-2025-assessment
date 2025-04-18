# advanced-bioinformatics-2023-assessment

# NGS Variant Calling Pipeline – Assignment Submission

This repository contains a fully automated Bash script designed to execute a standard NGS data analysis pipeline from raw Illumina paired-end reads to variant calling and annotation. This script was developed as part of an academic assignment for evaluating bioinformatics skills in a Linux environment.

---

## Assignment Context

This project was completed in response to the following assessment brief:

> "You are required to share a Bash script that runs the workflow and takes the provided sequencing data as input. The script will be evaluated for its ability to install tools, perform read alignment, variant discovery, and annotation using a standard command-line NGS pipeline on a Linux OpenStack instance."

- Single-sample analysis using paired-end reads from a HiSeq 2500 run
- All steps automated via Bash
- Script tested on an Ubuntu OpenStack virtual machine, met some issues with storage even using the virtual machine

---

## 1. Installation and Dependencies

The script uses **Miniconda** to create an isolated environment and installs all necessary tools from **Bioconda**:

```bash
# Create and activate conda environment
conda create -y -n ngs_pipeline
conda activate ngs_pipeline

# Install required tools
conda install -y -c bioconda fastqc=0.11.9 trimmomatic=0.39 bwa=0.7.17 samtools=1.13 \
    picard=2.26.0 freebayes=1.3.5 bedtools=2.30.0 bcftools=1.13 annovar snpeff=5.0
