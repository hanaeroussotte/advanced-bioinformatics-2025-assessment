# advanced-bioinformatics-2025-assessment

## 🧬 Advanced Bioinformatics – Assignment Submission

This repository contains two Bash scripts designed to execute a standard NGS data analysis pipeline and an R Markdown file exploring other related bioinformatics exercises. It was developed as part of the **Advanced Bioinformatics 2025** assignment to demonstrate command-line workflow design and scripting.

---

## 🗂️ Contents

- `Hanae_Roussotte_AF52866_Adv.BioinformaticsAssignment_BWA.sh`: Main pipeline using **BWA-MEM**, with detailed comments and rationale.
- `Hanae_Roussotte_AF52866_Adv.BioinformaticsAssignment_Bowtie.sh`: Alternative pipeline using **Bowtie**, with a streamlined version of the same workflow.
- `Hanae_Roussotte_AF52866_Adv.BIoinformaticsAssignment_RFiles`: Folder containing the R Markdown document and its rendered HTML output.

The answer sheet (PDF) has been submitted separately via **KEATS**.

---

## ⚙️ Installation and Dependencies

The pipeline scripts use **Miniconda** and **Bioconda** to install all required NGS tools, ensuring reproducibility in a Linux environment. R and RStudio were used for the RNA-seq and ChIP-seq analysis, included as R Markdown (`.Rmd`) and rendered HTML output.

Example conda setup for the Bash pipeline:

```bash
conda create -y -n ngs_pipeline
conda activate ngs_pipeline

conda install -y -c bioconda fastqc trimmomatic bwa samtools picard freebayes bedtools bcftools annovar snpeff

