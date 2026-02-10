# iMAP Part 09 — Statistical Analysis

This repository contains the **statistical analysis and inference stage** of the iMAP (Integrated Microbiome Analysis Pipeline) framework.

Part 09 focuses on formal statistical testing of processed 16S-based microbiome data, building directly on the exploratory insights established in Part 08.
The emphasis is on **assessing evidence, quantifying uncertainty, and supporting biologically meaningful conclusions**, rather than applying tests in isolation.

This component represents the transition from exploratory interpretation to statistical decision-making and downstream claims, and is intended to be used
after data structure, variability, and experimental context have been examined.

For pipeline structure and upstream context, see the **iMAP Project Overview** and the **Exploratory Analysis (Part 08)** component.

## Related iMAP Repositories

Below is a summary of the iMAP repositories, showing how this statistical analysis stage fits within the full pipeline. Each repository represents a distinct step in the workflow;
together they document the complete analytical structure of iMAP.

<br>

| IMAP-Repo                                      | Description                                               | GH-Pages                                             |
| :---------------------------------------------- | :--------------------------------------------------------- | :-----------------------------------------------------: |
| [OVERVIEW](https://github.com/tmbuza/imap-project-overview/) | IMAP project overview                                     | [Link](https://tmbuza.github.io/imap-project-overview/) |
| [PART 01](https://github.com/tmbuza/imap-essential-software/) | Software requirements for microbiome data analysis with Snakemake workflows | [Link](https://tmbuza.github.io/imap-essential-software/) |
| [PART 02](https://github.com/tmbuza/imap-sample-metadata/) | Downloading and exploring microbiome sample metadata from SRA Database | [Link](https://tmbuza.github.io/imap-sample-metadata/) |
| [PART 03](https://github.com/tmbuza/imap-download-sra-reads/) | Downloading and filtering microbiome sequencing data from SRA database | [Link](https://tmbuza.github.io/imap-download-sra-reads/) |
| [PART 04](https://github.com/tmbuza/imap-read-quality-control/) | Quality control of microbiome next-generation sequencing reads | [Link](https://tmbuza.github.io/imap-read-quality-control/) |
| [PART 05](https://github.com/tmbuza/imap-bioinformatics-mothur/) | Microbial profiling using MOTHUR and Snakemake workflows | [Link](https://tmbuza.github.io/imap-mothur-bioinformatics/) |
| [PART 06](https://github.com/tmbuza/imap-bioinformatics-qiime2/) | Microbial profiling using QIIME2 and Snakemake workflows | [Link](https://tmbuza.github.io/imap-qiime2-bioinformatics/) |
| [PART 07](https://github.com/tmbuza/imap-data-processing/) | Processing output from 16S-based microbiome bioinformatics pipelines | [Link](https://tmbuza.github.io/imap-data-preparation/) |
| [PART 08](https://github.com/tmbuza/imap-exploratory-analysis/) | Exploratory analysis of processed 16S-based microbiome data | [Link](https://tmbuza.github.io/imap-data-exploration/) |
| [PART 09](https://github.com/tmbuza/imap-statistical-analysis/) | Statistical analysis of processed 16S-based microbiome data | [Link](https://tmbuza.github.io/imap-statistical-analysis/) |
| [PART 10](https://github.com/tmbuza/imap-machine-learning/) | Machine learning analysis of processed 16S-based microbiome data | [Link](https://tmbuza.github.io/imap-machine-learning/) |


## Session information

For a detailed overview of the tools and versions suitable for this guide, explore the [session information](session_info.txt).

## Citation
> Please consider citing the [iMAP article](https://rdcu.be/b5iVj) if you find any part of the iMAP practical user guides helpful in your microbiome data analysis.

Buza, T. M., Tonui, T., Stomeo, F., Tiambo, C., Katani, R., Schilling, M., … Kapur, V. (2019). iMAP: An integrated bioinformatics and visualization pipeline for microbiome data analysis. BMC Bioinformatics, 20. https://doi.org/10.1186/S12859-019-2965-4

## :tada: Raise awareness
> Please help increase awareness of freely available tools for microbiome data analysis.
See [Dimensions of the iMAP article](https://badge.dimensions.ai/details/id/pub.1117740326)
