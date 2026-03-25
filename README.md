# PA1874 Mutation Analysis and Visualisation Workflow

This repository contains the R workflow used to explore sequence variation around **PA1874** in a collection of *Pseudomonas aeruginosa* PES genomes.

The script brings together mutation calls from multiple sources, compares variant callers, maps SNP positions to relative and structural coordinates, and generates visual summaries focused on the PA1874 region.

## Project overview

This workflow was developed as part of a bioinformatics dissertation investigating the genetic and structural context of **PA1874**, a gene of interest associated with transmissibility in an epidemic PES strain of *Pseudomonas aeruginosa*.

The analysis focuses on:

- comparing mutation calls from **breseq**, **marginal breseq output**, and **snippy**
- examining mutation density across the full genome and within the **PA1874** region
- converting genomic SNP positions into **relative PA1874 coordinates**
- linking relative positions to **AlphaFold residue positions**
- summarising SNP appearances across genomes, years, and patients
- visualising docking results and RTX motif scan results relevant to PA1874 interpretation

## Main components of the workflow

### 1. Mutation density plots
The script generates sliding-window mutation density plots for:
- the full genome
- zoomed-in regions around PA1874
- selected genomes and selected timepoints

These plots were used to visually compare mutation patterns across variant callers and across isolates.

### 2. Relative SNP position mapping
Marginal SNP positions are extracted from HTML output files and filtered to the genomic coordinates of PA1874. These positions are then converted into relative positions within the gene, taking into account that PA1874 lies on the reverse strand.

### 3. AlphaFold position mapping
Relative PA1874 positions are linked to AlphaFold residue positions using a combined lookup table. This allows sequence variation to be interpreted in a structural context.

### 4. Breseq vs snippy comparison
The workflow includes both combined and per-genome comparisons of SNP calls from breseq and snippy, allowing agreement and disagreement between callers to be visualised directly.

### 5. Docking result visualisation
A final docking scatter plot is included to compare:
- largest cluster size
- best energy score

Thresholds based on negative controls are added to support interpretation of docking convergence and interaction strength.

### 6. Multi-patient SNP visualisation
The script builds a combined dataset of PA1874 SNP appearances across patients and years, and generates a multi-patient plot showing:
- genomic position
- collection year
- caller type
- patient identity

### 7. RTX motif scan summary
A summary bar plot is included to compare the number of RTX motifs detected across known adhesins and PA1874.

## Input files expected

The workflow expects files such as:

- `breseq_*.gd`
- `marginal_*.html`
- `snippy_*.csv`
- `snps_*.html`
- CSV mapping files for relative-to-AlphaFold coordinate conversion

Some parts of the workflow also assume specific folder structures such as:

- `marginal_htmls/`
- `marginal_relative_tables/`
- `marginal_with_alphafold/`
- genome-specific subfolders such as `A005/`, `A138/`, etc.

## Output generated

The script generates outputs including:

- mutation density plots
- caller comparison plots
- relative SNP coordinate tables
- AlphaFold-mapped SNP tables
- docking plots
- multi-patient PA1874 SNP plots
- RTX motif scan plot

## Packages used

Main R packages used in this workflow include:

- `dplyr`
- `readr`
- `stringr`
- `ggplot2`
- `ggrepel`
- `purrr`

## Author

Terenia Lee Leh Yi  
MSc Bioinformatics and Systems Biology  
University of Manchester
