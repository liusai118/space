# Cosmic Connections: Spaceflight Redefines Ageing-Associated Microbiota

## Metagenomics Analysis

- **metagenomics_all**: The process of analyzing the metagenomics data of both the ageing and spaceflight cohorts.
  - These scripts contain detailed parameters for the upstream analysis of metagenomics:
    - For the **spaceflight cohort**:
      ```bash
      space_oral.sh
      space_gut.sh
      space_skin.sh
      ```
    - For the **ageing cohort**:
      ```bash
      skin_metagenomics.sh
      gut_metagenomics.sh
      oral_metagenomics.sh
      ```
  - All `.txt` files contain sample IDs of the ageing and spaceflight cohorts.

## Results

### Figure 1:

- **Analysis**:
  - `scRNA.ipynb`: Used to analyze the raw single-cell RNA sequencing data.
  - `H5ad_to_seurat.R`: Used to convert the scRNA data and perform differential expression gene (DEG) analysis.

### 

- **Data**:
  - `MR_all.csv`: Contains all the results of the Mendelian Randomization (MR) analysis.
  - Detailed parameters for MR analysis are provided in the scripts.

###

- **Species Annotation**:
  - `taxonmy_prepare.R`: Used to prepare the species annotation information table for basic analysis.

### Figure 2

- **Species Annotation**:
  - `taxonmy_prepare.R`: Used to prepare the species annotation information table for basic analysis.

### 

- **Data and Analysis**:
  - The `WGCNA` folder contains the scripts and data for Weighted Gene Co-expression Network Analysis (WGCNA).
  - Data for machine learning were prepared by the WGCNA analysis.
  - The `ML` folder contains all machine learning scripts.

### 

- **Analysis**:
  - Upstream Metatranscriptome Analysis:
    ```bash
    Skintrans.sh
    oraltrans.sh
    ```
  - Downstream Analysis:
    - The `R` folder contains all scripts for downstream metatranscriptome analysis.
    - `space_all.csv`: Metadata for metatranscriptome samples.
    - `species.tsv` and `all.tsv`: Operational Taxonomic Units (OTUs) for the metagenomics data of the spaceflight cohort.
    - `all_species_metatranscript.tsv` and `all_all_transcript.tsv`: OTUs for the metatranscriptome data of the spaceflight cohort.
    - `eggnog.KO.raw.txt`: Results of aligning metatranscriptome data to the eggNOG database.
  - Differential Gene Expression Analysis:
    ```bash
    Skin_FC.sh
    oral_FC.sh
    ```


