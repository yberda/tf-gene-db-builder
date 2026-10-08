# TF Target Gene Databases — ReMap2022 × ENCODE rE2G/CREdb/GeneHancer

## Overview

This repository contains the R scripts used to build curated transcription
factor (TF) target gene databases by integrating two sources of regulatory
information: **ReMap2022** ChIP-seq peak data (~8,000 experiments across
hundreds of TFs) and **ENCODE rE2G**, **CREdb** or **GeneHancer** 
regulatory region–gene interaction predictions. Promoter regions (±1,000 bp 
around each TSS) are included in every general variant. For each TF, the
pipeline identifies which genes have a ReMap2022 peak overlapping a predicted
regulatory region (enhancer or promoter), then stores the result as a compact R
list: element `[[1]]` is a character vector of unique Entrez Gene IDs, and
element `[[2]]` is a named list where each entry is a ChIP-seq experiment ID
mapping to integer indices into `[[1]]`. Both general (pan-cell-type) and
cell-type-specific database variants are produced.

---

## Repository structure

```
tf-target-db/
├── R/
│   ├── 00_session_info.R           # Record package versions for reproducibility
│   ├── 01_build_remap.R            # Load ReMap2022 files → named GRanges list
│   ├── 02_build_promoters.R        # Define hg38 TSS ±1 kb promoter regions
│   ├── 03_filter_re2g_by_depth.R    # Filters rE2G ranges by recurrence
│   ├── 04_build_general_databases.R # General TF–gene databases 
│   └── 05_build_celltype_databases.R # Cell-type-specific TF–gene databases
├── data/                           # (not tracked by git)
│   ├── remap_hg38/                 # One .txt file per ReMap2022 experiment
│   ├── all_encode_re2g.txt         # Merged ENCODE rE2G table (pipe-delimited)
│   ├── ReMap2022+EnsTSS+GH.Rdata   # GeneHancer object from TFEA.ChIP
│   ├── generationDB_lps.R          # script to generate GeneHancer object
│   └── encode_re2g_celltype/       # Per-biosample ENCODE rE2G BED.gz files
├── out_rdata/                      # Generated outputs (not tracked by git)
│   ├── re2g_filtered_300.rds       # ENCODE rE2G recurrence filtered ranges (300)
│   ├── re2g_filtered_200.rds       # ENCODE rE2G recurrence filtered ranges (200)
│   ├── re2g_filtered_100.rds       # ENCODE rE2G recurrence filtered ranges (100)
│   ├── re2g_filtered_50.rds        # ENCODE rE2G recurrence filtered ranges (50)
│   └── re2g_cell_type/             # Per-biosample cell-type-specific databases
└── README.md
```

Only scripts are tracked by git.

---

## Data download and setup

### 1. ReMap2022 peak files

Download the per-experiment region files from the
[ReMap2022 portal](https://remap2022.univ-amu.fr/) (human hg38 assembly).
Place all peak files inside `data/remap_hg38/`.

### 2. ENCODE rE2G enhancer–gene predictions

Search all data files in ENCODE: "ENCODE-rE2G predictions of enhancer-gene 
regulatory interactions". Download them, merge in one table and save it as
`data/all_encode_re2g.txt`. The file uses `|` to separate the gene symbol from
a tab-delimited `<genomic range>\t<score>` field.

Per-biosample BED.gz files from the ENCODE portal are downloaded in
`data/encode_re2g_celltype/`. File names must contain the ENCODE accession
pattern `ENCFF###XXX` (e.g. `ENCFF123ABC.bed.gz`).

### 3. GeneHancer / CREdb object

The `ReMap2022+EnsTSS+GH.Rdata` file is distributed alongside the previous
[TFEA.ChIP package downloads](https://github.com/LauraPS1/TFEA.ChIP_downloads).
GeneHancer PMID reference: 28605766

CREdb data is found at Zenodo (zenodo.11407329), as indicated in the
associated paper (PMID: 39014503).

### 4. R packages

Install all dependencies once from CRAN and Bioconductor:

```r
# CRAN
install.packages(c("dplyr", "ggplot2", "stringr", "purrr", "pals"))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "GenomicRanges",
    "IRanges",
    "limma",
    "rtracklayer",
    "biomaRt",
    "TxDb.Hsapiens.UCSC.hg38.knownGene",
    "org.Hs.eg.db",
    "TFEA.ChIP",
    "babelgene"
))
```

---

## Outputs

### General (pan-cell-type) databases — `out_rdata/`

| File | Regulatory regions |
|------|--------------------|
| `encode_re2g_object_hg38.rda`           | rE2G + promoters |
| `encode_re2g-75_object_hg38.rda`        | rE2G score > 0.75 + promoters |
| `encode_re2g-5_object_hg38.rda`         | rE2G score > 0.50 + promoters |
| `encode_re2g-25_object_hg38.rda`        | rE2G score > 0.25 + promoters |
| `encode_re2g-300depth_object_hg38.rda`  | depth-filtered rE2G (≥300 reads) + promoters |
| `encode_re2g-200depth_object_hg38.rda`  | depth-filtered rE2G (≥200 reads) + promoters |
| `encode_re2g-100depth_object_hg38.rda`  | depth-filtered rE2G (≥100 reads) + promoters |
| `encode_re2g-50depth_object_hg38.rda`   | depth-filtered rE2G (≥50 reads) + promoters |
| `credb_object_hg38.rda`                 | CREdb + promoters |
| `genehancer_object.rda`                 | GeneHancer + promoters |

`genehancer_object.rda` is the GeneHancer database from 
`ReMap2022+EnsTSS+GH.Rdata`, saved under a name consistent with the rest of 
the database files.
Each `.rda` file contains a single object called `res` — load it with:

```r
load("out_rdata/encode_re2g-75_object_hg38.rda")  
gene_universe <- ChIPDB[[1]]   # character vector of Entrez Gene IDs
chip_index    <- ChIPDB[[2]]   # named list: experiment ID -> integer indices
```

### Cell-type-specific databases — `out_rdata/re2g_cell_type/`

One pair of `.rda` / `.rds` files per ENCODE biosample, named
`<cell_type>_<ENCFF_ID>_object_hg38.{rda,rds}`.

---

## Pre-built databases on Zenodo

The fully computed databases are available for direct download:

| Dataset | `.rda` DOI | Parquet DOI |
|---------|-----------|-------------|
| Curated general databases | [10.5281/zenodo.15667485](https://doi.org/10.5281/zenodo.23237047) | [10.5281/zenodo.23192462](https://doi.org/10.5281/zenodo.23237198) |
| Cell-type-specific — organs & tissues (Part 1) | [10.5281/zenodo.17192067](https://doi.org/10.5281/zenodo.17192067) | [10.5281/zenodo.23191430](https://doi.org/10.5281/zenodo.23191430) |
| Cell-type-specific — organs & tissues (Part 2) | [10.5281/zenodo.17192175](https://doi.org/10.5281/zenodo.17192175) | [10.5281/zenodo.23191960](https://doi.org/10.5281/zenodo.23191960) |
| Cell-type-specific — stem, progenitor & iPS cells | [10.5281/zenodo.17191689](https://doi.org/10.5281/zenodo.17191689) | [10.5281/zenodo.23192059](https://doi.org/10.5281/zenodo.23192059) |
| Cell-type-specific — cancer & immortalized cell lines | [10.5281/zenodo.17177454](https://doi.org/10.5281/zenodo.17177454) | [10.5281/zenodo.23192306](https://doi.org/10.5281/zenodo.23192306) |
| Cell-type-specific — connective, mesenchymal & muscle | [10.5281/zenodo.17191999](https://doi.org/10.5281/zenodo.17191999) | [10.5281/zenodo.23192216](https://doi.org/10.5281/zenodo.23192216) |
| Cell-type-specific — primary blood & immune cells | [10.5281/zenodo.17191614](https://doi.org/10.5281/zenodo.17191614) | [10.5281/zenodo.23192383](https://doi.org/10.5281/zenodo.23192383) |
| Cell-type-specific — epithelial & endothelial cells | [10.5281/zenodo.17191957](https://doi.org/10.5281/zenodo.17191957) | [10.5281/zenodo.23192127](https://doi.org/10.5281/zenodo.23192127) |

---

## Citation

If you use these databases or scripts, please cite this repository and the
underlying resources: ReMap2022, ENCODE rE2G, GeneHancer, CREdb, TFEA.ChIP, 
and the Zenodo deposits above.
