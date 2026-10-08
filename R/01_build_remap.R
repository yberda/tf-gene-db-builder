# =============================================================================
# 01_build_remap.R
# Purpose : Load ReMap2022 ChIP-seq peak files (All peaks, hg38) 
#           and convert them to a named list of GRanges objects
#           (one entry per TF experiment).
#
# Input   : data/remap_hg38/  – directory containing one .txt file per TF
#           experiment. Each file has four tab-separated columns:
#             chr | start | end | experiment_id
#
# Output  : out_rdata/remap_hg38.rds
#
# Notes   : ReMap2022 files can be downloaded from https://remap2022.univ-amu.fr/
#           Place the per-experiment region files under data/remap_hg38/.
# =============================================================================

library(GenomicRanges)   
library(dplyr)          
library(limma)          

# ----- 1. Locate input files -------------------------------------------------
remap_dir <- file.path("data", "remap_hg38")
filepath <- list.files(remap_dir, full.names = TRUE, pattern = "\\.txt$")

# ----- 2. Read each file into a data frame -----------------------------------
# Expected columns: chr, start, end, experiment_id
remap_raw <- lapply(filepath, read.table,
                    header    = FALSE,
                    sep       = "\t",
                    col.names = c("chr", "start", "end", "experiment_id"))

# ----- 3. Convert each data frame to a GRanges object -----------------------
# The experiment_id stored in the metadata comes from the first row of each
# file (all rows in a single file share the same ID).
granges_list <- lapply(remap_raw, function(df) {
  GRanges(
    seqnames     = df$chr,
    ranges       = IRanges(start = df$start, end = df$end),
    experiment_id = df$experiment_id[1]
  )
})

# ----- 4. Name list elements by experiment ID --------------------------------
# File names are expected to follow the pattern: <experiment_id>_regions.txt
experiment_ids <- limma::strsplit2(filepath, split = "/")[, ncol(limma::strsplit2(filepath, split = "/"))]
experiment_ids <- gsub("_regions\\.txt$", "", experiment_ids)

names(granges_list) <- experiment_ids

message("ReMap GRanges list built with ", length(granges_list), " experiments.")

# ----- 5. Save output --------------------------------------------------------
out_dir <- "out_rdata"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(granges_list, file = file.path(out_dir, "remap_hg38.rds"))
message("Saved: out_rdata/remap_hg38.rds")
