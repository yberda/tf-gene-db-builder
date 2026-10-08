# =============================================================================
# 05_build_celltype_databases.R
# Purpose : Build cell-type-specific TF target databases by processing each
#           ENCODE rE2G BED file independently. Each file corresponds to one
#           biosample; enhancer–gene links are extracted, filtered, and
#           overlapped with ReMap2022 ChIP-seq peaks via makeChIPGeneDB().
#
# Inputs  :
#   out_rdata/remap_hg38.rds               – from script 01
#   data/encode_re2g_celltype/             – directory of per-biosample
#                                            ENCODE rE2G BED.gz files)
#
# Outputs :
#   out_rdata/re2g_cell_type/<cell_type>_<ENCFF_ID>_object_hg38.rda
#   out_rdata/re2g_cell_type/<cell_type>_<ENCFF_ID>_object_hg38.rds
#   (one pair of files per biosample)
#
# Notes   :
#   - The script tracks which files have already been processed and skips them,
#     so it is safe to re-run after interruption.
#   - BED files use 0-based coordinates; start positions are incremented by 1
#     to convert to 1-based (R/Bioconductor convention).
# =============================================================================

library(limma)          
library(dplyr)           
library(GenomicRanges)   
library(TFEA.ChIP)      
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  
library(org.Hs.eg.db)   
library(stringr)         

# ----- 1. Load shared inputs -------------------------------------------------
remap <- readRDS(file.path("out_rdata", "remap_hg38.rds"))

# Canonical Entrez gene IDs from hg38 knownGene track
txdb      <- TxDb.Hsapiens.UCSC.hg38.knownGene
knownGene <- genes(txdb)$gene_id

# ----- 2. Identify files to process ------------------------------------------
re2g_dir <- file.path("data", "encode_re2g_celltype")

all_files <- list.files(re2g_dir, pattern = "\\.bed\\.gz$", full.names = TRUE)

# Extract ENCODE file IDs (e.g. ENCFF123ABC) from file paths
all_ids <- str_extract(all_files, "ENCFF\\d{3}[A-Z]{3}")

# Check which IDs have already been saved to skip them on re-run
out_cell_dir <- file.path("out_rdata", "re2g_cell_type")
dir.create(out_cell_dir, showWarnings = FALSE, recursive = TRUE)

done_files <- list.files(out_cell_dir, pattern = "\\.rds$", full.names = TRUE)
done_ids   <- str_extract(done_files, "ENCFF\\d{3}[A-Z]{3}")

pending_ids   <- setdiff(all_ids, done_ids)
files_to_run  <- all_files[all_ids %in% pending_ids]

message(length(all_files),   " total files found.")
message(length(done_ids),    " already processed – skipping.")
message(length(files_to_run), " files to process.")

# ----- 3. Define per-file processing function --------------------------------
process_file <- function(file, knownGene, remap, out_cell_dir) {

  file_id <- str_extract(file, "ENCFF\\d{3}[A-Z]{3}")
  message("\n--- Processing: ", basename(file), " ---")

  # 3a. Read the BED file (tab-separated, no header)
  df <- read.table(file, sep = "\t", stringsAsFactors = FALSE)

  if (nrow(df) == 0) {
    message("  Skipping: file is empty.")
    return(NULL)
  }

  # Expected BED columns used:
  #   V1  = chr
  #   V2  = start (0-based)
  #   V3  = end
  #   V7  = target gene (Ensembl gene ID)
  #   V10 = biosample / cell type name
  enh <- data.frame(
    gene      = df$V7,
    chr       = df$V1,
    start     = df$V2 + 1,   # convert 0-based BED → 1-based R convention
    end       = df$V3,
    cell_type = df$V10,
    stringsAsFactors = FALSE
  )

  # 3b. Map Ensembl IDs to Entrez Gene IDs
  entrez_ids <- mapIds(
    org.Hs.eg.db,
    keys      = enh$gene,
    column    = "ENTREZID",
    keytype   = "ENSEMBL",
    multiVals = "first"
  )
  enh$entrez_id <- as.character(unname(entrez_ids[enh$gene]))

  # 3c. Filter: remove unmapped genes and genes outside knownGene
  enh <- enh[!is.na(enh$entrez_id), ]
  enh <- enh[enh$entrez_id %in% knownGene, ]

  # 3d. Build GRanges object
  gr <- GRanges(
    seqnames  = enh$chr,
    ranges    = IRanges(start = enh$start, end = enh$end),
    gene_id   = enh$entrez_id,
    cell_type = enh$cell_type
  )

  # 3e. Overlap with ReMap2022 and build TF–gene database
  res <- makeChIPGeneDB(gr, remap)

  # 3f. Save outputs (both .rda and .rds for flexibility)
  cell_type_name <- unique(enh$cell_type)[1]
  # Sanitise cell type name for use in file paths
  cell_type_safe <- gsub("[^A-Za-z0-9_.-]", "_", cell_type_name)

  base_name <- paste0(cell_type_safe, "_", file_id, "_object_hg38")
  rds_path  <- file.path(out_cell_dir, paste0(base_name, ".rds"))
  rda_path  <- file.path(out_cell_dir, paste0(base_name, ".rda"))

  saveRDS(res, file = rds_path)
  
  ChIPDB <- res
  save(ChIPDB, file = rda_path)

  message("  Saved: ", rds_path)
  return(rds_path)
}

# ----- 4. Run across all pending files ---------------------------------------
results <- lapply(
  files_to_run,
  process_file,
  knownGene    = knownGene,
  remap        = remap,
  out_cell_dir = out_cell_dir
)

# Report
n_success <- sum(!sapply(results, is.null))
message("\nDone. ", n_success, " / ", length(files_to_run), " files processed successfully.")
message("Outputs in: ", out_cell_dir)
