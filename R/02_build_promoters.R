# =============================================================================
# 02_build_promoters.R
# Purpose : Define hg38 gene promoter regions (±1000 bp around each TSS) and
#           annotate them with Entrez Gene IDs. The resulting GRanges object
#           is used in subsequent scripts as the "baseline" regulatory region
#           for every gene (combined later with enhancer–gene links).
#
# Input   : None (all data come from Bioconductor annotation packages).
#
# Output  : out_rdata/prom-hg38_entrez.rds
#
# Annotation packages required (install once via BiocManager):
#   TxDb.Hsapiens.UCSC.hg38.knownGene
#   org.Hs.eg.db
# =============================================================================

library(GenomicRanges)                      
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  
library(org.Hs.eg.db)                       

# ----- 1. Load transcript annotation database --------------------------------
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# ----- 2. Define promoter regions --------------------------------------------
# ±1000 bp around each annotated transcription start site (TSS).
# 'promoters()' returns one range per transcript; transcripts sharing a gene
# may produce overlapping ranges – this is intentional and handled downstream.
promoters_gr <- promoters(txdb, upstream = 1000, downstream = 1000)

message("Total promoter ranges (all transcripts): ", length(promoters_gr))

# ----- 3. Map UCSC transcript IDs to Entrez Gene IDs ------------------------
# The names() of the promoters GRanges are UCSC Known Gene transcript IDs.
ucsc_ids <- names(promoters_gr)

entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys     = ucsc_ids,
  column   = "ENTREZID",  
  keytype  = "UCSCKG",     # source ID type (UCSC Known Gene)
  multiVals = "first"     
)

# Store as metadata column; drop the default tx_id / tx_name columns
mcols(promoters_gr)$tx_id   <- NULL
mcols(promoters_gr)$tx_name <- NULL
mcols(promoters_gr)$gene_id <- entrez_ids   # renamed for consistency

# ----- 4. Filter to canonical chromosomes ------------------------------------
# Retain only autosomes and sex chromosomes 
canonical_chrs <- c(paste0("chr", 1:22), "chrX", "chrY")

promoters_gr <- promoters_gr[seqnames(promoters_gr) %in% canonical_chrs]
message("After chromosome filter: ", length(promoters_gr))

# ----- 5. Remove promoters without an Entrez ID ------------------------------
n_na <- sum(is.na(promoters_gr$gene_id))
message("Promoters with NA Entrez ID (to be removed): ", n_na)

promoters_gr <- promoters_gr[!is.na(promoters_gr$gene_id)]
message("Final promoter count: ", length(promoters_gr))

# ----- 6. Save output --------------------------------------------------------
out_dir <- "out_rdata"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(promoters_gr, file = file.path(out_dir, "prom-hg38_entrez.rds"))
message("Saved: out_rdata/prom-hg38_entrez.rds")
