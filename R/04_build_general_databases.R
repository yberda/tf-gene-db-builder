# =============================================================================
# 04_build_general_databases.R
# Purpose : Build curated, general (not cell-type-specific) TF target
#           databases by overlapping ReMap2022 ChIP-seq peaks with regulatory
#           regions defined as:
#             promoters  – ±1000 bp around each TSS (from script 02)
#             enhancers  – ENCODE rE2G predictions filtered at several score
#                          thresholds, OR CRE-DB enhancer–gene links
#           For each combination the script calls TFEA.ChIP::makeChIPGeneDB(),
#           which returns a list with:
#             [[1]] character vector of unique Entrez Gene IDs
#             [[2]] named list: ChIP-seq experiment → integer indices into [[1]]
#
# Inputs  :
#   out_rdata/remap_hg38.rds                 – from script 01
#   out_rdata/prom-hg38_entrez.rds           – from script 02
#   data/all_encode_re2g.txt                 – ENCODE rE2G enhancer–gene table
#                                              (downloaded from ENCODE)
#   ReMap2022+EnsTSS+GH.Rdata               – GeneHancer object
#                                              (from TFEA.ChIP v1 GitHub)
#   out_rdata/re2g_filtered_300.rds  \
#   out_rdata/re2g_filtered_200.rds   |– depth-filtered rE2G GRanges
#   out_rdata/re2g_filtered_100.rds   |  (pre-computed; see README)
#   out_rdata/re2g_filtered_50.rds   /
#
# Outputs : one .rda file per database variant (see bottom of script)
# =============================================================================

library(limma)           
library(dplyr)           
library(GenomicRanges)   
library(TFEA.ChIP)      
library(TxDb.Hsapiens.UCSC.hg38.knownGene) 
library(org.Hs.eg.db)  

out_dir <- "out_rdata"

# ----- 0. Helper: save as .rda with a fixed object name ---------------------
save_db <- function(obj, path) {
  ChIPDB <- obj                          
  save(ChIPDB, file = path)
  message("Saved: ", path)
}

# ----- 1. Load shared inputs -------------------------------------------------
remap <- readRDS(file.path("out_rdata", "remap_hg38.rds"))
prom  <- readRDS(file.path("out_rdata", "prom-hg38_entrez.rds"))

# Canonical Entrez gene IDs from hg38 knownGene track
txdb      <- TxDb.Hsapiens.UCSC.hg38.knownGene
knownGene <- genes(txdb)$gene_id

# GeneHancer object (shipped with TFEA.ChIP)
load(file.path("data", "ReMap2022+EnsTSS+GH.Rdata"))  

# ----- 2. Load and process ENCODE rE2G data ----------------------------------
# The file uses '|' to separate the gene symbol from the genomic range + score,
# which are themselves tab-separated.
re2g_raw <- read.table(
  file.path("data", "all_encode_re2g.txt"),
  sep             = "|",
  stringsAsFactors = FALSE
)
colnames(re2g_raw) <- c("gene_symbol", "range_score")

# Split the second column into range and score
tmp <- limma::strsplit2(re2g_raw$range_score, "\t")
re2g_df <- data.frame(
  gene  = re2g_raw$gene_symbol,
  range = tmp[, 1],
  score = as.numeric(tmp[, 2]),
  stringsAsFactors = FALSE
)

# Map gene symbols to Entrez IDs
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys      = re2g_df$gene,
  column    = "ENTREZID",
  keytype   = "SYMBOL",
  multiVals = "first"
)

# Build a GRanges object with gene_id and score metadata
enh_all <- GRanges(
  re2g_df$range,
  gene_id = entrez_ids,
  score   = re2g_df$score
)

# Remove entries without an Entrez ID or outside the knownGene universe
enh_all <- enh_all[!is.na(enh_all$gene_id)]
enh_all <- enh_all[enh_all$gene_id %in% knownGene]


message("rE2G entries after gene filter: ", length(enh_all))


# ----- 3. Load and process CREdb data ----------------------------------

# Read enhancer-gene interaction
enh.gene.int <- read.table(gzfile('CREdb_data/enhancer_gene_interaction.csv.gz'),
                      sep = ',', header = TRUE)
enh.gene.int <- enh.gene.int[,c('enhancer_id', 'gene_id')]

# Read combined elements to add start and end of enhancer to the previous table
comb.elem <- read.table(gzfile('CREdb_data/combined_elements.csv.gz'),
                                          sep = ',', header = TRUE)
comb.elem <- comb.elem[,c('chromosome', 'start_position', 'end_position', 'element_type', 'element_id')]

# Merge two tables
all.info <- merge(enh.gene.int, comb.elem, by.x = 'enhancer_id', by.y = 'element_id', all.x = TRUE)
all.info %>% head

# Add entrez id
entrez_ids <- mapIds(org.Hs.eg.db, keys = all.info$gene_id, 
              column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
all.info$entrez_id <- entrez_ids
all.info <- all.info[!is.na(all.info$start_position), ]

# Now generate genomic ranges
credb <- GRanges(seqnames = all.info$chromosome,
                      ranges = IRanges(all.info$start_position, all.info$end_position),
                      gene_id = all.info$entrez_id)
                      
# Remove entries without an Entrez ID or outside the knownGene universe
credb <- credb[!is.na(credb$gene_id)]
credb <- credb[credb$gene_id %in% knownGene]


# ----- 4. Apply score thresholds ---------------------------------------------
enh_75  <- enh_all[enh_all$score > 0.75]
enh_50  <- enh_all[enh_all$score > 0.50]
enh_25  <- enh_all[enh_all$score > 0.25]

# ----- 5. Load depth-filtered rE2G variants ----------------------------------
re2g_300 <- readRDS(file.path("out_rdata", "re2g_filtered_300.rds"))
re2g_200 <- readRDS(file.path("out_rdata", "re2g_filtered_200.rds"))
re2g_100 <- readRDS(file.path("out_rdata", "re2g_filtered_100.rds"))
re2g_50  <- readRDS(file.path("out_rdata", "re2g_filtered_50.rds"))

# ----- 6. Combine enhancers + promoters into regulatory region sets ----------
# Each object = enhancer GRanges + promoter GRanges.

reg_re2g     <- c(enh_all,   prom)
reg_re2g_75  <- c(enh_75,    prom)
reg_re2g_50  <- c(enh_50,    prom)
reg_re2g_25  <- c(enh_25,    prom)
reg_depth300 <- c(re2g_300,  prom)
reg_depth200 <- c(re2g_200,  prom)
reg_depth100 <- c(re2g_100,  prom)
reg_depth50  <- c(re2g_50,   prom)
reg_credb    <- c(credb,     prom)  

# ----- 7. Overlap with ReMap2022 and build TF–gene databases -----------------
# makeChIPGeneDB() finds which ChIP-seq peaks overlap each regulatory region
# and records the corresponding gene index for each experiment.

message("Building database: rE2G + promoters …")
res_re2g <- makeChIPGeneDB(reg_re2g, remap)
save_db(res_re2g, file.path(out_dir, "encode_re2g_object_hg38.rda"))
# also save ranges
saveRDS(enh_all, "out_rdata/encodere2g_hg38_ranges.rds")

message("Building database: rE2G score > 0.75 + promoters …")
res_re2g_75 <- makeChIPGeneDB(reg_re2g_75, remap)
save_db(res_re2g_75, file.path(out_dir, "encode_re2g-75_object_hg38.rda"))

message("Building database: rE2G score > 0.50 + promoters …")
res_re2g_50 <- makeChIPGeneDB(reg_re2g_50, remap)
save_db(res_re2g_50, file.path(out_dir, "encode_re2g-5_object_hg38.rda"))

message("Building database: rE2G score > 0.25 + promoters …")
res_re2g_25 <- makeChIPGeneDB(reg_re2g_25, remap)
save_db(res_re2g_25, file.path(out_dir, "encode_re2g-25_object_hg38.rda"))

message("Building database: depth-filtered rE2G (≥300 reads) + promoters …")
res_depth300 <- makeChIPGeneDB(reg_depth300, remap)
save_db(res_depth300, file.path(out_dir, "encode_re2g-300depth_object_hg38.rda"))

message("Building database: depth-filtered rE2G (≥200 reads) + promoters …")
res_depth200 <- makeChIPGeneDB(reg_depth200, remap)
save_db(res_depth200, file.path(out_dir, "encode_re2g-200depth_object_hg38.rda"))

message("Building database: depth-filtered rE2G (≥100 reads) + promoters …")
res_depth100 <- makeChIPGeneDB(reg_depth100, remap)
save_db(res_depth100, file.path(out_dir, "encode_re2g-100depth_object_hg38.rda"))

message("Building database: depth-filtered rE2G (≥50 reads) + promoters …")
res_depth50 <- makeChIPGeneDB(reg_depth50, remap)
save_db(res_depth50, file.path(out_dir, "encode_re2g-50depth_object_hg38.rda"))

message("Building database: CREdb + promoters …")
res_credb <- makeChIPGeneDB(reg_credb, remap)
save_db(res_credb, file.path(out_dir, "credb_object_hg38.rda"))

# GeneHancer database (loaded from ReMap2022+EnsTSS+GH.Rdata)
# res object already exists inside that file as ChIPDB; re-save under a 
# consistent name
message("Saving original GeneHancer database …")
save(ChIPDB, file.path(out_dir, "genehancer_object.rda"))

message("All general databases saved to out_rdata/")
