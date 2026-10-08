# =============================================================================
# 03_filter_re2g_by_depth.R
# Purpose : Filter the full ENCODE rE2G enhancer–gene set to retain only
#           high-confidence regions based on "cross-biosample depth" — i.e.
#           how many rE2G biosample files place a peak over the same locus for
#           the same gene. Regions consistently supported by many biosamples
#           are more likely to be genuine regulatory elements.
#
#           The approach per gene:
#             1. Pool all rE2G ranges for that gene across biosamples.
#             2. Merge overlapping ranges (reduce).
#             3. Compute per-base coverage across the merged regions to find
#                how many original ranges pile up at each position.
#             4. Summarise each merged region by its maximum coverage value
#                (i.e. the peak depth at the most-supported position).
#             5. Retain only regions whose max coverage meets a threshold.
#           Four thresholds are saved: ≥50, ≥100, ≥200, ≥300 biosamples.
#
# Input   : out_rdata/encodere2g_hg38_ranges.rds
#             – the full rE2G GRanges object built upstream.
#
# Outputs :
#   out_rdata/re2g_filtered_50.rds
#   out_rdata/re2g_filtered_100.rds
#   out_rdata/re2g_filtered_200.rds
#   out_rdata/re2g_filtered_300.rds
#
# =============================================================================

library(GenomicRanges)  
library(IRanges)        
library(dplyr)           

# ----- 1. Load the full rE2G ranges -------------------------------
re2g_only <- readRDS(file.path("out_rdata", "encodere2g_hg38_ranges.rds"))

# Drop strand information — rE2G regions are not strand-specific
strand(re2g_only) <- "*"

# Remove names so c() operations downstream stay clean
names(re2g_only) <- NULL

# ----- 2. Split by gene_id ---------------------------------------------------
# The coverage computation is done per-gene so that overlapping ranges from
# different genes are never mixed.
df_re2g <- as.data.frame(re2g_only)
split_by_gene <- split(df_re2g, df_re2g$gene_id)

message("Number of genes: ", length(split_by_gene))

# ----- 3. Per-gene coverage computation --------------------------------------
# For each gene:
#   a) Build a GRanges from all its rE2G entries.
#   b) Reduce to non-overlapping merged regions.
#   c) Compute per-base coverage (how many original ranges overlap each base).
#   d) For each merged region, find the maximum coverage and the positions
#      that achieve it — this becomes the representative "peak" interval.
#   e) Return a GRanges with gene_id and max_coverage metadata.

cc <- 0L  # progress counter 

overlap_counts <- lapply(split_by_gene, function(gene_df) {

  cc <<- cc + 1L
  if (cc %% 1000L == 0L) {
    message("  Processing gene ", cc, " of ", length(split_by_gene))
  }

  # Build GRanges for this gene
  gr <- GRanges(
    seqnames = gene_df$seqnames,
    ranges   = IRanges(start = gene_df$start, end = gene_df$end),
    gene_id  = gene_df$gene_id
  )

  # Merge overlapping ranges -> non-redundant set of loci for this gene
  gr_reduced <- GenomicRanges::reduce(gr)

  # Per-base coverage across the original (unreduced) ranges
  cov_obj <- coverage(gr)

  output <- vector("list", length(gr_reduced))

  for (i in seq_along(gr_reduced)) {
    region <- gr_reduced[i]
    chr    <- as.character(seqnames(region))

    # Skip chromosomes not present in the coverage object
    if (!chr %in% names(cov_obj)) next

    # Extract per-base coverage vector over this merged region
    view     <- Views(cov_obj[[chr]],
                      start = start(region),
                      end   = end(region))
    cov_vals <- as.numeric(view[[1]])

    max_cov <- max(cov_vals)
    if (max_cov == 0L) next

    # Identify the contiguous span of bases that achieve max coverage.
    # Taking min/max of their positions gives a tight "peak" interval.
    max_pos   <- which(cov_vals == max_cov)
    new_start <- start(region) + min(max_pos) - 1L
    new_end   <- start(region) + max(max_pos) - 1L

    output[[i]] <- GRanges(
      seqnames     = seqnames(region),
      ranges       = IRanges(start = new_start, end = new_end),
      gene_id      = gene_df$gene_id[1],
      max_coverage = max_cov
    )
  }

  # Collapse non-NULL results into a single GRanges for this gene
  do.call(c, output[!sapply(output, is.null)])
})

message("Coverage computation complete.")

# ----- 4. Apply depth thresholds and save ------------------------------------
# Each threshold keeps only merged regions whose maximum coverage (number of
# biosample-level rE2G peaks piling up) meets the cutoff.
# The filtered objects contain *enhancer regions only*; promoters are added
# back in script 04 when each database variant is assembled.

out_dir <- "out_rdata"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

filter_and_save <- function(counts_list, min_cov, filename) {
  filtered <- lapply(counts_list, function(x) {
    if (is.null(x) || length(x) == 0) return(NULL)
    x[x$max_coverage >= min_cov]
  })
  filtered_gr <- do.call(c, filtered[!sapply(filtered, is.null)])
  names(filtered_gr) <- NULL
  message("Threshold ≥", min_cov, ": ", length(filtered_gr), " regions retained.")
  saveRDS(filtered_gr, file = file.path(out_dir, filename))
  message("Saved: ", file.path(out_dir, filename))
  invisible(filtered_gr)
}

filter_and_save(overlap_counts, 300L, "re2g_filtered_300.rds")
filter_and_save(overlap_counts, 200L, "re2g_filtered_200.rds")
filter_and_save(overlap_counts, 100L, "re2g_filtered_100.rds")
filter_and_save(overlap_counts,  50L, "re2g_filtered_50.rds")

message("All depth-filtered rE2G objects saved to out_rdata/")

