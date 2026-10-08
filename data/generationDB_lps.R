# Database generation script originally authored by Laura Puente Santamaría.
# Produces the "GeneHancer" database for use with the first version of 
# TFEA.ChIP, and serves as a reference example for generating the 
# "ReMap2022+EnsTSS+GH.Rdata" database.
#
# Original paper:
#
#   Puente-Santamaría L, Wasserman WW, Del Peso L.
#   TFEA.ChIP: a tool kit for transcription factor binding site enrichment
#   analysis capitalizing on ChIP-seq datasets.
#   Bioinformatics. 2019 Dec 15;35(24):5339–5340.
#   doi: 10.1093/bioinformatics/btz573 | PMID: 31347689

library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

#### enhancers GeneHancer ####
geneHancer <- read.table(
    file="~/Regulatory Networks/contestación paper/genehancer.csv",
    sep =",", header=T, stringsAsFactors = F)


enh <- geneHancer[ geneHancer$feature.name == "Enhancer",]

atr <- strsplit(enh$attributes, ";")

enh_atr <- lapply(
    seq_along(enh$chrom),
    function(i, enh){
        
        enhID <- substr( atr[[i]][1], 15, nchar(atr[[i]][1]) )
        enhChr <- enh$chrom[i]
        enhStart <- enh$start[i]
        enhEnd <- enh$end[i]
        enhScore <- enh$score[i]
        
        row_ <- atr[[i]][-1] # remove ID
        
        genes <- row_[ c(T,F) ]
        scores <- row_[ c(F,T) ]
        
        genes <- substr( genes, 16, nchar(genes) )
        scores <- as.numeric(substr( scores, 7, nchar(scores) ))
        
        res <- data.frame(
            Chr = enhChr, Start = enhStart, End = enhEnd,
            ID = enhID, EnhancerScore = enhScore,
            Gene = genes, GeneScore=scores, stringsAsFactors = F)
        return(res)
    },
    enh = enh
)

enh_atr <- do.call( rbind, enh_atr )

GH_hg19_coords <- read.table(
    file="~/Regulatory Networks/contestación paper/GeneHancer_hg19_coordinates.bed",
    sep="\t",header = F, stringsAsFactors = F)

tmp <- enh_atr
tmp$Start <- GH_hg19_coords[ match(tmp$ID, GH_hg19_coords$V4), ]$V2
tmp$End <- GH_hg19_coords[ match(tmp$ID, GH_hg19_coords$V4), ]$V3
tmp <- tmp[ !is.na(tmp$Start),]

enh_atr_ <- tmp[ tmp$GeneScore > 10 & tmp$EnhancerScore > 1, ]

enh_atr_$Gene <- toupper(enh_atr_$Gene)

symbols <- enh_atr_[ !grepl("^ENSG", enh_atr_$Gene), ]
ensembl <- enh_atr_[ grepl("^ENSG", enh_atr_$Gene), ]

symbols_ <- TFEA.ChIP::GeneID2entrez( symbols$Gene, return.Matrix = T )
ensembl_ <-  TFEA.ChIP::GeneID2entrez( ensembl$Gene, return.Matrix = T )

symbols$entrez <- symbols_$ENTREZ.ID
ensembl$entrez <- ensembl_$ENTREZ.ID

enh_atr_ <- rbind(symbols, ensembl)

enh_atr_ <- enh_atr_[ !is.na(enh_atr_$entrez),]

# enh_gr <- makeGRangesFromDataFrame(
#     enh_atr_10, seqnames.field = "Chr",
#     start.field = "Start", end.field = "End",
#     keep.extra.columns = T
# )
# save(enh_gr, compress = "bzip2",
#     file="~/Regulatory Networks/contestación paper/enhancers_GeneHancer.RData")

enh_gr <- GenomicRanges::makeGRangesFromDataFrame(
    enh_atr_, seqnames.field = "Chr",
    start.field = "Start", end.field = "End",
    keep.extra.columns = F
)
enh_gr$gene_id <- as.character(enh_atr_$entrez)
#### TFBS ####

tmp_dhs_2add<-read.table(
    file="~/Regulatory Networks/OpenDnasa-All/DNase clusters.txt",
    header = T,sep="\t",stringsAsFactors = F)

tmp_dhs_2add<-makeGRangesFromDataFrame(tmp_dhs_2add,
                                        ignore.strand=T,
                                        seqnames.field="chrom",
                                        start.field="chromStart",
                                        end.field="chromEnd")

Genes<-genes(TxDb.Hsapiens.UCSC.hg19.knownGene)

tmp_1 <- findOverlaps(tmp_dhs_2add, Genes, maxgap = 1000 )
qhits_1<-queryHits(tmp_1)
genesCercanos_1 <- subjectHits(tmp_1)
genes_1 <- Genes[ genesCercanos_1 ]$gene_id
DHS_1 <- tmp_dhs_2add[ qhits_1 ]
mcols(DHS_1)$gene_id <- genes_1

DHS_db_1 <- GRangesList( DHS_1, enh_gr )
DHS_db_1 <- unlist(DHS_db_1, recursive = TRUE, use.names = TRUE)

load("~/Regulatory Networks/TFBS updates/mar.2018/Encode+geo(all)_GR.Rdata")
GeneNames<-Genes$gene_id

TFBS_db1 <- TFEA.ChIP::GR2tfbs_db( DHS_db_1, all_GRs, outputAsVector = TRUE)
Mat01 <- TFEA.ChIP::makeTFBSmatrix(GeneList = GeneNames, id_db = TFBS_db1, geneSetAsInput = FALSE)

MetaData <- read.table(
    file = "~/Regulatory Networks/TFBS updates/mar.2018/Metadata.csv",
    header=T, sep=",", stringsAsFactors = F )
MetaData <- MetaData[ MetaData$Accession %in% names(TFBS_db1), 2:7]

# save
save(Mat01, MetaData, file = "~/Regulatory Networks/contestación paper/TFBS_GeneHancer/TSS1Kb+ds.RData", compress = "bzip2" )

