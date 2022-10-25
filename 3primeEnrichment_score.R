library(DescTools)
library(GenomicPlot)
library(BiocFileCache)
bfc <- BiocFileCache()
url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz"
gz <- bfcrpath(bfc, url)
gtffile <- gsub(".gz", "", gz)
if(!file.exists(gtffile)) gtffile <- R.utils::gunzip(gz, remove=FALSE)

treatBam <- "C:/GREENBLATT/Zuyao/RNAseq/Nonstranded/Greenblatt_1_DKO_2_ATCACG/Aligned.sortedByCoord.out.bam"
names(treatBam) <- "treat"
controlBam <- "C:/GREENBLATT/Zuyao/RNAseq/Nonstranded/Greenblatt_10_Scramble_Cri_2_GAGTGG/Aligned.sortedByCoord.out.bam"
names(controlBam) <- "control"

handleInputParams <- list(CLIP_reads=FALSE, fix_width=0, fix_point="start", norm=TRUE, useScore=FALSE,
                          outRle=FALSE, useSizeFactor=TRUE, genome="hg38")
threePrimeEnrichmentScore(treatBam, controlBam, gtfFile){
   txdb <- makeTxDbFromGFF(gtfFile)
   gene<- get_genomic_feature_coordinates(txdb, featureName="gene", longest=TRUE, protein_coding=TRUE)
   gene_gr <- gene$GRanges
   
   inputs <- handle_input(c(treatBam, controlBam), handleInputParams)
   
   names(inputs) <- c(names(treatBam), names(controlBam))
}