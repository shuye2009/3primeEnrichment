# compute 3prime enrichment score and plot non-pooled scores

scriptd <- dirname(this.path::this.path())
source(file.path(scriptd, "3primeEnrichment_score.R"))

bfc <- BiocFileCache()
url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz"
gz <- bfcrpath(bfc, url)
gtfFile <- gsub(".gz", "", gz)
if(!file.exists(gtfFile)) gtfFile <- R.utils::gunzip(gz, remove=FALSE)

if(1){
args <- commandArgs(trailingOnly = T)
wd <- args[1]
setwd(wd)

treatBam <- args[2] 
names(treatBam) <- args[3]
controlBam <- args[4] 
names(controlBam) <- args[5] 
genome <- args[6]
nc <- as.integer(args[7])
stranded <- as.logical(args[8])
}

if(0){
   
   wd <-  "C:/GREENBLATT/Zuyao/RNAseq/Nonstranded"
   setwd(wd)
   
   treatBam <- "C:/GREENBLATT/Zuyao/RNAseq/Nonstranded/Greenblatt_1_DKO_2_ATCACG/Aligned.sortedByCoord.out.bam"
   names(treatBam) <- "DKO_2"
   controlBam <- "C:/GREENBLATT/Zuyao/RNAseq/Nonstranded/Greenblatt_10_Scramble_Cri_2_GAGTGG/Aligned.sortedByCoord.out.bam"
   names(controlBam) <- "Scramble_2"
   genome <- "hg38"
   nc <- 20
   stranded <- FALSE
}


handleInputParams <- list(CLIP_reads=FALSE, fix_width=0, fix_point="start", norm=FALSE, useScore=FALSE,
                          outRle=FALSE, useSizeFactor=FALSE, genome=genome)

df <- threePrimeEnrichmentScore(treatBam, controlBam, gtfFile, stranded, handleInputParams, nc)
