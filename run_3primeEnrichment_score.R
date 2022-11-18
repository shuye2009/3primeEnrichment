# compute 3prime enrichment score and plot non-pooled scores

scriptd <- dirname(this.path::this.path())
source(file.path(scriptd, "3primeEnrichment_score.R"))

bfc <- BiocFileCache()
url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz"
gz <- bfcrpath(bfc, url)
gtfFile <- gsub(".gz", "", gz)
if(!file.exists(gtfFile)) gtfFile <- R.utils::gunzip(gz, remove=FALSE)

txdb <- makeTxDbFromGFF(gtfFile)

if(0){
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
   handleInputParams <- list(CLIP_reads=FALSE, fix_width=0, fix_point="start", norm=FALSE, useScore=FALSE,
                             outRle=FALSE, useSizeFactor=FALSE, genome=genome)
   
   df <- threePrimeEnrichmentScore(treatBam, controlBam, txdb, stranded, handleInputParams, nc)

}

if(1){
   genome <- "hg38"
   nc <- 40
   stranded <- FALSE
   wd <-  "C:/GREENBLATT/Zuyao/RNAseq/Nonstranded"
   setwd(wd)
   bam <- "Aligned.sortedByCoord.out.bam"
   dirs <- c("Greenblatt_5_RPRD1A_Cri_1_TAGCTT_L007", 
             "Greenblatt_6_RPRD1A_Cri_3_GGCTAC", 
             "Greenblatt_9_Scramble_Cri_1_2_CGTACG_L007",
             "Greenblatt_10_Scramble_Cri_2_GAGTGG",
             "Greenblatt_7_RPRD1B_Cri_5_GTGGCC",
             "Greenblatt_8_RPRD1B_Cri_6_GTTTCG",
             "Greenblatt_1_DKO_2_ATCACG",
             "Greenblatt_2_DKO_3_D2_TTAGGC")
   sampleNames <- c("RPRD1A_1", "RPRD1A_2", "Scramble_1", "Scramble_2", "RPRD1B_1", "RPRD1B_2", "DKO_1", "DKO_2")
   
   sampleBams <- file.path(wd, dirs, bam)
   names(sampleBams) <- sampleNames
   
   handleInputParams <- list(CLIP_reads=FALSE, fix_width=0, fix_point="start", norm=FALSE, useScore=FALSE,
                             outRle=FALSE, useSizeFactor=FALSE, genome=genome)
   
   comparisons <- list(c("RPRD1B_1", "Scramble_1"), c("RPRD1B_2", "Scramble_2"), c("RPRD1A_1", "Scramble_1"), c("RPRD1A_2", "Scramble_2"), c("DKO_1", "Scramble_1"), c("DKO_2", "Scramble_2"), c("DKO_2", "DKO_1"),  c("RPRD1A_2", "RPRD1A_1"), c("RPRD1B_2", "RPRD1B_1"),  c("Scramble_2", "Scramble_1"))
   
   
   for(i in 1:length(comparisons)){
      comparison <- comparisons[[i]]
      df <- threePrimeEnrichmentScore(treatBam=sampleBams[comparison[1]], controlBam=sampleBams[comparison[2]], txdb, stranded, handleInputParams, nc)
   }
   
   n_bins <- 20
   dataFiles_list <- list(RPRD1AScramble_rep2_vs_rep1_pooled=c("RPRD1A_2_vs_RPRD1A_1_3prime_enrichment_score.tab", "Scramble_2_vs_Scramble_1_3prime_enrichment_score.tab"),
                          RPRD1BScramble_rep2_vs_rep1_pooled =c("RPRD1B_2_vs_RPRD1B_1_3prime_enrichment_score.tab", "Scramble_2_vs_Scramble_1_3prime_enrichment_score.tab"),
                          DKOScramble_rep2_vs_rep1_pooled=c("DKO_2_vs_DKO_1_3prime_enrichment_score.tab","Scramble_2_vs_Scramble_1_3prime_enrichment_score.tab"),
                          RPRD1A_vs_Scramble_pooled=c("RPRD1A_1_vs_Scramble_1_3prime_enrichment_score.tab", "RPRD1A_2_vs_Scramble_2_3prime_enrichment_score.tab"),
                          RPRD1B_vs_Scramble_pooled=c("RPRD1B_1_vs_Scramble_1_3prime_enrichment_score.tab", "RPRD1B_2_vs_Scramble_2_3prime_enrichment_score.tab"),
                          DKO_vs_Scramble_pooled=c("DKO_1_vs_Scramble_1_3prime_enrichment_score.tab", "DKO_2_vs_Scramble_2_3prime_enrichment_score.tab"))
   for(dataName in names(dataFiles_list)){
      dataFiles <- dataFiles_list[[dataName]]
      plot_3prime_enrichment(datafm=NULL, dataFiles, n_bins, dataName, stranded)
   }
}

if(1){
   genome <- "hg38"
   nc <- 40
   stranded <- TRUE
   wd <-  "C:/GREENBLATT/Zuyao/RNAseq/Stranded"
   setwd(wd)
   bam <- "Aligned.sortedByCoord.out.bam"
   dirs <- list.dirs(full.names=FALSE, recursive = FALSE)
   sampleNames <- c("Scramble_1", "Scramble_2", "DKO_1", "RPRD1B_GFP_1", "RPRD1B_GFP_2", "RPRD1A_GFP_1", "RPRD1A_GFP_2", "DKO_2", "GFP_1", "GFP_2", "DKO_3")
   
   sampleBams <- file.path(wd, dirs, bam)
   names(sampleBams) <- sampleNames
   
   handleInputParams <- list(CLIP_reads=FALSE, fix_width=0, fix_point="start", norm=FALSE, useScore=FALSE,
                             outRle=FALSE, useSizeFactor=FALSE, genome=genome)
   
   comparisons <- list(c("RPRD1B_GFP_1", "GFP_1"), c("RPRD1B_GFP_2", "GFP_2"), c("DKO_1", "Scramble_1"), c("DKO_2", "Scramble_2"), c("DKO_3", "Scramble_1"), c("DKO_2", "DKO_1"), c("DKO_3", "DKO_1"), c("DKO_3", "DKO_2"), c("RPRD1A_GFP_2", "RPRD1A_GFP_1"), c("RPRD1B_GFP_2", "RPRD1B_GFP_1"), c("GFP_2", "GFP_1"), c("RPRD1A_GFP_1", "GFP_1"), c("RPRD1A_GFP_2", "GFP_2"),  c("Scramble_2", "Scramble_1"))
   
   
   for(i in 1:length(comparisons)){
      comparison <- comparisons[[i]]
      df <- threePrimeEnrichmentScore(treatBam=sampleBams[comparison[1]], controlBam=sampleBams[comparison[2]], txdb, stranded, handleInputParams, nc)
   }
   
   n_bins <- 20
   dataFiles_list <- list(RPRD1AGFP_rep2_vs_rep1_pooled=c("RPRD1A_GFP_2_vs_RPRD1A_GFP_1_3prime_enrichment_score.tab", "GFP_2_vs_GFP_1_3prime_enrichment_score.tab"),
                          RPRD1BGFP_rep2_vs_rep1_pooled =c("RPRD1B_GFP_2_vs_RPRD1B_GFP_1_3prime_enrichment_score.tab", "GFP_2_vs_GFP_1_3prime_enrichment_score.tab"),
                          DKOScramble_rep2_vs_rep1_pooled1=c("DKO_2_vs_DKO_1_3prime_enrichment_score.tab","Scramble_2_vs_Scramble_1_3prime_enrichment_score.tab"),
                          DKOScramble_rep2_vs_rep1_pooled2=c("DKO_3_vs_DKO_1_3prime_enrichment_score.tab","Scramble_2_vs_Scramble_1_3prime_enrichment_score.tab"),
                          DKOScramble_rep2_vs_rep1_pooled3=c("DKO_3_vs_DKO_2_3prime_enrichment_score.tab","Scramble_2_vs_Scramble_1_3prime_enrichment_score.tab"),
                          RPRD1A_vs_GFP_pooled=c("RPRD1A_GFP_1_vs_GFP_1_3prime_enrichment_score.tab", "RPRD1A_GFP_2_vs_GFP_2_3prime_enrichment_score.tab"),
                          RPRD1B_vs_GFP_pooled=c("RPRD1B_GFP_1_vs_GFP_1_3prime_enrichment_score.tab", "RPRD1B_GFP_2_vs_GFP_2_3prime_enrichment_score.tab"),
                          DKO_vs_Scramble_pooled1=c("DKO_1_vs_Scramble_1_3prime_enrichment_score.tab", "DKO_2_vs_Scramble_2_3prime_enrichment_score.tab"),
                          DKO_vs_Scramble_pooled2=c("DKO_3_vs_Scramble_1_3prime_enrichment_score.tab", "DKO_2_vs_Scramble_2_3prime_enrichment_score.tab"))
   for(dataName in names(dataFiles_list)){
      dataFiles <- dataFiles_list[[dataName]]
      plot_3prime_enrichment(datafm=NULL, dataFiles, n_bins, dataName, stranded)
   }
}

