## for plotting pooled 3primeErichment score files
scriptd <- dirname(this.path::this.path())
source(file.path(scriptd, "3primeEnrichment_score.R"))

if(0){
args <- commandArgs(trailingOnly = T)
wd <- args[1] # where the output should be
setwd(wd)

dataFiles <- unlist(strsplit(args[2], split=",", fixed=TRUE)) # score files, comma separated
dataName <- args[3] # "DKO_pooled_vs_Scramble_pooled"
n_bins <- as.integer(args[4])
stranded <- as.logical(args[5])

plot_3prime_enrichment(datafm=NULL, dataFiles, n_bins, dataName, stranded)
}

if(1){
  
   wd <-  "C:/GREENBLATT/Zuyao/RNAseq/Nonstranded"
   setwd(wd)
   
   n_bins <- 20
   stranded <- FALSE
   
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