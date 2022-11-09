## for plotting pooled 3primeErichment score files
scriptd <- dirname(this.path::this.path())
source(file.path(scriptd, "3primeEnrichment_score.R"))

args <- commandArgs(trailingOnly = T)
wd <- args[1] # where the output should be
setwd(wd)

dataFiles <- unlist(strsplit(args[2], split=",", fixed=TRUE)) # score files, comma separated
dataName <- args[3] # "DKO_pooled_vs_Scramble_pooled"
n_bins <- as.integer(args[4])
stranded <- as.logical(args[5])


plot_3prime_enrichment(datafm=NULL, dataFiles, n_bins, dataName, stranded)