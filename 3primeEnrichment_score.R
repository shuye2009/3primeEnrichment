library(DescTools)
library(plyranges)
library(GenomicPlot)
library(BiocFileCache)
library(foreach)
library(doParallel)
library(GenomeInfoDb)
library(fields)
library(dplyr)

threePrimeEnrichmentScore <- function(treatBam, controlBam, gtfFile, stranded=FALSE, handleInputParams=NULL, nc=20){
   txdb <- makeTxDbFromGFF(gtfFile)
   
   print("obtaining gene coordinates")
   gene <- get_genomic_feature_coordinates(txdb, featureName="gene", longest=TRUE, protein_coding=TRUE)
   gene_gr <- gene$GRanges
   gene_gr <- gene_gr[width(gene_gr)>100]
   print(head(gene_gr))
   
   print("obtaining transcripts lendth")
   tx_lens_df <- transcriptLengths(txdb)
   tx_lens <- tx_lens_df$tx_len
   names(tx_lens) <- tx_lens_df$gene_id
   
   system.time(inputs <- handle_input(c(treatBam, controlBam), handleInputParams))
   dataName <- paste(c(names(treatBam), names(controlBam)), collapse = "_vs_")
   print("Reading in data succeeded")
   
   if(!stranded){
      print(system.time(sample_cvg <- coverage(inputs[[1]]$query)))
      print(system.time(control_cvg <- coverage(inputs[[2]]$query)))
   
      print("starting calculate scores")
      scores <- NULL
      print(system.time({
         cl <- start_parallel(nc)
         #registerDoParallel(cl)
         print("exporting global environment variables")
         clusterExport(cl, c("seqnames", "runValue", "expand_gene", "calculate_score", "AUC"), envir=.GlobalEnv)
         print("exporting local environment variables")
         clusterExport(cl, c("gene_gr", "tx_lens", "sample_cvg", "control_cvg"), envir=environment())
         print("variables exported")
         if(0){
         scores <- foreach(i=1:length(gene_gr), .combine="rbind") %dopar% {
            library(GenomicRanges)
            asample_cvg <- sample_cvg[[seqnames(gene_gr[i])]]
            acontrol_cvg <- control_cvg[[seqnames(gene_gr[i])]]
            len <- tx_lens[names(gene_gr[i])]
            score <- calculate_score(gene_gr[i], len, asample_cvg, acontrol_cvg, plt=FALSE)
            score
         }
         }
         scores <- parLapply(cl, seq_along(gene_gr), function(i){
            x <- gene_gr[i]
            asample_cvg <- sample_cvg[[seqnames(x)]]
            acontrol_cvg <- control_cvg[[seqnames(x)]]
            len <- tx_lens[names(x)]
            score <- calculate_score(x, len, asample_cvg, acontrol_cvg, plt=FALSE)
            score
            
         })
         stop_parallel(cl)
         scores <- bind_rows(scores)
      }))
      
      df <- data.frame(Name=names(gene_gr), Length=width(gene_gr))
      df <- cbind(df, scores)
      colnames(df) <- c("Name", "length", "Tx_length", "Score", "Count", "Count_ctr", "Depth")
      
   }else{
      sample_grl <- split(inputs[[1]]$query, strand(inputs[[1]]$query))
      control_grl <- split(inputs[[2]]$query, strand(inputs[[2]]$query))
      system.time(sample_cvg_pos <- coverage(sample_grl[["+"]]))
      system.time(sample_cvg_neg <- coverage(sample_grl[["-"]]))
      system.time(control_cvg_pos <- coverage(control_grl[["+"]]))
      system.time(control_cvg_neg <- coverage(control_grl[["-"]]))
      
      scores <- NULL
      
      print(system.time({
         cl <- start_parallel(nc)
         registerDoParallel(cl)
         
         print("exporting global environment variables")
         clusterExport(cl, c("seqnames", "runValue", "expand_gene", "calculate_score", "AUC"), envir=.GlobalEnv)
         print("exporting local environment variables")
         clusterExport(cl, c("gene_gr", "tx_lens", "sample_cvg_pos", "control_cvg_pos", "sample_cvg_neg", "control_cvg_neg"), envir=environment())
         print("variables exported")
         
         scores <- foreach(i=1:length(gene_gr), .combine="rbind") %dopar% {
            library(GenomicRanges)
            asample_cvg_pos <- sample_cvg_pos[[seqnames(gene_gr[i])]]
            acontrol_cvg_pos <- control_cvg_pos[[seqnames(gene_gr[i])]]
            asample_cvg_neg <- sample_cvg_neg[[seqnames(gene_gr[i])]]
            acontrol_cvg_neg <- control_cvg_neg[[seqnames(gene_gr[i])]]
            
            len <- tx_lens[names(gene_gr[i])]
            
            if(runValue(strand(gene_gr[i])) == "+"){
               sense_scores <- calculate_score(gene_gr[i], len, asample_cvg_pos, acontrol_cvg_pos, plt=FALSE)
               antisense_scores <- calculate_score(gene_gr[i], len, asample_cvg_neg, acontrol_cvg_neg, plt=FALSE, countPerBase=1)
            }else{
               antisense_scores <- calculate_score(gene_gr[i], len, asample_cvg_pos, acontrol_cvg_pos, plt=FALSE, countPerBase=1)
               sense_scores <- calculate_score(gene_gr[i], len, asample_cvg_neg, acontrol_cvg_neg, plt=FALSE)
            }
            c(sense_scores, antisense_scores)
         }
         stop_parallel(cl)
      }))
      df <- data.frame(Name=names(gene_gr), Length=width(gene_gr))
      df <- cbind(df, scores)
      colnames(df) <- c("Name", "length", "Tx_length", "Score", "Count", "Count_ctr", "Depth", "Antisense_tx_length", "Antisense_score", "Antisense_Count", "Antisense_count_ctr", "Antisense_depth")
      
   }
   write.table(df, paste0(dataName, "_3prime_enrichment_score.tab"), sep="\t", row.names=FALSE, quote=FALSE)
   plot_3prime_enrichment(datafm=df, dataFiles=NULL, n_bins=20, dataName=dataName, stranded=stranded)
   return(df)
}

expand_gene <- function(gr){
   s <- runValue(strand(gr))
   chr <- runValue(seqnames(gr))
   st <- start(gr)
   en <- end(gr)
   
   if(s == "+"){
      ranges <- IRanges(c(seq(en, st, -1), st), en)
      rlen <- length(ranges)
      grx <- GRanges(seqnames=Rle(chr, rlen),
                     ranges=ranges,
                     strand=rep(s, rlen))
   }else if(s == "-"){
      ranges <- IRanges(st, c(seq(st, en, 1), en))
      rlen <- length(ranges)
      grx <- GRanges(seqnames=Rle(chr, rlen),
                     ranges=ranges,
                     strand=rep(s, rlen))
   }else stop("strand is NOT accepted")
   return(grx)
}

calculate_score <- function(agene_gr, len, asample_cvg, acontrol_cvg, plt=FALSE, countPerBase=2){
   grx <- expand_gene(agene_gr)
   names(len) <- NULL
  
   sample_v <- Views(asample_cvg, start(grx), end(grx))
   counts <- viewSums(sample_v)
   control_v <- Views(acontrol_cvg, start(grx), end(grx))
   ctr_counts <- viewSums(control_v)
      
   enrichment_score <- 0
   depth <- max(ctr_counts)/len
   if(depth > countPerBase){
      
      suppressWarnings(auc <- AUC(ctr_counts, counts, method="trapezoid"))
      triangle <- as.double((counts[length(counts)]))*(ctr_counts[length(ctr_counts)])/2
      enrichment_score <- round(auc/triangle, digits=3)
      if(plt){
         plot(ctr_counts, counts, type="l")
         abline(0, (counts[length(counts)])/(ctr_counts[length(ctr_counts)]))
         text((ctr_counts[length(ctr_counts)])/10, (counts[length(counts)])/1.1, enrichment_score)
      }
   }
   
   return(c(len, enrichment_score, max(counts), max(ctr_counts), depth))
}


plot_3prime_enrichment <- function(datafm=NULL, dataFiles=NULL, n_bins=20, dataName="Test", stranded=FALSE){
   if(!is.null(dataFiles)){
      print(dataFiles)
      datafm <- lapply(dataFiles, read.delim, header=TRUE)
      datafm <- data.frame(bind_rows(datafm))
   }
   datafm <- datafm %>%
      filter(!is.na(Score)) %>%
      filter(Score > 0)
   data_sorted <- datafm[order(datafm[,2]),]
   head(data_sorted)
   size <- dim(data_sorted)[1]
   step = ceiling(size/(n_bins+1))
   bin <- seq(1, size, step)
   breaks <- data_sorted[bin,2]
   length(breaks)
   average_bin <- stats.bin(data_sorted[,2], data_sorted[,4], breaks=breaks)
   mean_score <- average_bin$stats["mean",]
   length(mean_score)
   
   pdf(paste0(dataName, "_3prime_enrichment.pdf"))
   plot(data_sorted[,2],data_sorted[,4], log="x", type="p",  pch=20, xlab="Gene Length",ylab="3 prime enrichment score", main=dataName,col="black")
   breaks <-breaks[2:length(breaks)] - step/2
   mean_score <- mean_score[1:length(mean_score)]
   lines(breaks, mean_score, type="b", col="red")
   abline(h=1, col="cyan")
   r<-cor(breaks, mean_score)
   rsq <- r*r
   r_2 = sprintf("%.3f", rsq)
   usr <- par( "usr" )
   mtext(paste("R-squared = ",r_2, sep=""), side=3, col="blue" )
   
   if(stranded){
      datafm <- datafm %>%
         filter(!is.na(Antisense_score)) %>%
         filter(Antisense_score > 0)
      data_sorted <- datafm[order(datafm[,2]),]
      head(data_sorted)
      size <- dim(data_sorted)[1]
      step = ceiling(size/(n_bins+1))
      bin <- seq(1, size, step)
      breaks <- data_sorted[bin,2]
      length(breaks)
      average_bin <- stats.bin(data_sorted[,2], data_sorted[,9], breaks=breaks)
      mean_score <- average_bin$stats["mean",]
      length(mean_score)
      plot(data_sorted[,2],data_sorted[,9], log="x", type="p",  pch=20, xlab="Gene Length",ylab="3 prime enrichment score", main=paste0(dataName, "_antisense"),col="black")
      breaks <-breaks[2:length(breaks)] - step/2
      mean_score <- mean_score[1:length(mean_score)]
      lines(breaks, mean_score, type="b", col="red")
      abline(h=1, col="cyan")
      r<-cor(breaks, mean_score)
      rsq <- r*r
      r_2 = sprintf("%.3f", rsq)
      usr <- par( "usr" )
      mtext(paste("R-squared = ",r_2, sep=""), side=3, col="blue" )
   }
   dev.off()
}
