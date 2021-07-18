library(gtools)
load("~/herit/GWAS/tair10.rda")
load("~/herit/GWAS/genes_high_all.rda")
tair10_filtered <- tair10[which(tair10$Gene %in% genes_high_all),]
thresh <- 0.05/2500000
i <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

cat(i)


if(i %in% which(genes_high_all %in% tair10$Gene)){
  results <- read.csv(paste("~/RNA/RNA_boxcox/",i,results.csv,sep="")
  RNA_results <- results[which(results$pval < thresh),]
  if(nrow(RNA_results) != 0){
 	 for(j in 1:nrow(RNA_results)){
   	 RNA_results$tair_match[j] <- paste(tair10_filtered$Gene[which(tair10_filtered$Chr == RNA_results[j,1] 
                                                       & tair10_filtered$Start-50000 < RNA_results[j,2] 
                                                       & tair10_filtered$Stop+50000 > RNA_results[j,2])], collapse = ", ")
 	 if(tair10_filtered$Gene[i] %in% strsplit(as.character(RNA_results$tair_match[j]), ", ")[[1]]){
   	 RNA_results$cis[j] <- "y"
 	 }
 	 else
   	 RNA_results$cis[j] <- "n"
 	 }
	}
  save(RNA_results, file=paste("~/RNA/RNA_box_bon/",i,".rda",sep=""))
}

