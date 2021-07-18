library(gtools)
load("~/herit/GWAS/tair10.rda")
load("~/herit/GWAS/genes_high_all.rda")
tair10_filtered <- tair10[which(tair10$Gene %in% genes_high_all),]
i <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))


if(i %in% which(genes_high_all %in% tair10$Gene)){
	thresh <- read.csv(paste("~/RNA/RNA/", i, "/perm_results.csv",sep=""))[5,3]
  results <- read.csv(paste("~/RNA/RNA/",i,"/results.csv",sep="")
  RNA_results_perm <- results[which(results$pval < thresh),]
  if(nrow(RNA_results_perm) != 0){
 	 for(j in 1:nrow(RNA_results_perm)){
   	 RNA_results_perm$tair_match[j] <- paste(tair10_filtered$Gene[which(tair10_filtered$Chr == RNA_results_perm[j,1] 
                                                       & tair10_filtered$Start-50000 < RNA_results_perm[j,2] 
                                                       & tair10_filtered$Stop+50000 > RNA_results_perm[j,2])], collapse = ", ")
 	 if(tair10_filtered$Gene[i] %in% strsplit(as.character(RNA_results_perm$tair_match[j]), ", ")[[1]]){
   	 RNA_results_perm$cis[j] <- "y"
 	 }
 	 else
   	 RNA_results_perm$cis[j] <- "n"
 	 }
	}
  save(RNA_results_perm, file=paste("~/RNA/RNA_results/RNA_results_perm",i,".rda",sep=""))
}

