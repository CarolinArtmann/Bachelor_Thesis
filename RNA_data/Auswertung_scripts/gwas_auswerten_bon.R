load("/storage/evolgen/data/tair10.rda")
load("~/RNA/genes_high_all.rda")
tair10_filtered <- tair10[which(tair10$Gene %in% genes_high_all),]
thresh <- 0.05/2500000

for(i in seq(1,100)){


if(i %in% which(genes_high_all %in% tair10$Gene)){
  results <- read.csv(paste("~/RNA/GWAS_perm/",i,"/results.csv",sep=""))
  RNA_results <- results[which(results$pval < thresh),]
  if(nrow(RNA_results) != 0){
 	 for(j in 1:nrow(RNA_results)){
   	 RNA_results$tair_match[j] <- paste(tair10_filtered$Gene[which(tair10_filtered$Chr == RNA_results[j,1] 
                                                       & tair10_filtered$Start-100000 < RNA_results[j,2] 
                                                       & tair10_filtered$Stop+100000 > RNA_results[j,2])], collapse = ", ")
 	 if(tair10_filtered$Gene[i] %in% strsplit(as.character(RNA_results$tair_match[j]), ", ")[[1]]){
   	 RNA_results$cis[j] <- "y"
 	 }
 	 else
   	 RNA_results$cis[j] <- "n"
 	 }
	}
  save(RNA_results, file=paste("~/RNA/RNA_bon/",i,".rda",sep=""))
}

}
RNA_bon <- mixedsort(list.files())
RNA_bon <- lapply(RNA_bon, load,.GlobalEnv)
