RNA_boxcox_perm_graph <- list()


# cis_hits
for(i in seq(1,100)){
  RNA_boxcox_perm_graph[[i]] <- data.frame(chr = rep(0), pos=rep(0))
  if(nrow(RNA_boxcox_perm[[i]]) == 0)
    RNA_boxcox_perm_graph[[i]][1,] <- c(NA,NA)
  else {cis <- RNA_boxcox_perm[[i]][,8] 
  if("y" %in% cis)
    RNA_boxcox_perm_graph[[i]][1,] <- RNA_boxcox_perm[[i]][which(cis == "y")[1], c(1,2)]
  else 
    RNA_boxcox_perm_graph[[i]][1,] <- c(NA,NA)
  }
}

#trans_hits

for(i in seq(1,100)){
  if(nrow(RNA_boxcox_perm[[i]]) == 0)
    RNA_boxcox_perm_graph[[i]][1,] <- c(NA,NA)
  else{
    trans <- subset(RNA_boxcox_perm[[i]], RNA_boxcox_perm[[i]][,8] == "n")
    if(nrow(trans) == 0)
       RNA_boxcox_perm_graph[[i]] <- rbind(RNA_boxcox_perm_graph[[i]], c(NA,NA))
    else RNA_boxcox_perm_graph[[i]] <- rbind(RNA_boxcox_perm_graph[[i]], trans[1,c(1,2)])
    if(nrow(trans) > 1){
      for(j in seq(1, nrow(trans)-1)){
        if(trans$chr[j] == trans$chr[j+1] & trans$pos[j+1] < trans$pos[j] + 100000)
          test <- 1
        else RNA_boxcox_perm_graph[[i]] <- rbind(RNA_boxcox_perm_graph[[i]], trans[j+1,c(1,2)]) 
      }
    }
  }
}

###1 ist cis

for(i in seq(1,100)){
  RNA_boxcox_perm_graph[[i]] <-  cbind(RNA_boxcox_perm_graph[[i]], tair = rep(tair10_graph$mean[i]))
  RNA_boxcox_perm_graph[[i]] <-  cbind(RNA_boxcox_perm_graph[[i]], cistrans = rep(0))
  RNA_boxcox_perm_graph[[i]][1,4] <- 1
}

RNA_boxcox_perm_graph_df <- do.call(rbind.data.frame, RNA_boxcox_perm_graph)

RNA_boxcox_perm_graph_df <- na.omit(RNA_boxcox_perm_graph_df)

RNA_boxcox_perm_graph_df$chr <- factor(RNA_boxcox_perm_graph_df$chr, levels = c(5,4,3,2,1))

