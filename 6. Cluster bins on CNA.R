# Clustering of bins based on CNAs

# Pvclust was previously used to compute a euclidean distance matrix and cluster with boostrapping. 
# Using data showing 1=aneuploid, 0=diploid, averaged across MR samples
load("~/Documents/PhD Barts/CNA AdVsCar/Predicting ITH/HPC/output/hiclust_diploid.aneu.rda")
hiclust <- hiclust_diploid.aneu

# Extract significant clusters into a dataframe
sig.hclust <- pvpick(hiclust, alpha=0.95, pv="au") 

# For each cluster extract the corresponding bins
list <- list()
for ( i in 1:length(sig.hclust$clusters) ) { 
  list[[i]] <- data.frame(bin = unlist(sig.hclust$clusters[[i]]), cluster = i)
}

# Create a dataframe with columns: bin | cluster
sig.hclust <- do.call('rbind', list)
sig.hclust$bin <- sub('.', '', sig.hclust$bin)

# Convert dataframe to numeric
sig.hclust <- data.frame(apply(sig.hclust, 2, function(x) as.numeric(as.character(x))))

# Add chromosome data based on bin number
sig.hclust$chr <- carcinoma$start.stop[match(sig.hclust$bin, carcinoma$start.stop$bin), 2]


# Create dataframe showing which chromosomes are within each cluster

# For each cluster count how many chromosomes are within, and identify chromosome
list1 <- list()
list2 <- list()
for ( i in 1:max(sig.hclust$cluster) ) {
  data <- sig.hclust[which(sig.hclust$cluster==i),]
  list1[[i]] <- paste(unique(data$chr), collapse = ',')
  list2[[i]] <- length(unique(data$chr))
}
cluster.info <- data.frame(cluster = 1:max(sig.hclust$cluster),
                           no.ofChr = do.call('rbind',list2),
                           chr = do.call('rbind',list1) )

# Add cluster group for facet
cluster.info$chrGroup <- as.numeric(as.character(cluster.info$chr))
cluster.info$chrGroup[is.na(cluster.info$chrGroup)] <- 23


# Add % of chromosome represented by each cluster

# For each cluster holding bins from one chromosome, claculate % of chromosome represented
list <- list()
for ( i in 1:length(sig.hclust$cluster) ) {
  chr <- as.numeric(cluster.info[which(cluster.info$cluster==i),]$chr)
  chrlength <- carcinoma$chr.end$end[chr] - carcinoma$chr.end$start[chr]
  number <- nrow(sig.hclust[which(sig.hclust$cluster==i),])
  list[[i]] <- number / chrlength
}

# Add % of chromosome into cluster.info dataframe
cluster.info <- cbind(cluster.info, data.frame(pcChrPerClus = unlist(list)))
cluster.info$pcChrPerClus <- cluster.info$pcChrPerClus*100


# Calculate the % aneuploidy per cluster per patient
list <- list()
l <- 1

# For a given patient
for ( k in 1:carcinoma$noPatients ) { 
  
  # For a given cluster, calculate the pc aneuploidy across the bins in the cluster
  for ( i in 1:max(sig.hclust$cluster) ) { 
    bins <- sig.hclust$bin[which(sig.hclust$cluster == i)]
    data <- matrices.list$diploid.aneu[bins, k]
    list[[l]] <- (sum(data[data==1])) / (length(data))
    l <- l + 1
  }
}

# Unlist into a dataframe
pcClusAneu <- data.frame(matrix(unlist(list), ncol = carcinoma$noPatients), stringsAsFactors = FALSE)
colnames(pcClusAneu) <- carcinoma$patients

# Add to cluster.info
cluster.info <- cbind(cluster.info, 
                      pcClusAneu = (apply(pcClusAneu, 1, mean))*100)


# Use % aneuploidy per cluster to predict diversity per patient

# Transpose to make patient per row, and a cluster per col
clustReg.in <- t(pcClusAneu)

# Add in ITH per patient
clustReg.in <- cbind(data.frame(ITH = carcinoma$ith[rep(seq_len(carcinoma$noPatients), each = 2), 5]), clustReg.in)

# To store output
clustReg.out <- data.frame(cluster=1:max(cluster.info$cluster), coeff=NA, pval=NA, adjR2=NA, Modelpval=NA)

# For a each cluster, consider if % aneuploidy can predict ITH
for ( k in 2:ncol(clustReg.in) ) { 
  data <- clustReg.in[,c(1,k)]
  data <- na.omit(data)
  colnames(data) <- c('ITH','predictor')
  
  # Regression can only be run when there are >1 factors present
  if ( length(unique(as.character(data$predictor))) == 1 ) { 
    next
  }
  reg <- lm(ITH ~ ., data)
  
  coeffs <- data.frame(t(summary(reg)$coefficients[,1])) 
  pvals <- data.frame(t(summary(reg)$coefficients[,4]))
  adjR2 <- round(summary(reg)$adj.r.squared,3)
  pval <- signif(lmp(reg),3)
  
  clustReg.out[(clustReg.out$cluster==(k-1)),c(2:5)] <- c(coeffs$predictor,pvals$predictor,adjR2,pval)
}

# Add col to indicate if aneuploidy in cluster is predictive (p<0.05)
clustReg.out$sig <- 'insig'
clustReg.out$sig[which(clustReg.out$Modelpval <= 0.05)] <- 'sig'

# Add FDR
clustReg.out$Bonferroni = p.adjust(clustReg.out$Modelpval, method = "BH")

# Designate significance of each cluster
clustReg.out$sigFDR <- 'insig'
clustReg.out$sigFDR[which(clustReg.out$BH <= 0.05)] <- 'sig'

# add coeffs and pvalx
cluster.info <- cbind(cluster.info, clustReg.out[,c(2:8)])
