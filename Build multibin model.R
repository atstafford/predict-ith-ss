# Multi-bin combinations to predict genome wide diversity
# Requires the 'cluster bins on CNA.R' to be run in order to calculate clusters

# First charcterise correlations between aberrations in individual bins and genome CNA diversity

# We will use 3 of the matrices: diploid.aneu (dip=0, aneu=1), diploid.gain (dip=0 gain=1), diploid.loss (dip=0 loss=1)
uniReg.in.list <- car.matrices[c(1,2,3)]

uniReg.in.list <- lapply(uniReg.in.list, function(x) {
  # Transpose to make a bin per column
  x <- t(x)
  
  # Convert to character and then a factor
  x <- data.frame(apply(x, 2, as.character), row.names = rownames(x), check.names = FALSE)
  x <- data.frame(lapply(x, factor, levels=c(0,1), labels=c('diploid','CNA')), row.names = rownames(x), check.names = FALSE)
  
  # # Add PIC.frac (stored in col 5 of ith df) to be predicted in first column
  x <- cbind(PIC.frac = carcinoma$ith[rep(seq_len(carcinoma$noPatients), each = 2), 5], x)
  x
})

# Run regression for each bin, and collect the coefficients for a gain, a loss, or noCNA at that bin

# Create list of 3 dataframes to store univariate regression output
uniReg.out.list <- list(data.frame(CNA='aneuploid', chr=carcinoma$start.stop$chr, bin=carcinoma$start.stop$bin, coeff=NA, pval=NA),
                        data.frame(CNA='gain', chr=carcinoma$start.stop$chr, bin=carcinoma$start.stop$bin, coeff=NA, pval=NA),
                        data.frame(CNA='loss', chr=carcinoma$start.stop$chr, bin=carcinoma$start.stop$bin, coeff=NA, pval=NA))
names(uniReg.out.list) = c('diploid.aneu', 'diploid.gain', 'diploid.loss')

# For each input matrix
for ( i in 1:length(uniReg.in.list) ) { 
  
  # Skipt first column as it holds ITH
  for ( k in 2:ncol(uniReg.in.list[[i]]) ) {
    
    # Creating a working dataframe holding ITH and the predictor bin
    data <- uniReg.in.list[[i]][,c(1,k)]
    data <- na.omit(data)
    colnames(data) <- c('ITH','predictor')
    
    # Set diploid as baseline
    data <- data %>%
      mutate(predictor = relevel(predictor, ref = 'diploid')) 
    
    # Regression can only be run when there are >1 factors present
    if ( length(unique(as.character(data$predictor))) == 1 ) { 
      next
    }
    
    # Run regression
    reg <- lm(ITH ~ ., data)
    
    # Extract the coefficient and pvalue
    coeffs <- data.frame(t(summary(reg)$coefficients[,1])) 
    names(coeffs) <- unlist(reg$xlevels)
    pvals <- data.frame(t(summary(reg)$coefficients[,4]))
    names(pvals) <- unlist(reg$xlevels)
    
    if ('CNA' %in% names(coeffs)) {
      uniReg.out.list[[i]][(uniReg.out.list[[i]]$bin==(k-1)),c(4:5)] <- c(coeffs$CNA,pvals$CNA)
    }
  }
}


# Add frequency of diploid, gain, loss

# For each regression output matrix
for ( i in 1:length(uniReg.out.list) ) {
  
  freq <- list()
  l <- 1
  
  # Pull the frequency of a diploid, gain, or loss occuring in at least one sample and divde by noPatient to get fraction
  for ( k in 1:nrow(uniReg.out.list[[i]]) ) {
    bin <- uniReg.out.list[[i]]$bin[k]
    status <- uniReg.out.list[[i]]$CNA[k]
    
    if ( status == 'aneuploid') {
      freq[[l]] <- (carcinoma$GaLo.clo$gain[bin] + carcinoma$GaLo.clo$loss[bin]) / carcinoma$noPatients
      l <- l + 1
    }
    else if ( status == 'gain') {
      freq[[l]] <- carcinoma$GaLo.clo$gain[bin] / carcinoma$noPatients
      l <- l + 1
    }
    else if ( status == 'loss') {
      freq[[l]] <- carcinoma$GaLo.clo$loss[bin] / carcinoma$noPatients
      l <- l + 1
    }
  }
  
  uniReg.out.list[[i]]$CNA.freq <- unlist(freq)
}


# Add adjusted p-value
uniReg.out.list <- lapply(uniReg.out.list, function(x) {
  x$sig <- 'insig'
  x$sig[which(x$pval<=0.05)] <- 'sig'
  x$pcorrect = p.adjust(x$pval, method = "BH")
  x$adjsig <- 'insig'
  x$adjsig[which(x$pcorrect<=0.05)] <- 'sig'
  x
})

# Plot coeffs and freq of loss/gain/diploid
text <- element_text(size=24, colour='black')
text.bold <- element_text(size=24, colour='black',face='bold')

# Store the 3 figures in a list
plot.list <- list()
l <- 1

for ( i in 1:length(uniReg.out.list) ) {
  # Set colours for: dark background, insig bars, sig bars, light background
  if ( i == 1 ) { col <- c("white",alpha("#009999", 0.2),"#CC9966"); labs <- c(1:18,'\n19','20','\n21','22'); tit <- 'Aneuploidy' } 
  if ( i == 2 ) { col <- c("white",alpha("#CC0033", 0.2),"#CC0033","#CC9966"); labs <- c(1:18,'\n19','20','\n21','22'); tit <- 'Gains only' } 
  if ( i == 3 ) { col <- c("white",alpha("#330099", 0.2),"#330099","#CC9966"); labs <- c(1:18,'\n19','20','\n21','22'); tit <- 'Losses only' }
  
  plot.list[[l]] <- ggplot(uniReg.out.list[[i]]) +
    geom_bar(aes(x=bin, y=coeff, fill=factor(adjsig)), stat="identity", width = 1) + 
    geom_rect(data=carcinoma$chr.end, aes(NULL,NULL,xmin=start,xmax=end,fill=col),ymin=-Inf,ymax=Inf,alpha=0.1) + 
    scale_fill_manual(values = col) + 
    
    geom_line(aes(x=bin, y=CNA.freq/3.33), group=1, color='#000033', linetype=1, size=0.7) + 
    geom_hline(yintercept = 0, size=0.3, color="black") + 
    
    ggtitle(tit) +
    scale_x_continuous(expand = c(0,0), name="chromosome", breaks=carcinoma$chr.mid, labels = labs) +
    scale_y_continuous(breaks = c(seq(-0.1,0.3,0.1)), limits = c(-0.14,0.34), 
                       sec.axis = sec_axis(trans = ~.*3.33)) +
    
    theme(plot.margin = unit(c(t=0,r=0.7,b=-0.2,l=0.7), "cm"),
          panel.background = element_blank(),
          plot.title = text.bold,
          
          axis.title = element_blank(),
          axis.text.x = text,
          axis.text.y = text,
          axis.line = element_line(colour = "black", size = 0.5),
          axis.ticks.length=unit(0.2, "cm"),
          legend.position = "none") 
  l <- l + 1
}

# Prepare common legend 
coeff.legend <- cowplot::get_legend(
  ggplot(uniReg.out.list$diploid.gain) +
    geom_bar(aes(x=bin, y=coeff, fill=factor(sig)), stat="identity") + 
    geom_line(aes(x=bin, y=500*CNA.freq, colour='line')) + 
    scale_fill_manual(values = c(alpha("#CC0033", 0.2),"#CC0033"), labels = c('p>0.05', 'p\u22640.05')) +
    labs(colour = NULL, fill = NULL) +
    scale_color_manual(values = '#000066', labels = 'Freq.') + 
    guides(fill = guide_legend(direction = 'horizontal', label.position = "left", label.hjust = 1),
           colour = guide_legend(direction = 'horizontal', label.position = "left", label.hjust = 4)) +
    theme(plot.margin=margin(0,0,0,0,"cm"),
          legend.position = "bottom",legend.direction="horizontal",
          legend.margin = margin(grid::unit(c(0,0,0,0),"cm")),
          legend.text = element_text(size = 24, margin = margin(l = 20, unit = "pt")),
          legend.key.height = grid::unit(0.8,"cm"),
          legend.key.width = grid::unit(0.8,"cm")) )

# Assemble plots
coeff.plot <- cowplot::plot_grid(as_ggplot(coeff.legend), plot.list[[1]], plot.list[[2]], plot.list[[3]], ncol = 1, align = "v", rel_heights = c(0.2,1,1,1))

coeff.plot <- annotate_figure(coeff.plot, 
                              left = text_grob('Coefficient', rot = 90, size = 24),
                              right = text_grob('Frequency of aberration', rot = 270, size = 24),
                              bottom = text_grob('Chromosome', size = 24))


# Perform stepwise selection to maximise Akaike information criteria, using one bin per cluster

# Bind the three dataframes holding the univariate regression data
candidate.bins <- do.call('rbind',uniReg.out.list[c(2:3)])

# Extract bins that are significant (p<0.05) before the FDR was applied as so few bins was sig after FDR
candidate.bins <- candidate.bins[which(candidate.bins$pval <= 0.05),]
candidate.bins$cluster <- sig.hclust[match(candidate.bins$bin, sig.hclust$bin), 2]
rownames(candidate.bins) <- 1:nrow(candidate.bins)

# See if there are any bins where multiple CNA statuses are sig
candidate.bins$bin[duplicated(candidate.bins$bin)] 

# Assign unique clusters to bins with no cluster
new.clust <- max(candidate.bins$cluster[is.finite(candidate.bins$cluster)]) + 1
for (i in 1:nrow(candidate.bins) ) {
  if ( is.na(candidate.bins$cluster[i]) ) {
    candidate.bins$cluster[i] <- new.clust
    new.clust <- new.clust + 1
  }
}


# Where multiple signiciant bins fall within the same cluster, select only the bin with the largets effect size (coefficient)

# Select one representive bin per cluster
list <- list()
for ( i in 1:max(candidate.bins$cluster) ) {
  data <- candidate.bins[which(candidate.bins$cluster == i),]
  list[[i]] <- data[which.max(abs(data$coeff)),]
}
representitive.bins <- do.call('rbind',list)


# Prepare input for multivariate regression

# The input matrix requires data on loss/gain/diploid, with a bin per column
multiReg.in <- t(car.raw[,-c(1:3)])

# Convert matrix data to character then factor
multiReg.in <- data.frame(apply(multiReg.in, 2, as.character), check.names = FALSE)
multiReg.in <- data.frame(lapply(multiReg.in, factor, levels=c(2,1,3), labels=c('diploid','loss','gain')), check.names = FALSE)

# Keep only the representitive bins (columns)
multiReg.in <- multiReg.in[,c(representitive.bins$bin)]

# Add ITH to be predicted in first column
multiReg.in <- cbind(ITH = carcinoma$ith[rep(seq_len(carcinoma$noPatients), each = 2), 5], multiReg.in)
rownames(multiReg.in) <- carcinoma$samples


# Perform stepAIC

# Create the full model with all representitive bins
full.model <- lm(ITH ~., data = multiReg.in)

# Perform forwards and backwards selection
step.model <- stepAIC(full.model, direction = "both", trace = FALSE)


# Extract step.model info

# Pull the bin names which step.AIC selected
predictors <- data.frame(bin = rownames(summary(step.model)$coefficients))

# Edit the names to keep only the bin number, and remove duplicates
predictors$bin <- as.numeric(str_extract_all(predictors$bin, "[0-9]+"))

# Add in chr and start/stop codon
predictors <- merge(unique(predictors), carcinoma$start.stop)

# Add cluster info for each bin
predictors$cluster <- candidate.bins[match(predictors$bin, candidate.bins$bin), 10] 

# Pull the adjusted R2 and pval of the step.model
adjR2 <- round(summary(step.model)$adj.r.squared,3)
pval <- signif(lmp(step.model),3)


# Predict using training (carcinoma)

# List actual patient ITH, repeat each value twice as there are two samples
actualIth <- carcinoma$ith[rep(seq_len(carcinoma$noPatients), each = 2), 5]

# Predict ITH using step.model
predictedIth <- predict(step.model, multiReg.in) 

# Enter actual and predicted ITH into a comparative dataframe
comparison <- as.data.frame(cbind(actualIth, predictedIth))
comparison$patient <- sub('\\..*', '', rownames(comparison))

# Add two labels for plot
label1 <- paste('R[adj]^2 ==', adjR2)
label2 <- paste('p.value ==', pval)

# Plot
text <- element_text(size=24, colour='black')
text.bold <- element_text(size=24, colour='black',face='bold')

plot <- ggplot(data = comparison, aes(x = actualIth, y = predictedIth)) +
  geom_line(aes(group = patient), size=0.2, colour = "#003366") +
  geom_point(fill = "#003366", colour = "#003366", size = 5) +
  geom_line(aes(x = actualIth, y = actualIth), linetype = "dashed") +
  xlab("Patient CNA diversity") +
  ylab("Patient CNA diversity (predicted)") +
  xlim(0,0.5) +
  ylim(0,0.5) +
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = text,
        axis.text = text,
        legend.position = "none") +
  annotate('text', x = c(0), y = c(0.500,0.475), label = c(label1, label2), parse=TRUE, size = 8, hjust = 0, vjust = 1)


# Model diagnostics

# The residuals vs fitted plot 
# Residuals are evenly spread around the 0 line. This suggests that the assumption that the relationship is linear is reasonable.
# Linearity holds wells, as red line is close to dashed line. This suggests that the variances of the error terms are equal.
# There are three outliers. This suggests that the assumption of homoscedasticity is reasonable.
# Theres no 'cone'.
plot(step.model,1)

# Scale location: Homogeneity of variance
# This plot test the linear regression assumption of equal variance (homoscedasticity)
# Residuals should have equal variance along the regression line. 
# This plot shows if residuals are spread equally along the ranges of predictors
# Good if you see a horizontal line with equally spread points. 
# The residuals appears to have a non-constant variance, but the line is horizontal
# 'Not plotting observations with leverage one'
plot(step.model,3)

# Normal Q-Q plot: Normality of residuals
# The normal probability plot of residuals should approximately follow a straight line.
# All the points fall approximately along this reference line, so we can assume normality.
plot(step.model,2)

# Residuals vs leverage: Outliers and high levarage points
# Can be used to find influential cases in the dataset
# An influential case may or may not be an outlier and the purpose of this chart is to identify cases that have high influence in the model. 
# Outliers will tend to exert leverage and therefore influence on the model.
# An influential case will appear in the top right or bottom left of the chart inside a red line which marks Cook’s Distance
# There are none
plot(step.model,5)


## What is the difference in predicted ITH between testing sample 1 and sample 2

# Calculate the difference in ITH between the two samples
predDiff <- list()
i <- 1
while ( i < nrow(comparison) ) {
  if ( comparison$predictedIth[i] >= comparison$predictedIth[i+1] ) {
    predDiff[[i]] <- comparison$predictedIth[i] - comparison$predictedIth[i+1]
    i <- i + 2
  }
  else if ( comparison$predictedIth[i] < comparison$predictedIth[i+1] ) {
    predDiff[[i]] <- comparison$predictedIth[i+1] - comparison$predictedIth[i]
    i <- i + 2
  }
}

predDiff <- do.call('rbind',predDiff)

mean(predDiff)
max(predDiff)
min(predDiff)

# Compare the error of sample 1 and 2

# Find the difference between predicted and actual
comparison$error <- NA
for ( i in 1:nrow(comparison) ) {
  comparison$error[i] <- comparison$actualIth[i] - comparison$predictedIth[i]
}

# Take the last number from the sample ID to identify is as 1 or 2
comparison$sample <- str_sub(carcinoma$samples,5,5)

# To store the error values for each patient
error <- data.frame(carcinoma$patients, sample1 = NA, sample2 = NA)

# To reduce lines crossing over, sample1 will have the highest error, and sample2 the lowest
for ( i in 1:carcinoma$noPatients ) {
  patient <- error$carcinoma.patients[i]
  data <- comparison[which(comparison$patient==patient),]
  error$sample1[i] <- max(data$error)
  error$sample2[i] <- min(data$error)
}

# Plot
plot <- ggplot(data = error) +
  geom_segment(aes(x=0, xend=1, y=sample1, yend=sample2), size=0.4, color='#003366') +
  geom_point(aes(x=rep.int(0,nrow(error)), y=sample1), color='#003366', size=5) +
  geom_point(aes(x=rep.int(1,nrow(error)), y=sample2), color='#003366', size=5) +
  scale_x_continuous(name = NULL, breaks = c(0,1), labels = c('Sample 1','Sample 2')) +
  scale_y_continuous(name = 'Actual - Predicted', limits = c(-0.13, 0.3)) +
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = text,
        axis.text = text,
        legend.position = "none")