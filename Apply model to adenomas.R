# Application of the step.model to adenomas
# We have 19 adenomas from the orinigal cohort and an additional 5 from the validation cohort

# Prepare input for predicting adenomas in original cohort

# Requires bin per col
multiReg.in <- t(ad.raw[,-c(1:3)])

# Change to factor
multiReg.in <- data.frame(apply(multiReg.in, 2, as.character), check.names = FALSE)
multiReg.in <- data.frame(lapply(multiReg.in, factor, levels=c(2,1,3), labels=c('diploid','loss','gain')), check.names = FALSE)

# Keep only the representitive bins
multiReg.in <- multiReg.in[,c(representitive.bins$bin)]

# Add ITH to be predicted in first column
PIC.frac <- data.frame(PIC.frac = adenoma$ith$PIC.frac)
multiReg.in <- cbind(ITH = adenoma$ith[rep(seq_len(adenoma$noPatients), each = 2), 5], multiReg.in)
rownames(multiReg.in) <- colnames(ad.raw[-c(1:3)])

# Predict
actualIth <- adenoma$ith[rep(seq_len(adenoma$noPatients), each = 2), 5]
predictedIth <- predict(step.model, multiReg.in) 
comparison_t <- as.data.frame(cbind(actualIth, predictedIth))
comparison_t$patient <- sub('\\..*', '', rownames(comparison_t))

# Load the 10 test validation set adenoma patients
test_ad1 <- read.table("~/Documents//PhD Barts/CNA AdVsCar/Data/testSet_ad/Polyp.08.WGS.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_ad2 <- read.table("~/Documents//PhD Barts/CNA AdVsCar/Data/testSet_ad/Polyp.02.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_ad3 <- read.table("~/Documents//PhD Barts/CNA AdVsCar/Data/testSet_ad/Polyp.05.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_ad4 <- read.table("~/Documents//PhD Barts/CNA AdVsCar/Data/testSet_ad/Polyp.09.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_ad5 <- read.table("~/Documents//PhD Barts/CNA AdVsCar/Data/testSet_ad/Polyp.03.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")

testraw.list <- list(test_ad1,test_ad2,test_ad3,test_ad4,test_ad5)
rm(test_ad1,test_ad2,test_ad3,test_ad4,test_ad5)

# Convert rawdata to numeric
testraw.list <- lapply(testraw.list, function(x) {
  x <- apply(x, 2, function(x) as.numeric(as.character(x)))
  x
})

x <- testraw.list[[1]]

# Generate absolute copy number (minor + major)
testcn.list <- testraw.list
testcn.list <- lapply(testcn.list, function(x) {
  x <- x[,c(1,2,4)]
  colnames(x) <- c('chr','start','stop')
  x <- as.data.frame(x)
  x
})

for ( j in 1:length(testraw.list)) {
  
  i <- 5
  l <- 1
  abscn <- list()
  
  while ( i < ncol(testraw.list[[j]]) ) {
    for ( k in 1:nrow(testraw.list[[j]]) ) {
      abscn[[l]] <- sum(testraw.list[[j]][k,i], testraw.list[[j]][k,i+1])
      l <- l + 1
    }
    i <- i + 2
  }
  
  testcn.list[[j]] <- cbind(testcn.list[[j]], matrix(unlist(abscn), nrow = nrow(testraw.list[[j]])) )
}

# Assess ploidy and recentre if average CN is 2.8 or more
testcn.list <- ploidyRecentre(testcn.list)

# Bin cn data to match training dataset
testcnBinned.list <- alignBins(bins = adenoma$start.stop, cn.list = testcn.list)

# If absolute copy number >=3, its a gain, if <=2 its a loss
testcnBinned.list <- lapply(testcnBinned.list, function(x) {
  x[,-c(1:4)] <- ifelse(x[,-c(1:4)] >= 3, 3, ifelse(x < 2, 1, 2))
  x
})

# Rename columns with sampleID and combine all patients into one dataframe
sample.count <- list()
for ( i in 1:length(testcnBinned.list) ) {
  for ( k in 5:ncol(testcnBinned.list[[i]]) ) {
    colnames(testcnBinned.list[[i]])[k] <- paste(i, k-4, sep = '.')
    sample.count[[i]] <- ncol(testcnBinned.list[[i]]) - 4 # to be used later
  }
}

testcnBinned <- testcnBinned.list %>% purrr::reduce(full_join, by = c("chr","start","stop","bin"))

# Calculate PIC.frac using all samples for each patient
pic <- list()
for ( i in 1:length(testcnBinned.list)) {
  # Set working data as the cols holding the samples
  wd <- testcnBinned.list[[i]][,-c(1:4)]
  
  # Record the number of samples
  upto <- ncol(wd)
  
  # Use PIC function on wd
  pic[[i]] <- PIC(wd, upto, c(1:upto))
  pic[[i]] <- na.omit(pic[[i]])
  
  # Define the max possible diversity given the number of sample, as maxPIC*number of bins
  max.ITH <- max.pics[[upto]] * length(pic[[i]])
  
  # Store as a dataframe
  pic[[i]] <- data.frame(patient = i, PIC.frac = sum(pic[[i]], na.rm = TRUE)/max.ITH,  samplesUsed = 'max')
}

testPIC.frac <- do.call('rbind',pic)


# Prepare input for prediction

# Transpose as input table requires cols=bins and convert to character
test.input <- data.frame(t(testcnBinned[,-c(1:4)]), check.names = FALSE)
test.input <- ifelse(test.input == 3, "gain", ifelse(test.input == 1, "loss", "diploid"))

# Convert to dataframe and change to factor
test.input <- data.frame(test.input, check.names = FALSE)
test.input <- data.frame(lapply(test.input, as.factor), row.names = rownames(test.input), check.names = FALSE)

# Use the model to predict diversity in test dataset
actualith.test <- data.frame(lapply(data.frame(testPIC.frac$PIC.frac), rep, sample.count),  row.names = rownames(test.input))
predictedith.test <- data.frame(predict(step.model, test.input))
comparison_v <- as.data.frame(cbind(actualith.test, predictedith.test)) 
comparison_v <- na.omit(comparison_v)
colnames(comparison_v) <- c("actualIth","predictedIth")
comparison_v$patient <- sub('\\..*', '', rownames(comparison_v))

# Combine results from the two cohorts
comparison_t$group <- 'training'
comparison_v$group <- 'validation'
x <- rbind(comparison_t, comparison_v)

# Plot
text <- element_text(size=24, colour='black')
text.bold <- element_text(size=24, colour='black',face='bold')

plot <- ggplot(data = x, aes(x = actualIth, y = predictedIth)) +
  geom_line(aes(group = patient), size=0.2, colour = "#003366") +
  geom_point(aes(color=factor(group)), size = 5) +
  geom_line(aes(x = actualIth, y = actualIth), linetype = "dashed") +
  scale_color_manual(values = c('#003366','#336699')) + 
  xlab("Patient CNA diversity") +
  ylab("Patient CNA diversity (predicted)") +
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = text,
        axis.text = text,
        legend.position = "none") +
  annotate('text', x = c(0), y = c(0.500,0.475), label = c(label1, label2), parse=TRUE, size = 8, hjust = 0, vjust = 1)