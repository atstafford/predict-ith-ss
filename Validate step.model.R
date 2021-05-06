# Validate step.model
# Requires the 'Buildmultibin model.R' to be run in order to create step.model

# Prepare validation dataset
# In order to use upcoming functions, data must be in format of a list with a dataframe per patient, with columns: 
# chr | start | stop | sample1 | sample2...

# Read in the datasets
test_car1 <- read.table("~/Documents//PhD Barts/CNA AdVsCar/Data/testSet_car/Set.01.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car2 <- read.table("~/Documents//PhD Barts/CNA AdVsCar/Data/testSet_car/Set.02.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car3 <- read.table("~/Documents//PhD Barts/CNA AdVsCar/Data/testSet_car/Set.03.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car4 <- read.table("~/Documents//PhD Barts/CNA AdVsCar/Data/testSet_car/Set.04.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car5 <- read.table("~/Documents//PhD Barts/CNA AdVsCar/Data/testSet_car/Set.05.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car6 <- read.table("~/Documents//PhD Barts/CNA AdVsCar/Data/testSet_car/Set.06.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car7 <- read.table("~/Documents//PhD Barts/CNA AdVsCar/Data/testSet_car/Set.07.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car8 <- read.table("~/Documents//PhD Barts/CNA AdVsCar/Data/testSet_car/Set.08.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car9p <- read.table("~/Documents//PhD Barts/CNA AdVsCar/Data/testSet_car/Set.09.Proximal.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car9d <- read.table("~/Documents//PhD Barts/CNA AdVsCar/Data/testSet_car/Set.09.Distal.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")
test_car10 <- read.table("~/Documents//PhD Barts/CNA AdVsCar/Data/testSet_car/Set.10.penalty0.95.baf.gt.txt", header=T, skip=0, sep="\t")

# Load patient datasets into a list
testraw.list <- list(test_car1,test_car2,test_car3,test_car4,test_car5,test_car6,test_car7,test_car8,test_car9p,test_car9d,test_car10)
rm(test_car1,test_car2,test_car3,test_car4,test_car5,test_car6,test_car7,test_car8,test_car9p,test_car9d,test_car10)

# Convert data to numeric
testraw.list <- lapply(testraw.list, function(x) {
  x <- apply(x, 2, function(x) as.numeric(as.character(x)))
  x
})

# For each patient, create dataframe to hold absolute copy number per chr/start/stop
testcn.list <- lapply(testraw.list, function(x) {
  x <- x[,c(1,2,4)]
  colnames(x) <- c('chr','start','stop')
  x <- as.data.frame(x)
  x
})

# Generate absolute copy number (minor + major)
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

# Bin CN data to match training dataset
testcnBinned.list <- alignBins(bins = carcinoma$start.stop, cn.list = testcn.list)

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

# Calculate PIC

# Define the maximum pic score depending on the number of samples (up to 13)
max.pic1 <- PIC <- 1 - ( (1/1)^2 + (0/1)^2 + (0/1)^2 )
max.pic2 <- PIC <- 1 - ( (1/2)^2 + (1/2)^2 + (0/2)^2 )
max.pic3 <- PIC <- 1 - ( (1/3)^2 + (1/3)^2 + (1/3)^2 )
max.pic4 <- PIC <- 1 - ( (1/4)^2 + (1/4)^2 + (2/4)^2 )
max.pic5 <- PIC <- 1 - ( (1/5)^2 + (2/5)^2 + (2/5)^2 )
max.pic6 <- PIC <- 1 - ( (2/6)^2 + (2/6)^2 + (2/6)^2 )
max.pic7 <- PIC <- 1 - ( (2/7)^2 + (2/7)^2 + (3/7)^2 )
max.pic8 <- PIC <- 1 - ( (2/8)^2 + (3/8)^2 + (3/8)^2 )
max.pic9 <- PIC <- 1 - ( (3/9)^2 + (3/9)^2 + (3/9)^2 )
max.pic10 <- PIC <- 1 - ( (3/10)^2 + (3/10)^2 + (4/10)^2 )
max.pic11 <- PIC <- 1 - ( (3/11)^2 + (4/11)^2 + (4/11)^2 )
max.pic12 <- PIC <- 1 - ( (4/12)^2 + (4/12)^2 + (4/12)^2 )
max.pic13 <- PIC <- 1 - ( (4/13)^2 + (4/13)^2 + (5/13)^2 )

# Combine into list
max.pics <- list(max.pic1,max.pic2,max.pic3,max.pic4,max.pic5,max.pic6,max.pic7,max.pic8,max.pic9,max.pic10,max.pic11,max.pic12,max.pic13)

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

# Load actual ITH into dataframe. Repeat each value for the number of sample there are
actualith.test <- data.frame(lapply(data.frame(testPIC.frac$PIC.frac), rep, sample.count),  row.names = rownames(test.input))

# Use step.model to predict
predictedith.test <- data.frame(predict(step.model, test.input))

# Bind together predicted and actual
comparison <- as.data.frame(cbind(actualith.test, predictedith.test)) 
comparison <- na.omit(comparison)
colnames(comparison) <- c("actual","predicted")

# Add in patient names
comparison$patient <- sub('\\..*', '', rownames(comparison))

# Create labels to show adjR2 and pval
adjR2 <- round(summary(step.model)$adj.r.squared,3)
pval <- signif(lmp(step.model),3)
label1 <- paste('R[adj]^2 ==', adjR2)
label2 <- paste('p.value ==', pval)

# Plot predicted and actual
plot <- ggplot(data = comparison, aes(x = actual, y = predicted)) +
  geom_line(aes(group = patient), size=0.2, colour = "#003366") +
  geom_point(fill = "#003366", colour = "#003366", size = 5) +
  geom_line(aes(x = predicted, y = predicted), linetype = "dashed") +
  xlab("Patient CNA diversity") +
  ylab("Patient CNA diversity (predicted)") +
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = text,
        axis.text = text,
        legend.position = "none") +
  annotate('text', x = c(0), y = c(0.600,0.575), label = c(label1, label2), parse=TRUE, size = 8, hjust = 0, vjust = 1)
