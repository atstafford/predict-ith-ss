# Download libraries
Packages <- c("tidyr", "ggplot2","ggrepel","ggpubr","cowplot","reshape2","ggpmisc","pvclust","MASS","stringr","dplyr","readxl","corrplot","survival","survminer")
lapply(Packages, library, character.only = TRUE)

# The original dataset from Cross et al has CN data across 2694 bins for 19 adenomas and 81 carcinomas
# Two sample from each lesion we taken
# Columns are: chr | start | stop | sample1 | sample2...
# Rows will contain the raw copy number assignment per bin (1=loss, 2=diploid, 3=gain). 
# From columns 4 onwards, sample identifiers must be constructed as: patientName.sampleNum

# Load data
rawdata <- read.table("~/Documents//PhD Barts/CNA AdVsCar/Data/sWGS Cross rename with type.txt", header=T, skip=2, sep="\t")

# Split into 2 separate raw data files (ad.raw and car.raw) 
ad.raw<- rawdata[,c(1:3,which(sub("^([[:alpha:]]*).*", "\\1", colnames(rawdata))=="A"))]
car.raw <- rawdata[,c(1:3,which(sub("^([[:alpha:]]*).*", "\\1", colnames(rawdata))=="C"))]

# Use dataSetup function to create necessary dataframes for analysis
adenoma <- dataSetup(rawdata = ad.raw)
carcinoma <- dataSetup(rawdata = car.raw)

# Use genMatrices to generate four matrices from the raw data (where 1=loss, 2=noCNA, 3=gain).
ad.matrices <- genMatrices(ad.raw)
car.matrices <- genMatrices(car.raw)




