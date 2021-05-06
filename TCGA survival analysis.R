# Fraction of TCGA patients with high/low predicted CNA diversity by stage

# Load TCGA data
COAD <- read.table("~/Documents//PhD Barts/CNA AdVsCar/Data/TCGA/TCGA_mastercalls.abs_segtabs.fixed.txt", header=T, skip=0, sep="\t")

# Load clinical data for CRC (COAD) and remove duplicates
COAD_clinical <- read_excel("~/Documents/PhD Barts/CNA AdVsCar/Data/TCGA/Clinical COAD.xlsx")
COAD_clinical <- COAD_clinical[!duplicated(COAD_clinical$case_submitter_id), ]

# Keep only COAD
COAD <- COAD[which((str_sub(COAD$Sample,1,str_length(COAD$Sample)-3)) %in% COAD_clinical$case_submitter_id),c(1:4,9)]

# Remove patients missing CN data
COAD <- COAD[!is.na(COAD$Modal_Total_CN),]

# Rename columns
colnames(COAD) <- c('sample','chr','start','stop','cn')


# Separate patients into separate dataframes held in a list

# Extract sample IDs (patient)
names <- unique(COAD$sample)
COAD.list <- list()

# Pull cn data into element corresponding to sample ID
for ( i in 1:length(names) ) {
  COAD.list[[i]] <- COAD[which(COAD$sample == names[i]),]
}

# Assess ploidy and recentre if average CN is 2.8 or more
COAD.list <- ploidyRecentre(COAD.list)

# Bin CN data to match training dataset
COAD.binned <- alignBins(bins = carcinoma$start.stop, cn.list = COAD.list)

# If absolute copy number >=3, its a gain, if <=2 its a loss
COAD.binned <- lapply(COAD.binned, function(x) {
  x[,-c(1:4)] <- ifelse(x[,-c(1:4)] >= 3, 3, ifelse(x < 2, 1, 2))
  x
})

# Rename columns with sampleID and combine all patients into one dataframe
COAD.binned <- COAD.binned %>% purrr::reduce(full_join, by = c("chr","start","stop","bin"))
colnames(COAD.binned)[-c(1:4)] <- c(1:(ncol(COAD.binned)-4)) 


# Predict

# Transpose as input table requires cols=bins
TCGA.input <- data.frame(t(COAD.binned[,-c(1:4)]), check.names = FALSE)

# Convert to characters
TCGA.input <- ifelse(TCGA.input == 3, "gain", ifelse(TCGA.input == 1, "loss", "diploid"))

# Convert to dataframe
TCGA.input <- data.frame(TCGA.input, check.names = FALSE)

# Change to factor
TCGA.input <- data.frame(lapply(TCGA.input, as.factor), row.names = rownames(TCGA.input), check.names = FALSE)

# Predict using step.model
TCGA.predict <- data.frame(case_submitter_id = str_sub(names,1,str_length(names)-3), predictedITH = predict(step.model, TCGA.input))

# Remove NAs
TCGA.predict <- TCGA.predict[!is.na(TCGA.predict$predictedITH),]

# Merge predicted ITH with clinical data
TCGA.predict <- merge(TCGA.predict, COAD_clinical, by = 'case_submitter_id')


# How many high/low patients by stage

# Keep only those samples whose ITH could be predicted
IthByStage <- TCGA.predict[!is.na(TCGA.predict$predictedITH),]

# Remove samples missing stage data
IthByStage <- IthByStage[which(IthByStage$ajcc_pathologic_stage!="'--"),]

# How many samples lack alterations in any predictive bins 
length(x[which(x$predictedITH!=coef(step.model)["(Intercept)"]),])

# Divide into low and high diversiyt by median
IthByStage$divGroup <- 'low' 
IthByStage$divGroup[IthByStage$predictedITH>= median(x$predictedITH) ] <- "High"

# Change stage names to remove A/B
colnames(IthByStage)[28] <- 'stage'
for ( i in 1:nrow(IthByStage) ) {
  IthByStage$stage[i] <- ifelse(IthByStage$stage[i]=="Stage IA", "Stage I", 
                                ifelse(IthByStage$stage[i]=="Stage IIA", "Stage II",
                                       ifelse(IthByStage$stage[i]=="Stage IIB", "Stage II",
                                              ifelse(IthByStage$stage[i]=="Stage IIIA", "Stage III",
                                                     ifelse(IthByStage$stage[i]=="Stage IIIB", "Stage III",
                                                            ifelse(IthByStage$stage[i]=="Stage IIIC", "Stage III",
                                                                   ifelse(IthByStage$stage[i]=="Stage IVA", "Stage IV",
                                                                          ifelse(IthByStage$stage[i]=="Stage IVB", "Stage IV",IthByStage$stage[i]))))))))
}

# Perform chi squared test
chisq <- chisq.test(table(IthByStage$divGroup, IthByStage$stage))
corrplot(chisq$residuals, is.cor = FALSE)

# Pull diversity by stage into summary table
IthByStage_summary <- as.data.frame(table(IthByStage$stage, IthByStage$divGroup))

# Plot
text <- element_text(size=24, colour='black')
text.bold <- element_text(size=24, colour='black',face='bold')

divByStagePlot <- ggplot(IthByStage_summary, aes(fill=Var2, y=Freq, x=Var1)) + 
  geom_bar(position="fill", stat="identity", width=0.8) +
  scale_fill_manual(values = c('#00CC66','#FF9933'), labels = c('High','Low')) + 
  ylab('Fraction') +
  labs(fill = 'Predicted CNA diversity') +
  guides(fill = guide_legend(direction = 'horizontal', label.position = "left", label.hjust = 1)) +
  scale_x_discrete(labels=c('Stage I\n(n=61)','Stage II\n(n=133)','Stage III\n(n=99)','Stage IV\n(n=47)')) +
  theme(plot.background=element_blank(),
        plot.margin=margin(t=0.4,r=0,b=1,l=0,"cm"),
        panel.background = element_blank(),
        legend.position = "top",
        legend.text = text,
        legend.title = text.bold,
        plot.title = element_blank(),
        axis.title.y = text,
        axis.title.x = element_blank(),
        axis.text = text,
        axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks.length=unit(0.2, "cm"),)

# Survival of high/low diversity patients (KM) 

# Calculate serial time between death and censor

# Convert year of birth/death to numeric
IthByStage$year_of_death <- as.numeric(as.character(IthByStage$year_of_death))
IthByStage$year_of_birth <- as.numeric(as.character(IthByStage$year_of_birth))

# Convert age at diagnosis from days to years
IthByStage$age_at_diagnosis <- as.numeric(as.character(IthByStage$age_at_diagnosis))
IthByStage$age_at_diagnosis <- IthByStage$age_at_diagnosis/365

# Calculate the year at diagnosis
IthByStage$year_at_diagnosis <- IthByStage$year_of_birth + floor(IthByStage$age_at_diagnosis)

# Remove patients who are dead with no year of death
IthByStage <- IthByStage[-(which(IthByStage$vital_status=='Dead' & is.na(IthByStage$year_of_death)==TRUE)),]

# Calculate serial time (yr death - yr diagnosis)
IthByStage$serialTime <- IthByStage$year_of_death - IthByStage$year_at_diagnosis

# If no year of death given, assume survival until 2012
IthByStage$serialTime[is.na(IthByStage$serialTime)] <- 2012 - IthByStage$year_at_diagnosis

# Censor code to show 1=dead, 0=censor (no year of death)
IthByStage$censored <- 1
IthByStage$censored[is.na(IthByStage$year_of_death)] <- 0

# Fit survival data using the Kaplan-Meier method
stage = c("Stage I", "Stage II", "Stage III", "Stage IV")
plot.list <- list()
for (i in 1:length(stage)) {
  surv_object <- Surv(time = IthByStage$serialTime[which(IthByStage$stage==stage[i])], 
                      event = IthByStage$censored[which(IthByStage$stage==stage[i])])
  fit <- survfit(surv_object ~ divGroup, data = IthByStage[which(IthByStage$stage==stage[i]),])
  plot.list[[i]] <- ggsurvplot(fit, data = IthByStage[which(IthByStage$stage==stage[i]),], pval = TRUE, 
                               title=stage[i], legend="none",
                               font.main=24, font.x=NA, font.y=NA, font.tickslab=24, font.legend=24, pval.size=9, 
                               palette = c('#00CC66','#FF9933'), legend.title='Predicted CNA diversity', legend.labs=c('high','low'))
}

# Plot
plot <- plot_grid(plot.list[[1]]$plot, plot.list[[2]]$plot, plot.list[[3]]$plot, plot.list[[4]]$plot, nrow = 2)

plot <- plot_grid(get_legend(divByStagePlot), plot, nrow = 2, rel_heights = c(1,10))

plot <- annotate_figure(plot, 
                        left = text_grob('Survival probability', rot = 90, size = 24),
                        bottom = text_grob('Time (years)', size = 24))