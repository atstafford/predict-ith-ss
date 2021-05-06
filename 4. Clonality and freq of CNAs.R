# This analysis assess the clonality and frequency of alterations and adenomas and carcinomas

# Bind together the adenoma and carcinoma dataframes holding the fraction of alterations at each bin
CloFreq <- rbind( cbind(data.frame(type='ad'),adenoma$GaLo.clo[,c(1,6,7)]), 
                  cbind(data.frame(type='car'),carcinoma$GaLo.clo[,c(1,6,7)]))

# Convert freq into fraction of patients
CloFreq[which(CloFreq$type=='ad'),c(3,4)] <- CloFreq[which(CloFreq$type=='ad'),c(3,4)] / adenoma$noPatients
CloFreq[which(CloFreq$type=='car'),c(3,4)] <- CloFreq[which(CloFreq$type=='car'),c(3,4)] / carcinoma$noPatients

# Make losses negative to create a mirror plot
CloFreq$loss <- CloFreq$loss*-1 

# Gather into long form
CloFreq <- gather(CloFreq,CNA,freq,-bin,-type) 

# Pull the % subclonal for a gain or loss at each bin from GaLo.status for adenoma and carcinoma
CloFreq$pcSubclonal <- NA

for ( i in 1: nrow(CloFreq) ) {
  type <- CloFreq$type[i]
  bin <- CloFreq$bin[i]
  CNA <- CloFreq$CNA[i]
  
  # GaLo.status dataframe holds the clonality for gains in col 2
  if ( type == 'ad' & CNA == 'gain' ) {
    CloFreq$pcSubclonal[i] <- adenoma$GaLo.clo$pcSubclonal.gain[bin] *100
  }
  if ( type == 'car' & CNA == 'gain' ) {
    CloFreq$pcSubclonal[i] <- carcinoma$GaLo.clo$pcSubclonal.gain[bin] *100
  }
  
  # GaLo.status dataframe holds the clonality for losses in col 3
  if ( type == 'ad' & CNA == 'loss' ) {
    CloFreq$pcSubclonal[i] <- adenoma$GaLo.clo$pcSubclonal.loss[bin] *100
  }
  if ( type == 'car' & CNA == 'loss' ) {
    CloFreq$pcSubclonal[i] <- carcinoma$GaLo.clo$pcSubclonal.loss[bin] *100
  }
}

# Add labels identifying cancer-associated genes. The x coordinate is the gene location
annotation_ad <- data.frame( #common CRC genes
  x = c(960,1253,1311,1479,1808,2035,2360,2387,2480),
  label = c('APC','EGFR','MET','MYC','CCND1','RB1','TP53','ERBB2','SMAD4'),
  y = NA)
annotation_car <- data.frame( #common CRC genes
  x = c(960,1253,1311,1479,1808,2035,2360,2387,2480),
  label = c('APC','EGFR','MET','MYC','CCND1','RB1','TP53','ERBB2','SMAD4'),
  y = NA)

# Position label along y axis
for ( i in 1:nrow(annotation_ad) ) {
  bin <- annotation_ad$x[i]
  data <- CloFreq[which(CloFreq$type=='ad' & CloFreq$bin==bin),]
  annotation_ad$y[i] <- data$freq[which.max(abs(data$freq))]
}
for ( i in 1:nrow(annotation_car) ) {
  bin <- annotation_car$x[i]
  data <- CloFreq[which(CloFreq$type=='car' & CloFreq$bin==bin),]
  annotation_car$y[i] <- data$freq[which.max(abs(data$freq))]
}

# Create plots
fig.option <- c('ad','car')
figs <- list()
text <- element_text(size=24, colour='black')
text.bold <- element_text(size=24, colour='black',face='bold')

for ( i in 1:length(fig.option)) {
  if (fig.option[i] == 'ad') { gtit <- 'Adenoma'; annotation <- annotation_ad; xtext <- text }
  if (fig.option[i] == 'car') { gtit <- 'Carcinoma'; annotation <- annotation_car; xtext <- text }
  
  figs[[i]] <- ggplot() + 
    geom_bar(data=CloFreq[which(CloFreq$type==fig.option[i]),], aes(fill=pcSubclonal, y=freq, x=bin),
             position="stack", stat="identity", width = 1) +
    scale_fill_distiller(palette = "PRGn") +
    
    ylab("fraction of patients with aberration") +
    ggtitle(gtit) +
    scale_x_continuous(expand = c(0,0), name="chromosome", breaks=carcinoma$chr.mid, labels = c(1:18,'\n19','20','\n21','22')) +
    scale_y_continuous(limits = c(-0.8,0.8), breaks = c(seq(-0.8,0.8,0.2)), 
                       labels = c(0.8,0.6,0.4,0.2,0,0.2,0.4,0.6,0.8)) +
    
    geom_hline(yintercept = 0, size=0.3, color="black") +
    geom_rect(data = carcinoma$chr.end[which(carcinoma$chr.end$col=='W'),], 
              aes(NULL,NULL,xmin=start, xmax=end),
              fill = alpha("#CC9966", 0.1),
              ymin = -1,
              ymax = 1) +
    
    theme(plot.margin = unit(c(t=0,r=0.2,b=-0.2,l=0.7), "cm"),
          panel.background = element_blank(),
          plot.title = text.bold,
          axis.title = element_blank(),
          axis.text.y = text,
          axis.text.x = xtext,
          axis.line = element_line(colour = "black", size = 0.5),
          axis.ticks.length=unit(0.2, "cm"),
          legend.position = "none") +
    
    geom_point(data = annotation, aes(x = x, y = y),
               shape = 18, size = 0) +
    geom_text_repel(data = subset(annotation, y >= 0),
                    aes(x = x, y = y, label = label),
                    size = 7,
                    nudge_y = 0.2,direction = 'x',
                    arrow = arrow(length = unit(0.015, "npc"))) +
    geom_text_repel(data = subset(annotation, y < 0),
                    aes(x = x, y = y, label = label),
                    size = 7,
                    nudge_y = -0.2,direction = 'x',
                    arrow = arrow(length = unit(0.015, "npc"))) +
    annotate('text', y=c(0.8,-0.8), x=c(50,50), 
             label=c('italic(Gains)','italic(Losses)'), size=7, parse=TRUE, hjust = 0)
}

# Create legend
galo.legend <- cowplot::get_legend(figs[[1]] + 
                                     guides(fill = guide_colorbar(title="% subclonal", label.position = "bottom",
                                                                  title.position = "left", title.vjust = 0.9)) +
                                     theme(plot.margin = unit(c(0,0,-100,0), "cm"),
                                           legend.position = "bottom",legend.direction="horizontal",
                                           legend.title = text.bold,
                                           legend.margin = margin(grid::unit(c(0,0,-100,0),"cm")),
                                           legend.text = text,
                                           legend.key.height = grid::unit(0.8,"cm"),
                                           legend.key.width = grid::unit(1.4,"cm")) )

galo.legend <- as_ggplot(galo.legend)

# Plot together
plot <- cowplot::plot_grid(galo.legend,
                           annotate_figure(plot_grid(figs[[1]], figs[[2]], ncol = 1, align = 'v', rel_heights = c(1,1)),
                                           left = text_grob('Fraction of samples with CNA', size = 24, rot = 90),
                                           bottom = text_grob('Chromosome', size = 24)),
                           ncol = 1, align = 'v', rel_heights = c(0.2, 1))

# Is clonality correlated with frequency in adenomas
x <- CloFreq
x$freq[which(x$freq<0)] <- x$freq[which(x$freq<0)]*-1
x <- x[which(x$type=='ad'),]
cor.test(x$freq, x$pcSubclonal, type = 'pearson')

# Is clonality correlated with frequency in adenomas
x <- CloFreq
x$freq[which(x$freq<0)] <- x$freq[which(x$freq<0)]*-1
x <- x[which(x$type=='car'),]
cor.test(x$freq, x$pcSubclonal, type = 'pearson')
