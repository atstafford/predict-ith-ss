# Heatmaps showing bin-bin pairwise correlation, with histograms of clonal/subclonal CNAs

# Prepare two histograms to display frequency of clonal and subclonal gains and losses
# Prepare long form of GaLo.Clo data for histograms
data <- gather(carcinoma$GaLo.clo,variable,freq,-bin,-chr1,-chr2)

# Prepare x axis histogram showing frequency of clonal gains and losses
text <- element_text(size=24, colour='black')
text.bold <- element_text(size=24, colour='black',face='bold')

bp.x <- ggplot(data[ which(data$variable == 'clonal.gain' | data$variable == 'clonal.loss'), ], aes(fill=variable, x=bin, y=freq)) + 
  geom_bar(position="stack", stat="identity", width = 1) +
  ggtitle('Clonal CNAs') +
  scale_fill_manual(values = c('#CC0033','#330099'),
                    labels = c('Gain','Loss')) + 
  
  facet_grid(~ chr1, space="free", scales="free") +
  
  scale_x_continuous(expand = c(0,0), name=NULL, breaks=NULL, labels = NULL) +
  scale_y_continuous(expand = c(0,0), name='Clonal \nfrequency', limits = c(0,0.8), breaks = c(0,0.4,0.8), labels = c(0,0.4,0.8)) +
  
  theme(plot.background=element_blank(),
        plot.margin=margin(t=0.4,r=0,b=1,l=0,"cm"),
        panel.border = element_rect(colour = 'black', size = 0.5, fill=NA),
        
        panel.spacing = unit(0.2, "lines"),
        panel.background = element_blank(),
        strip.text = element_blank(),
        
        legend.position = "none",
        
        plot.title = element_blank(),
        axis.title.y = text,
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = text,
        axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks.length=unit(0.2, "cm"),)

# Prepare x axis histogram showing frequency of subclonal gains and losses
bp.y <- ggplot(data[which(data$variable == 'subclonal.gain' | data$variable == 'subclonal.loss') , ], aes(fill=variable, x=bin, y=freq)) + 
  geom_bar(position="stack", stat="identity", width = 1) +
  scale_fill_manual(values = c('#CC0033','#330099'),
                    labels = c('gain','loss')) + 
  ggtitle('Subclonal CNAs') +
  
  facet_grid(chr2 ~., space="free", scales="free") +
  coord_flip() +
  
  scale_x_continuous(expand = c(0,0), name=NULL, breaks=NULL, labels = NULL) +
  scale_y_continuous(expand = c(0,0), name='Subclonal \nfrequency', limits = c(0,0.8), breaks = c(0,0.4,0.8), labels = c(0,0.4,0.8)) +
  
  theme(plot.background=element_blank(),
        plot.margin=margin(t=0,r=0.6,b=0,l=1,"cm"),
        panel.border = element_rect(colour = 'black', size = 0.5, fill=NA),
        
        panel.spacing = unit(0.2, "lines"),
        panel.background = element_blank(),
        strip.text = element_blank(),
        
        legend.position = "none",
        
        plot.title = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = text,
        axis.text.x = text,
        axis.text.y = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks.length=unit(0.2, "cm"),)


# Prepare data for heatmaps (x4) showing genome-wide correlations between:
# diploid/aneuploid, diploid/gain, diploid/loss, loss/gain, and loss/diploid/gain

# Having two multiregion samples per patient will drive up the correlation. 
# Therefore, we will consider the bin-bin correlations using the average copy number per patient
cor.inputs <- list()

# For each matrix 
for ( d in 1:length(car.matrices) ) {
  i <- 1
  l <- 1
  list <- list()
  
  # For a given bin
  for ( k in 1:nrow(car.matrices[[d]]) ) { 
    i <- 1
    
    # For a each patient, average the copy number between samples
    while ( i < ncol(car.matrices[[d]]) ) { 
      list[[l]] <- ((car.matrices[[d]][k,i]) + (car.matrices[[d]][k,i+1])) / 2 
      i <- i + 2
      l <- l + 1
    }
  }
  cor.inputs[[d]] <- data.frame(t(matrix(unlist(list), ncol=carcinoma$noPatients)), check.names = FALSE)
}

# Add list titles
names(cor.inputs) = c('diploid.aneu', 'diploid.gain', 'diploid.loss', 'loss.gain', 'loss.dip.gain')

# Convert to numeric matrices
cor.inputs <- lapply(cor.inputs, function(x) data.matrix(x,rownames.force = NA))


# Generate bin-bin corrleation matrices

# Calculate pairwise correlations
cormat.list <- lapply(cor.inputs, function(x) cor(x, method = 'pearson', use = 'pairwise.complete.obs'))

# Melt into long form for geom_tile
meltedcormat.list <- cormat.list
meltedcormat.list <- lapply(meltedcormat.list, function(x) melt(x, na.rm = FALSE))
meltedcormat.list <- lapply(meltedcormat.list, "colnames<-", c('first_bin','second_bin','correlation'))

# add 2x chr data for faceting
meltedcormat.list <- lapply(meltedcormat.list, function(x) {
  x <- cbind(x, chr1 = as.numeric(carcinoma$start.stop[match(x$first_bin, carcinoma$start.stop$bin), 2])) 
  x <- cbind(x, chr2 = as.factor(carcinoma$start.stop[match(x$second_bin, carcinoma$start.stop$bin), 2]))
  x$chr2 <- factor(x$chr2, levels = rev(c(seq(1,22,1))) ) # to order xaxis facets
  x
} )


# Prepare heatmaps 
text <- element_text(size=24, colour='black')
text.bold <- element_text(size=24, colour='black',face='bold')
hm.list <- list()

for ( i in 1:length(meltedcormat.list) ) {
  hm.list[[i]] <- ggplot(meltedcormat.list[[i]], aes(first_bin, second_bin, fill = correlation)) +
    geom_tile() +
    scale_fill_distiller(palette ="Spectral", na.value="grey85", limits = c(-1,1), breaks=c(-1,0,1), guide=guide_colorbar()) +
    facet_grid(chr2 ~ chr1, scales='free', space = 'free') +
    
    scale_x_continuous(name='Chromosome', breaks = c(round(adenoma$chr.mid, digits = 0)),  
                       labels = c(1:18,'\n19','20','\n21','22')) + 
    scale_y_continuous(name='Chromosome', breaks = c(round(adenoma$chr.mid, digits = 0)),
                       labels = c(1:18,'19    ','20','21    ','22')) +
    
    theme(plot.background=element_blank(),
          plot.margin=margin(t=0,r=0,b=0,l=0,"cm"),
          
          panel.spacing = unit(0.2, "lines"),
          panel.background = element_blank(),
          panel.border=element_blank(),
          strip.text = element_blank(),
          
          legend.position = "none",
          
          axis.text = text,
          axis.title = text,
          axis.ticks = element_line(size=0.4),
          axis.line = element_line(colour = "black", size = 0.5),
          axis.ticks.length=unit(0.2, "cm"),)
}

# Prepare legend for heatmaps and histograms 
bp.legend <- cowplot::get_legend(bp.x + guides(fill=guide_legend(label.position = 'bottom', label.vjust = 1,
                                                                 title.position = 'left', title.vjust = 0.9)) +
                                   theme(plot.margin = unit(c(-100,-100,-100,-100), "cm"),
                                         legend.position = "bottom",legend.direction="horizontal",
                                         legend.title = element_blank(),
                                         legend.margin = margin(grid::unit(c(t=-100,r=-100,b=-100,l=-100),"cm")),
                                         legend.text = text,
                                         legend.key.height = grid::unit(0.8,"cm"),
                                         legend.key.width = grid::unit(1.4,"cm")))

hm.legend <- cowplot::get_legend(hm.list[[1]] + 
                                   guides(fill = guide_colourbar(title = 'Correlation', label.position = "bottom", label.vjust = 2,
                                                                 title.position = 'top', title.hjust = 0.5, title.vjust = -1)) +
                                   theme(plot.margin = unit(c(-120,-100,-100,-100), "cm"),
                                         legend.position = "bottom",legend.direction="horizontal",
                                         legend.title = text.bold,
                                         legend.margin = margin(grid::unit(c(-120,-100,-100,-100),"cm")),
                                         legend.text = text,
                                         legend.key.height = grid::unit(0.8,"cm"),
                                         legend.key.width = grid::unit(1,"cm")) )

legend <- plot_grid(as_ggplot(bp.legend), as_ggplot(hm.legend), ncol=1, align = 'v', axis = 'lr')

# Assemble plots
# Each matrix will generate a different plot, stored in plot.left with the top histogram
plot.left <- list()
for ( i in 1:length(hm.list) ) {
  plot.left[[i]] <- plot_grid(bp.x, hm.list[[i]], nrow = 2, align = 'v', axis = 'l', rel_heights = c(2,10))
}

# The right of the plot will be constant, holding the left histogram and legend
plot.right <- plot_grid(legend, bp.y, nrow = 2, align = 'v', axis = 'lr', rel_heights = c(2,10))

# Bring left and right plots together
plot.list.final <- list()
for ( i in 1:length(plot.list) ) {
  plot.list[[i]] <- plot_grid(plot.left[[i]], plot.right, ncol = 2, rel_widths = c(10,2), align = 'h', axis = 'b')
}                        