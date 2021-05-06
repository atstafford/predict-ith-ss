# Predicting diversity in adenomas and carcinomas from PGA

# Pull together PGA data for adenoma and carcinoma samples, and melt into long form
pgaPredict.list <- list(adenoma$pga, carcinoma$pga)
names(pgaPredict.list) <- c('adenoma','carcinoma')
pgaPredict.list <- lapply(pgaPredict.list, function(x) {
  x$sample <- rownames(x)
  x <- gather(x,CNA,proportion,-sample, -ith)
  x
})

# Construct plot for adenoma and carcinoma
plot.list <- list()
for ( i in 1:length(pgaPredict.list) ) {
  
  if (i == 1) {tit <- 'Adenoma'} #adenoma
  if (i == 2) {tit <- 'Carcinoma'} #carcinoma
  
  my.formula <- y ~ x
  
  plot.list[[i]] <- ggplot(data=pgaPredict.list[[i]], aes(y=ith, x=proportion, colour=CNA)) + 
    geom_smooth(method = "lm", se=FALSE, formula = my.formula) +
    geom_point() +
    ggpmisc::stat_poly_eq(formula = my.formula, 
                          aes(label = paste(stat(adj.rr.label), "*\", \"*", stat(p.value.label), "*\"\"", sep = "")),
                          parse = TRUE, label.x = 'right', label.y = 'top', size=7) +
    ggtitle(tit) +
    scale_fill_manual(values = c("#009999",'#CC0033','#330099'), labels = c('Aneu.','Gain','Loss')) + 
    scale_colour_manual(values = c("#009999",'#CC0033','#330099'), labels = c('Aneu.','Gain','Loss')) + 
    scale_y_continuous(limits = c(0,0.6)) +
    scale_x_continuous(limits = c(0,0.6), breaks = c(seq(0,0.6,0.1)), labels = c(seq(0,0.6,0.1))) +
    theme(plot.margin=margin(t=0,r=0.7,b=0.5,l=0.7,"cm"),
          panel.background = element_blank(),
          plot.title = text.bold,
          axis.line = element_line(colour = "black", size = 0.5),
          axis.title = element_blank(),
          axis.text = text,
          axis.ticks.length=unit(0.2, "cm"),
          legend.position = "none")  
}

# Prepare a common legend 
pga.legend <- cowplot::get_legend(plot.list[[i]] + 
                                    guides(fill=guide_legend(label.position = 'bottom', 
                                                             title.hjust = 0.5)) +
                                    theme(plot.margin = unit(c(0,0,-100,0), "cm"),
                                          legend.position = "top", legend.direction="horizontal",
                                          legend.title = text.bold,
                                          legend.margin = margin(grid::unit(c(0,0,-100,0),"cm")),
                                          legend.text = text,
                                          legend.key.height = grid::unit(0.8,"cm"),legend.key.width = grid::unit(2,"cm")) ) 

# Assemble plots
pgaPredict.plot <- cowplot::plot_grid(plot.list[[1]], plot.list[[2]], nrow = 2, align = "v")

pgaPredict.plot <- annotate_figure(pgaPredict.plot, 
                                   left = text_grob('Patient CNA diversity', rot = 90, size = 24),
                                   bottom = text_grob('Proportion of genome altered', size = 24))

pgaPredict.plot <- cowplot::plot_grid(as_ggplot(pga.legend), pgaPredict.plot, ncol = 1, align = "v", rel_heights = c(0.2,1))
