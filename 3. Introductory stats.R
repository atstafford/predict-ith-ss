# To assess any significant difference in the PGA, average diversity, and clonality of alteration between ads and cars

# t-test to compare the average PGA between adenomas and carcinomas:
t.test(carcinoma$pga$prop.aneu, adenoma$pga$prop.aneu,
       alternative = 'two.sided', var.equal = FALSE)

# t-test and violin to compare the average PIC frac between adenomas and carcinomas:
t.test(carcinoma$ith$PIC.frac, adenoma$ith$PIC.frac,
       alternative = 'two.sided', var.equal = FALSE)

x <- rbind(carcinoma$ith,adenoma$ith)
x$type <- c(rep('car',81),rep('ad',19))
ggplot(x, aes(x=type, y=PIC.frac, colour=type)) +
  geom_violin(trim=FALSE, size=1) +
  geom_jitter(shape=16, size = 5, position = position_jitter(0.1)) +
  stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.05, colour='black', size = 1)+
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = text,
        axis.text = text,
        legend.position = "none") 

# t-test to compare the number of clonal CNAs between adenomas and carcinomas, as a proportion of PGA
t.test(carcinoma$patientClo$clonal/carcinoma$patientClo$CNA, adenoma$patientClo$clonal/adenoma$patientClo$CNA,
       alternative = 'two.sided',var.equal = FALSE)

# t-test to compare the number of subclonal CNAs between adenomas and carcinomas, as a proportion of PGA
t.test(carcinoma$patientClo$subclonal/carcinoma$patientClo$CNA, adenoma$patientClo$subclonal/adenoma$patientClo$CNA,
       alternative = 'two.sided',var.equal = FALSE)

# t-test to compare the number of clonal and subclonal CNAs in adenoma patients
t.test(adenoma$patientClo$subclonal/adenoma$patientClo$CNA,adenoma$patientClo$clonal/adenoma$patientClo$CNA,
       alternative = 'two.sided',var.equal = FALSE)

# t-test to compare the number of clonal and subclonal CNAs in carcinoma patients
t.test(carcinoma$patientClo$subclonal/carcinoma$patientClo$CNA,carcinoma$patientClo$clonal/carcinoma$patientClo$CNA,
       alternative = 'two.sided',var.equal = FALSE)