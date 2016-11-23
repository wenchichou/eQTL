x=55;m=444;n=1248322;k=18768
(1-phyper(x, m, n, k, lower.tail = T))
1-0.9999999999999999

library(ggplot2)
require(cowplot)
library(gtable)
library(grid)
library(GGally)


d<-read.csv("~/Desktop/selectedEnrichment.txt",sep="\t")
head(d)
colnames(d)[13]="neglogP"
d$GWASpvalueRangeLarge=-log10(d$GWASpvalueRangeLarge)
str(d)
table(d$GWASpvalueRangeSmall)
ggplot(d, aes(x=GWASpvalueRangeLarge, y=neglogP, color=GWAS)) +  geom_line()+ geom_point(aes(size=x))+scale_x_continuous(breaks=rev(c(4,6,8,10,20)), trans="reverse")+facet_wrap(~GWAS)
  
ggplot(d, aes(x=GWASpvalueRangeLarge, y=neglogP, color=GWAS)) +  geom_line()+ geom_point()+scale_x_continuous(breaks=rev(c(4,6,8,10,20)))+facet_wrap(~GWAS)


  
  
  geom_path(data=metadata, aes(group = ID), cex=0.3) + geom_point(aes(colour=factor(linkedToUTI), shape=factor(takeAntibiotics)), group=1, cex=4) + scale_x_date(date_labels = "%b %y", date_breaks = "1 month")+ xlab("") + ylab("Subject ID") + facet_wrap(~ studyGroup, ncol = 1) + scale_y_reverse(breaks = pretty(metadata$ID, n = 31)) + theme_light() + scale_size(name="relative abundance %\nof E coli") + scale_colour_hue(name="UTI or UTI follow-ups", labels=c("UTI", "healthy")) + scale_shape_discrete(name="took antibiotics?", labels=c("No", "Yes"))



