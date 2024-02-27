###############Pangenome shells

#read in break down of gene shells for all genomes
library(ggplot2)
library(ggpubr)
summary_statistics <- read.delim("~/GitHub/Pseudo_fluor/Pangenome_results/PFC_panaroo_results//summary_statistics.txt", header=FALSE)
all<-ggplot(summary_statistics[-5,], aes(x="", y=V3, fill=V2))+
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  ylab("")+
  theme_bw()+
  xlab("")

#read in shell break down of Sixty lake gene shells
sixty_sum_stats<-read.delim("~/GitHub/Pseudo_fluor/Pangenome_results/PFC_panaroo_results_sixty/summary_statistics.txt", header=F)
six<-  ggplot(sixty_sum_stats[-5,], aes(x="", y=V3, fill=V2))+
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  ylab("")+
  theme_bw()+
  xlab("")

#read in shell break down of Conness Pond gene shells
coness_sum_stats<-read.delim("~/GitHub/Pseudo_fluor/Pangenome_results/PFC_panaroo_results_coness/summary_statistics.txt", header=F)
con<-ggplot(coness_sum_stats[-5,], aes(x="", y=V3, fill=V2))+
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  theme_bw()+
   ylab("")+
  xlab("")

#plot all together
ggarrange(all, con, six, ncol=3, nrow=1, labels=c("A- all isolates", "B- Conness Pond", "C- Sixty Lake"), common.legend = T)

