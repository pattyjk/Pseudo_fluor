################Antismash analysis
###AntiSmash
anti_out<-read.delim("~/Github/Pseudo_fluor/AntiSmash_output/antismash_out.txt")
library(dplyr)
library(plyr)

#summarize
anti_sum<-ddply(anti_out, c("Genome", "Pond", "Type"), summarize, no_genes=length(Type))
ggplot(anti_sum, aes(Genome, no_genes, fill=Type))+
  geom_bar(stat='identity')+
  theme_bw()+
  ylab("Number of Genes")+
  xlab("")+
  #theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust=0.5))+
  coord_flip()+
  scale_fill_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
  facet_wrap(~Pond, scales='free_y')


#append Bd inhibition data to antismash data
bd_data<-read.delim("~/Github/Pseudo_fluor/bd_inhibition_data.txt")

#remove controls from data (blanks have been subtracted out already)
bd_data<-bd_data[-which(bd_data$Group =='Heat killed' | bd_data$Group =='Tryptone' | bd_data$Group =='Bd_Prova' | bd_data$Group =='Bd_Tryptone'),]

#calculate total number of BCG per genome
anti_sum2<-ddply(anti_sum, c("Genome"), summarize, n_clust=length(no_genes))

#calculate summary stats for BD inhibition data
bd_sum<-ddply(bd_data, c("Bd_strain", "Group"), summarize, mean=mean(Percent_inhibition), sd=sd(Percent_inhibition), n=length(Percent_inhibition), se=sd/n)
anti_sum_bd<-merge(anti_sum2, bd_sum, by.x='Genome', by.y='Group')

#plot clusters and Bd data
library(ggpubr)
ggplot(anti_sum_bd, aes(n_clust, mean, color=Bd_strain))+
  geom_point()+
  ylab("Percent Bd inhibition")+
  xlab("Number of BGC")+
  geom_smooth(method='lm')+
  stat_cor(method='spearman')