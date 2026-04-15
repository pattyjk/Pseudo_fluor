##################analayze Bd inhibition data
#load in Bd-inhibition data
bd_data<-read.delim("~/Documents/Github/Pseudo_fluor/bd_inhibition_data.txt")

#remove controls from data (blanks have been subtracted out already)
bd_data<-bd_data[-which(bd_data$Group =='Heat killed' | bd_data$Group =='Tryptone' | bd_data$Group =='Bd_Prova' | bd_data$Group =='Bd_Tryptone'),]

#correct p-values
bd_data$p.adj<-p.adjust(bd_data$P_value, method = 'hochberg')

#write table with corrected p-values
write.table(bd_data, '~/Documents/Github/Pseudo_fluor/bd_inhib_data_corrcted_p.txt', quote=F, row.names = F, sep='\t')
bd_data<-read.delim('~/Documents/Github/Pseudo_fluor/bd_inhib_data_corrcted_p.txt', header=T)

#plot data
library(ggplot2)
ggplot(bd_data, aes(Group, Percent_inhibition, fill=Bd_strain))+
  geom_boxplot()+
  coord_flip()+
  scale_fill_manual(values=c("#4C72B0", "#55A868"))+
  xlab("")+
  facet_wrap(~Pond, scales='free_y')+
  theme_bw()+
  ylab("Percent inhibition of Bd")

#plot histogram
ggplot(bd_data, aes(Percent_inhibition, fill=Pond))+
  geom_histogram()+
  scale_fill_manual(values=c("#4C72B0", "#55A868"))+
  facet_wrap(~Bd_strain)+
  ylab("Frequency")+
  xlab("Percent Bd inhibition")
