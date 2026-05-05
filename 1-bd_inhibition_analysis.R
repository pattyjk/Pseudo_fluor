##################analayze Bd inhibition data
#load in Bd-inhibition data
bd_data<-read.delim("~/Documents/Github/Pseudo_fluor/bd_inhib_data_corrcted_p.txt")

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
