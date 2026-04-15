###look at ANI, from fastANI

##Look at ANI
#read in dataframe from fastANI
PF_cali_ani_output <- read.delim("Pseudo_fluor/PFC_ANI.txt")
PF_cali_ani_output <- PF_cali_ani_output[-which(PF_cali_ani_output$Sp1 == '630A' | PF_cali_ani_output$Sp1 == '629A2A' | PF_cali_ani_output$Sp2 == '630A' | PF_cali_ani_output$Sp2 == '629A2A'),]
#load ggplot2
library(ggplot2)

#boxplot for all ANI
pond_ani<-
  
  ggplot(PF_cali_ani_output, aes(Lake1, ANI, fill=Comparison))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c("grey", 'orange'))+
  ylab("Average Nucleotide Identity")+
  xlab("")

mean(PF_cali_ani_output$ANI)
sd(PF_cali_ani_output$ANI)

t.test(PF_cali_ani_output$ANI ~ PF_cali_ani_output$Comparison)
#t = -5.8484, df = 347.07, p-value = 1.147e-08


ggarrange(dendrogram, pond_ani, labels=c("A", "B"), ncol=2)