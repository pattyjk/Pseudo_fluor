###look at ANI, from fastANI

##Look at ANI
#read in dataframe from fastANI
PF_cali_ani_output <- read.delim("~/Documents/GitHub/Pseudo_fluor/pfcali_ani_output.txt", header=T)
#load ggplot2
library(ggplot2)

#make new column if they comparision is within a pond or between ponds
PF_cali_ani_output$Comparison<-ifelse(PF_cali_ani_output$Lake1 == PF_cali_ani_output$Lake2, "Within a lake", 'Between lakes')


#boxplot for all ANI


ggplot(PF_cali_ani_output, aes(Lake1, ANI, fill=Comparison))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values=c("#4C72B0", "#55A868"))+
  ylab("Average Nucleotide Identity")+
  xlab("")

mean(PF_cali_ani_output$ANI)
#[1] 99.19844
sd(PF_cali_ani_output$ANI)
#[1] 0.4609539

t.test(PF_cali_ani_output$ANI ~ PF_cali_ani_output$Comparison)
#t = -5.8484, df = 347.07, p-value = 1.147e-08

