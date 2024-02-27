###Calculating biomodality for Bd inhibition data

###############Look at AmphiBac database
#read in data & library
library(ggplot2)
library(tidyverse)
amphibac<-read.delim("Pseudo_fluor/amphibac_database.txt", header=T)

#remove NAs
amphibac2<-amphibac %>% drop_na(Proportional.Growth.Bd)

#plot all
ggplot(amphibac2, aes(as.numeric(Proportional.Growth.Bd)))+
  geom_histogram(bins = 50)+
  #coord_cartesian(xlim=c(0,2))+
  facet_wrap(~O, scales = 'free')

#subset to only taxa with high number of isolates (>50)
amphibac_subset<-amphibac2[which(amphibac2$O == 'o__Pseudomonadales'| amphibac2$O =='o__Enterobacteriales' |
                                   amphibac2$O =='o__Aeromonadales' |
                                   amphibac2$O =="o__Sphingobacteriales"  |
                                   amphibac2$O =="o__Sphingomonadales" |
                                   amphibac2$O =="o__Burkholderiales"  |
                                   amphibac2$O =="o__Flavobacteriales" |
                                   amphibac2$O =="o__Caulobacterales" |
                                   amphibac2$O =="o__Actinomycetales"| amphibac2$O =="o__Bacillales"),]

amphibac_subset$Proportional.Growth.Bd<-as.numeric(amphibac_subset$Proportional.Growth.Bd)                                
amphibac_subset<- amphibac_subset %>% drop_na(Proportional.Growth.Bd)
amphibac_subset<- amphibac_subset[-which(amphibac_subset$Proportional.Growth.Bd >2),]

#plot those of interest
ggplot(amphibac_subset, aes(as.numeric(Proportional.Growth.Bd), fill=O))+
  geom_histogram(bins = 50)+
  facet_wrap(~O, scales = 'free_y')+
  xlab("Proportional Bd Growth")+
  ylab("Count")

#boxplot to see range
ggplot(amphibac_subset, aes(O, as.numeric(Proportional.Growth.Bd), fill=O))+
  geom_boxplot()+
  facet_wrap(~O, scales = 'free')

#fit distribution
library(fitdistrplus)
library(stats4)
library(MASS)

pseudo<-amphibac_subset[which(amphibac_subset$O == 'o__Pseudomonadales'),]

x11()
plotdist(pseudo$Proportional.Growth.Bd, histo=T, demp=T)
descdist(pseudo$Proportional.Growth.Bd, boot = 1000)

#calculate diptest for biomodality
library(diptest)
#https://cran.r-project.org/web/packages/diptest/diptest.pdf
set.seed(515)
dip.test(pseudo$Proportional.Growth.Bd, simulate.p.value = T, B=10000)
#D = 0.021419, p-value = 0.0309

dip(pseudo$Proportional.Growth.Bd, full.result = T)
#n = 596.  Dip statistic, D_n = 0.02141937 = 25.53189/(2n)
#Modal interval [xL, xU] = [x[29], x[80]] = [0, 0]

######look at own data
bd_data<-read.delim('Pseudo_fluor/bd_inhib_data_corrcted_p.txt', header=T)

#separate to 197/423
bd_197<-bd_data[which(bd_data$Bd_strain=='Jel197'),]
bd_423<-bd_data[which(bd_data$Bd_strain=='JEL423'),]

#plot histogram
ggplot(bd_data, aes(Percent_inhibition, fill=Bd_strain))+
  geom_histogram()+
  ylab("Count")+
  xlab("Percent Bd inhibition")+
  scale_fill_manual(values=c("grey", 'black'))

dip.test(bd_197$Percent_inhibition, simulate.p.value = T, B=10000)
#D = 0.094202, p-value = 1e-04

dip.test(bd_423$Percent_inhibition, simulate.p.value = T, B=10000)
#D = 0.044946, p-value = 0.4206

dip.test(bd_data$Percent_inhibition, simulate.p.value = T, B=10000)
#D = 0.065823, p-value = 1e-04

#globally there is biomodality, but not for 423

##calculate dip to get quantiles, use the full dataset
dip(bd_data$Percent_inhibition, full.result = T)
#n = 132.  Dip statistic, D_n = 0.0658227 = 17.37719/(2n)
#Modal interval [xL, xU] = [x[21], x[68]] = [36.88545, 52.19846]

dip(bd_197$Percent_inhibition, full.result = T)
#n = 66.  Dip statistic, D_n = 0.09420202 = 12.43467/(2n)
#Modal interval [xL, xU] = [x[1], x[23]] = [41.82438, 50.92756]

#not much difference between estimates for 197/whole dataset

#calculate means for Bd inhibition data for 197 to bin into new groups
library(plyr)
bd_summary<-ddply(bd_197, c('Group', 'Pond'), summarize, mean=mean(Percent_inhibition), sd=sd(Percent_inhibition), n=length(Percent_inhibition), se=sd/n)

#lower group (mildly inhibitory) < 41%
#upper group (inhibitory) >50%