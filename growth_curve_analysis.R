###growth curve analysis
#load in data
gcr<-read.delim("~/Documents/GitHub/Pseudo_fluor/cali_pf_growth_data.txt", header=T)

#remove isolate column name
gcr<-gcr[,-2]

#reshape data to be in correct format
library(reshape2)
gcr_cast<-dcast(gcr, time ~ Isolate2, value.var = 'OD', fun.aggregate = sum)

#analyze growth curves
library(growthcurver) 
gc_out <- SummarizeGrowthByPlate(gcr_cast, plot_file = '~/Documents/GitHub/Pseudo_fluor/pf_curves.pdf', plot_fit = T)

#t_gen is the fastest possible generation time
#the population size at the beginning of the growth curve is given by N 0. 
#the carrying capacity, is given by K
#The intrinsic growth rate of the population, r

##get means of the 4 metric above for each isolate
#reshape data for calculations
gc_out2<-gc_out[,c(1:4,6)]
gc_out2<-melt(gc_out2)
gc_out2$sample<-gsub("\\..*","", gc_out2$sample)

#calulate means
library(plyr)
gc_sum<-ddply(gc_out2, c('sample', 'variable'), mean=mean(value), value.var = 'value', fun.aggregate = sum)

#reshape data again
gc_sum2<-dcast(gc_sum, sample ~ variable, value.var = 'value', fun.aggregate = sum)

#write to file
write.table(gc_sum2, 'grow_curve_metrics.txt', sep='\t', quote=F, row.names = F)

#read in Bd inhibition data
bd_data<-read.delim('~/Documents/GitHub/Pseudo_fluor/bd_inhibition_data.txt', header=T)

#merge datasets
gc_bd<-merge(gc_sum2, bd_data, by.x='sample', by.y='Group')

#split by Bd strain
gc_bd_split<-split(gc_bd, gc_bd$Bd_strain)

#look at correlations for carrying capacity
cor.test(gc_bd_split$Jel197$k, gc_bd_split$Jel197$Percent_inhibition, method = 'spearman')
#S = 23151, p-value = 0.3973, rho=.1175462
plot(gc_bd_split$Jel197$k, gc_bd_split$Jel197$Percent_inhibition)

cor.test(gc_bd_split$JEL423$k, gc_bd_split$JEL423$Percent_inhibition, method = 'spearman')
#S = 17428, p-value = 0.01308, rho =0.3356826
plot(gc_bd_split$JEL423$k, gc_bd_split$JEL423$Percent_inhibition)

#look at correlations for growth rate
cor.test(gc_bd_split$Jel197$r, gc_bd_split$Jel197$Percent_inhibition, method = 'spearman')
#S = 24083, p-value = 0.5554, rho = 0.08203028
cor.test(gc_bd_split$JEL423$r, gc_bd_split$JEL423$Percent_inhibition, method = 'spearman')
#S = 20392, p-value = 0.1055, rho = 0.2227191

#look at correlations for generation time
cor.test(gc_bd_split$Jel197$t_gen, gc_bd_split$Jel197$Percent_inhibition, method = 'spearman')
#S = 28237, p-value = 0.5834, rho = -0.07630191 
cor.test(gc_bd_split$JEL423$t_gen, gc_bd_split$JEL423$Percent_inhibition, method = 'spearman')
#S = 16851, p-value = 0.007922, rho=0.3576


#make plot for SI
library(ggplot2)
library(ggpubr)

rate<-ggplot(gc_bd, aes(sample, r))+
  ylab("Intrinsic growth rate")+
  xlab("")+
  geom_bar(stat='identity')+
  theme_bw()+
  coord_flip()

gen_time<- ggplot(gc_bd, aes(sample, t_gen))+
  ylab("Generation time")+
  xlab("")+
  geom_bar(stat='identity')+
  theme_bw()+
  coord_flip()

OD<- ggplot(gc_bd, aes(sample, k))+
  ylab("Carrying capacity")+
  xlab("")+
  geom_bar(stat='identity')+
  theme_bw()+
  coord_flip()

ggarrange(rate, gen_time, OD, labels=c("A", "B", "C"))
