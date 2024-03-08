###growth curve analysis
#load in data
gcr<-read.delim("Pseudo_fluor/cali_pf_growth_data.txt", header=T)

#remove isolate column name
gcr<-gcr[,-2]

#reshape data to be in correct format
library(reshape2)
gcr_cast<-dcast(gcr, time ~ Isolate2, value.var = 'OD', fun.aggregate = sum)

#analyze growth curves
library(growthcurver) 
gc_out <- SummarizeGrowthByPlate(gcr_cast, plot_file = 'pf_curves.pdf', plot_fit = T)

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

#read in data with Bd inhibition data
growth_curve<-read.delim("Pseudo_fluor/grow_curve_metrics_with_inhib.txt", header=T)

#look at correlations for carrying capacity
cor.test(growth_curve$k, growth_curve$Inhib_197, method = 'spearman')
#t = -0.10882, df = 19, p-value = 0.9145, cor = -0.02495685
cor.test(growth_curve$k, growth_curve$Inhib_423, method = 'spearman')
#t = 0.66906, df = 19, p-value = 0.5115, cor = 0.1517166

#look at correlations for growth rate
cor.test(growth_curve$r, growth_curve$Inhib_197, method = 'spearman')
#t = 0.89151, df = 19, p-value = 0.3838, cor = 0.200378
cor.test(growth_curve$r, growth_curve$Inhib_423, method = 'spearman')
#t = 1.199, df = 19, p-value = 0.2453, cor = 0.2652281

#look at correlations for generation time
cor.test(growth_curve$t_gen, growth_curve$Inhib_197, method = 'spearman')
#t = 1.5268, df = 19, p-value = 0.1433, cor = -0.07795774
cor.test(growth_curve$t_gen, growth_curve$Inhib_423, method = 'spearman')
#t = 1.199, df = 19, p-value = 0.2453, cor = 0.330582

#make plot for SI
library(ggplot2)
library(ggpubr)

rate<-ggplot(growth_curve, aes(sample, r))+
  ylab("Intrinsic growth rate")+
  xlab("")+
  geom_bar(stat='identity')+
  theme_bw()+
  coord_flip()

gen_time<- ggplot(growth_curve, aes(sample, t_gen))+
  ylab("Generation time")+
  xlab("")+
  geom_bar(stat='identity')+
  theme_bw()+
  coord_flip()

OD<- ggplot(growth_curve, aes(sample, k))+
  ylab("Carrying capacity")+
  xlab("")+
  geom_bar(stat='identity')+
  theme_bw()+
  coord_flip()

ggarrange(rate, gen_time, OD, labels=c("A", "B", "C"))
