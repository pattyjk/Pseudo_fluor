###Calculate pangenome stats with micropan
#load packages
library(tidyverse)
library(micropan)

#read in data
pf.genes<-read.delim('~/Documents/GitHub/Pseudo_fluor/panaroo_results/gene_presence_absence2.Rtab', header=T, row.names=1)
panmat.blast <-t(pf.genes)

#calculate pangenome size
heaps.est <- heaps(panmat.blast, n.perm = 500)
print(heaps.est)
#Intercept       alpha 
#1105.797493    1.202049
#alpha is above 1, so pangenome is closed

print(chao(panmat.blast))
#9878

#rarefaction of pangenome, takes a while so set and do other stuff
rar.tbl <- rarefaction(panmat.blast, n.perm = 10000)

#aggregate data
rar.tbl2<-rar.tbl %>% 
  gather(key = "Permutation", value = "Clusters", -Genome)

#summarize and replot
library(plyr)
rar.sum<-ddply(rar.tbl2, c("Genome"), summarize, mean=mean(Clusters), sd=sd(Clusters))

ggplot(rar.sum, aes(Genome, mean))+
  geom_point()+
  ylab("Number of genes observed")+
  xlab("Number of genomes sampled")+
geom_smooth()+
  theme_bw()

#calculate average pangenome fluidity
pg_fluid<-fluidity(panmat.blast, n.sim = 10000)
#Mean        Std 
#0.07576517 0.03672605 
#For example, a genomic fluidity of 0.1 represents that a pair of genomes have on average 10% unique genes 
#and share 90% of their genes.
