###Calculate pangenome with micropan 
#load packages
library(tidyverse)
library(micropan)

#setwd
setwd("C:/Users/patty/OneDrive/Documents/GitHub/Pseudo_fluor/")

#load genome table
gnm.tbl <- read.delim("C:/Users/patty/OneDrive/Documents/GitHub/Pseudo_fluor/PFC_genome_table.txt")

#create new folder for BLAST results and faa files
dir.create("blast")
dir.create('faa')

#prep files for analysis
#this takes the name of the genome and adds them to every sequence and then adds a sequence name to each (and to file name)
for(i in 1:nrow(gnm.tbl)){
  panPrep(file.path("C:/Users/patty/OneDrive/Documents/GitHub/Pseudo_fluor/PFC_aggregated_prot/", str_c(gnm.tbl$File[i], ".faa")),
          gnm.tbl$genome_id[i],
          file.path("faa", str_c(gnm.tbl$genome_id[i], ".faa")))
}

#read in protein files and BLASTp them
faa.files<-list.files("C:/Users/patty/OneDrive/Documents/GitHub/Pseudo_fluor/faa", pattern = "\\.faa$", full.names = T)
blastpAllAll(faa.files, out.folder = "blast", verbose=T, threads=2)

#get list of BLAST files
blast.files<-list.files("C:/Users/patty/OneDrive/Documents/GitHub/Pseudo_fluor/blast/", pattern = "txt$", full.names = T)
#get distances
dst.tbl <- bDist(blast.files = blast.files, e.value=0.05, verbose=T)

#visualize distances (optional)
ggplot(dst.tbl) +
  geom_histogram(aes(x = Distance), bins = 100)

#hierarchical cluster data
clst.blast <- bClust(dst.tbl, linkage = "complete", threshold = 0.75)
clst.blast2<-as.data.frame(clst.blast)

#construct pangenome matrix
panmat.blast <- panMatrix(clst.blast)

#visualize the number of clusters per genomes
tibble(Clusters = as.integer(table(factor(colSums(panmat.blast > 0),
                                          levels = 1:nrow(panmat.blast)))),
       Genomes = 1:nrow(panmat.blast)) %>% 
  ggplot(aes(x = Genomes, y = Clusters)) +
  geom_col() +
  ylab("Genes")+
  xlab("Number of Genomes")

#calculate pangenome size
heaps.est <- heaps(panmat.blast, n.perm = 500)
print(heaps.est)
# Intercept      alpha 
#811.142958   1.141321
#alpha is above 1, so pangenome is closed
print(chao(panmat.blast))
#9302

fitted <- binomixEstimate(panmat.blast, K.range = 3:8)
print(fitted$BIC.tbl)

ncomp <- 4
fitted$Mix.tbl %>% 
  filter(Components == ncomp) %>% 
  mutate(Single = Mixing.proportion * Detection.prob) %>%
  ggplot() +
  geom_col(aes(x = "", y = Single, fill = Detection.prob)) +
  coord_polar(theta = "y") +
  labs(x = "", y = "", title = "Average genome gene-family distribution",
       fill = "Detection\nprobability") +
  scale_fill_gradientn(colors = c("pink", "orange", "green", "cyan", "blue"))

#view clustered genomes (manhattan unweighted)
library(ggdendro)
d.man <- distManhattan(panmat.blast)
ggdendrogram(dendro_data(hclust(d.man, method = "average")),
             rotate = T, theme_dendro = FALSE) +
  labs(x = "", y = "Unweighted Manhattan distance")

#view clustered genomes (manhattan weighted)
pm <- panmat.blast                                                   # make a copy
rownames(pm) <- gnm.tbl$GenBank_ID[match(rownames(pm), gnm.tbl$genome_id)] # new rownames
weights <- geneWeights(pm, type = "shell")
dendrogram<-distManhattan(pm, weights = weights) %>% 
  hclust(method = "average") %>% 
  dendro_data() %>% 
  ggdendrogram(rotate = TRUE, theme_dendro = FALSE) +
  labs(x = '', y = "Weighted Manhattan distance")


#rarefaction of pangenome
rar.tbl <- rarefaction(panmat.blast, n.perm = 10000)

#plot all curves
library(tidyverse)
rar.tbl %>% 
  gather(key = "Permutation", value = "Clusters", -Genome) %>% 
  ggplot(aes(x = Number of Genomes, y = Genes, group = Permutation)) +
  geom_line()

#aggregate data
rar.tbl2<-rar.tbl %>% 
  gather(key = "Permutation", value = "Clusters", -Genome)

#plot
ggplot(rar.tbl2, aes(Genome, Clusters))+
  geom_point()+
  ylab("Number of genes")+
  ylab("Number of genomes")+
  geom_smooth()+
  theme_bw()


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
#0.07105227 0.03486338
#For example, a genomic fluidity of 0.1 represents that a pair of genomes have on average 10% unique genes and share 90% of their genes.

#calculate gene weights
weighty<-geneWeights(panmat.blast)
hist(weighty)



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
