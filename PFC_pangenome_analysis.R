##Look at ANI
#read in dataframe from fastANI
setwd("Pseudo_fluor/")
PF_cali_ani_output <- read.delim("PFC_ANI.txt")
PF_cali_ani_output <- PF_cali_ani_output[-which(PF_cali_ani_output$Sp1 == '630A' | PF_cali_ani_output$Sp1 == '629A2A' | PF_cali_ani_output$Sp2 == '630A' | PF_cali_ani_output$Sp2 == '629A2A'),]
#load ggplot2
library(ggplot2)

#heatmap of ANI
ggplot(PF_cali_ani_output, aes(Sp1, Sp2, fill=ANI))+
  geom_tile()+
  theme_bw()+
  ylab("")+
  xlab("")

#boxplot
ggplot(PF_cali_ani_output, aes(Comparison, ANI, fill=Lake1))+
  geom_boxplot()+
  theme_bw()+
  ylab("Average Nucleotide Identity")+
  xlab("")

t.test(PF_cali_ani_output$ANI ~ PF_cali_ani_output$Comparison)
#t = -5.8484, df = 347.07, p-value = 1.147e-08

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
panmat.blast2<-as.data.frame(panmat.blast)

panmat.blast3<-as.matrix(panmat.blast)

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
             rotate = TRUE, theme_dendro = FALSE) +
  labs(x = "Genomes", y = "Unweighted Manhattan distance")

#view clustered genomes (manhattan weighted)
pm <- panmat.blast                                                   # make a copy
rownames(pm) <- gnm.tbl$File[match(rownames(pm), gnm.tbl$genome_id)] # new rownames
weights <- geneWeights(pm, type = "shell")
distManhattan(pm, weights = weights) %>% 
  hclust(method = "average") %>% 
  dendro_data() %>% 
  ggdendrogram(rotate = TRUE, theme_dendro = FALSE) +
  labs(x = "Genomes", y = "Weighted Manhattan distance")


#rarefaction of pangenome
rar.tbl <- rarefaction(panmat.blast, n.perm = 10000)

#plot it
rar.tbl %>% 
  gather(key = "Permutation", value = "Clusters", -Genome) %>% 
  ggplot(aes(x = Genome, y = Clusters, group = Permutation)) +
  geom_line()

#calculate average pangenome fluidity
pg_fluid<-fluidity(panmat.blast, n.sim = 10000)
#Mean        Std 
#0.07105227 0.03486338
#For example, a genomic fluidity of 0.1 represents that a pair of genomes have on average 10% unique genes and share 90% of their genes.

#calculate gene weights
weighty<-geneWeights(panmat.blast)
hist(weighty)

##########################
###Read in panaroo pan matrix
library(readr)
gene_presence_absence <- read_delim("Pseudo_fluor/PFC_panaroo_results/gene_presence_absence.Rtab", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

#read in KEGG annotations of genes
kegg_annotations <- read.delim("~/GitHub/Pseudo_fluor/kegg_annotations.txt")
kegg_annotations$KO<-trimws(kegg_annotations$KO)

#read in KEGG DB
ko_db<-read.delim("Pseudo_fluor/full_kegg.txt", header=T)
ko_db$KO<-trimws(ko_db$KO)

#create table of core genes (genes >19 per row)
gene_presence_absence$sum<-rowSums(gene_presence_absence[,-1])
core_genes<-gene_presence_absence[-which(gene_presence_absence$sum<=18),]
shell_genes<-gene_presence_absence[which(gene_presence_absence$sum<=18),]
shell_genes$Type<-"Shell"
core_genes$Type<-"Core"


#reshape data
library(dplyr)
library(plyr)
library(reshape2)
core_m<-melt(core_genes[,-21])
shell_m<-melt(shell_genes[,-21])

#append KO information
core_m<-merge(core_m, kegg_annotations, by='Gene', all.y=F)
shell_m<-merge(shell_m, kegg_annotations, by='Gene', all.y=F)

#append KO DB annotations
core_m<-merge(core_m, unique(ko_db), by.x='KO', by.y='KO', all.y=F, all.x=T)
shell_m<-merge(shell_m, unique(ko_db), by.x='KO', by.y='KO', all.y=F, all.x=T)

#summarize Level 3 kegg for each
core_m<-ddply(core_m, c("Level2", 'Type'), summarize, sum=sum(value))
names(core_m)<-c("Level3", "Type", "Total_core")
shell_m<-ddply(shell_m, c("Level2", 'Type'), summarize, sum=sum(value))
names(shell_m)<-c("Level3", "Type", "Shell_core")

#merge data sets
gene_sets<-merge(core_m, shell_m, by='Level3')
gene_sets$Total<-gene_sets$Total_core+gene_sets$Shell_core

ggplot(gene_sets, aes(Level3, Total_core))+
  geom_bar(stat='identity')+
  scale_y_log10()+
  coord_flip()
ggplot(gene_sets, aes(Level3, Shell_core))+
  geom_bar(stat='identity')+
  scale_y_log10()+
  coord_flip()

##########################################
#read in break down of gene shells
summary_statistics <- read.delim("~/GitHub/Pseudo_fluor/PFC_panaroo_results/summary_statistics.txt", header=FALSE)
ggplot(summary_statistics[-5,], aes(x="", y=V3, fill=V2))+
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)


