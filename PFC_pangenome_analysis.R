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

#construct pangenome matrix
panmat.blast <- panMatrix(clst.blast)

#visualize the number of clusters per genomes
tibble(Clusters = as.integer(table(factor(colSums(panmat.blast > 0),
                                          levels = 1:nrow(panmat.blast)))),
       Genomes = 1:nrow(panmat.blast)) %>% 
  ggplot(aes(x = Genomes, y = Clusters)) +
  geom_col() + labs(title = "Number of clusters found in 1, 2,...,all genomes")

#calculate pangenome size
heaps.est <- heaps(panmat.blast, n.perm = 500)
print(heaps.est)
# Intercept      alpha 
#811.142958   1.141321
print(chao(panmat.blast))
#9302

fitted <- binomixEstimate(panmat.blast, K.range = 3:8)
print(fitted$BIC.tbl)

ncomp <- 4
fitted$Mix.tbl %>% 
  filter(Components == ncomp) %>% 
  ggplot() +
  geom_col(aes(x = "", y = Mixing.proportion, fill = Detection.prob)) +
  coord_polar(theta = "y") +
  labs(x = "", y = "", title = "Pan-genome gene family distribution",
       fill = "Detection\nprobability") +
  scale_fill_gradientn(colors = c("pink", "orange", "green", "cyan", "blue"))

#view clustered genomes (manhattan unweighted)
library(ggdendro)
d.man <- distManhattan(panmat.blast)
ggdendrogram(dendro_data(hclust(d.man, method = "average")),
             rotate = TRUE, theme_dendro = FALSE) +
  labs(x = "Genomes", y = "Unweighted Manhattan distance", title = "Pan-genome dendrogram")

#view clustered genomes (manhattan weighted)
pm <- panmat.blast                                                   # make a copy
rownames(pm) <- gnm.tbl$File[match(rownames(pm), gnm.tbl$genome_id)] # new rownames
weights <- geneWeights(pm, type = "shell")
distManhattan(pm, weights = weights) %>% 
  hclust(method = "average") %>% 
  dendro_data() %>% 
  ggdendrogram(rotate = TRUE, theme_dendro = FALSE) +
  labs(x = "Genomes", y = "Weighted Manhattan distance", title = "Pan-genome dendrogram")


#rarefaction of pangenome
rar.tbl <- rarefaction(panmat.blast, n.perm = 10000)

#plot it
rar.tbl %>% 
  gather(key = "Permutation", value = "Clusters", -Genome) %>% 
  ggplot(aes(x = Genome, y = Clusters, group = Permutation)) +
  geom_line()

###Read in panaroo pan matrix
library(readr)
gene_presence_absence <- read_delim("PFC_panaroo_results/gene_presence_absence.Rtab", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

colSums(gene_presence_absence[,-1])
