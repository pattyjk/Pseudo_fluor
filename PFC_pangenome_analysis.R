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
dim(gene_presence_absence)
#[1] 7779   21

#calculate genes for shells
.15*20
#[1] 3

.95*20
#[1] 19

#read in KEGG annotations of genes
kegg_annotations <- read.delim("~/GitHub/Pseudo_fluor/kegg_annotations.txt")
kegg_annotations$KO<-trimws(kegg_annotations$KO)

#read in KEGG DB
ko_db<-read.delim("Pseudo_fluor/full_kegg.txt", header=T)
ko_db$KO<-trimws(ko_db$KO)

#create table of core genes (genes >19 per row), shell genes (>15%-95%), unique genes (<15%)
gene_presence_absence$sum<-rowSums(gene_presence_absence[,-1])
core_genes<-gene_presence_absence[which(gene_presence_absence$sum>=19),]
shell_genes<-gene_presence_absence[which(gene_presence_absence$sum<19 & gene_presence_absence$sum>3),]
unique_genes<-gene_presence_absence[which(gene_presence_absence$sum<=3),]
shell_genes$Type<-"Shell"
core_genes$Type<-"Core"
unique_genes$Type<-'Unique'

#see how many genes each has
dim(core_genes)
#[1] 5147   22

dim(shell_genes)
#[1] 1770   22

dim(unique_genes)
#[1] 862  22

#remove genes that are not annotated in KEGG
shell_genes<-shell_genes[-grep("group_", shell_genes$Gene),]
core_genes<-core_genes[-grep("group_", core_genes$Gene),]
unique_genes<-unique_genes[-grep("group_", unique_genes$Gene),]

#see how many genes each has after removing non annotated genes
dim(core_genes)
#[1] 3080   22

dim(shell_genes)
#[1] 268  22

dim(unique_genes)
#[1] 79 22

#reshape data
library(dplyr)
library(plyr)
library(reshape2)
core_m<-melt(core_genes[,-21])
shell_m<-melt(shell_genes[,-21])
unique_m<-melt(unique_genes[,-21])

#append KO information
core_m<-merge(core_m, kegg_annotations, by='Gene', all.y=F)
shell_m<-merge(shell_m, kegg_annotations, by='Gene', all.y=F)
unique_m<-merge(unique_m, kegg_annotations, by='Gene', all.y=F)

#append KO DB annotations
core_m<-merge(core_m, unique(ko_db), by.x='KO', by.y='KO', all.y=F, all.x=T)
shell_m<-merge(shell_m, unique(ko_db), by.x='KO', by.y='KO', all.y=F, all.x=T)
unique_m<-merge(unique_m, unique(ko_db), by.x='KO', by.y='KO', all.y=F, all.x=T)

#calculate L3 number of genes for each
core_m_count_l3<-ddply(core_m, c("Level3"), summarize, no_genes_core=length(unique(Gene.x)))
shell_m_count_l3<-ddply(shell_m, c("Level3"), summarize, no_genes_shell=length(unique(Gene.x)))
unique_m_count_l3<-ddply(unique_m, c("Level3"), summarize, no_genes_unique=length(unique(Gene.x)))

#merge datasets
l3_counts<-merge(shell_m_count_l3, core_m_count_l3, by='Level3')
l3_counts<-merge(l3_counts, unique_m_count_l3, by='Level3')

#reshape data to plot
l3_m<-melt(l3_counts)
ggplot(l3_m, aes(Level3, value, fill=variable))+
  geom_bar(stat='identity')+
  theme_bw()+
  ylab("Number of Genes")+
  coord_flip()+
  xlab("")

#calculate total genes for each
l3_counts$total<-l3_counts$no_genes_shell+l3_counts$no_genes_core+l3_counts$no_genes_unique
l3_counts$coreshell<-l3_counts$no_genes_shell+l3_counts$no_genes_core
l3_counts$coreunique<-l3_counts$no_genes_core+l3_counts$no_genes_unique
l3_counts$shellunique<-l3_counts$no_genes_shell+l3_counts$no_genes_unique

#calculate the number of genes not found in each
l3_counts$not_shell<-sum(l3_counts$no_genes_shell)-l3_counts$no_genes_shell
l3_counts$not_core<-sum(l3_counts$no_genes_core)-l3_counts$no_genes_core
l3_counts$not_unique<-sum(l3_counts$no_genes_unique)-l3_counts$no_genes_unique

#for each L3 KEGG category, calculate Fisher Exact Test for comparisons between shell, core, and unique categories
# Create an empty vector to store the p-values and LOR
p_values_coreunique <- numeric(nrow(l3_counts))
p_values_coreshell <- numeric(nrow(l3_counts))
p_values_shellunique <- numeric(nrow(l3_counts))

lor_coreunique <- numeric(nrow(l3_counts))
lor_coreshell <- numeric(nrow(l3_counts))
lor_shellunique <- numeric(nrow(l3_counts))

# Loop through each row and apply Fisher exact test to compare core/unique
set.seed(505)
for (i in 1:nrow(l3_counts)) {
  contingency_table <- matrix(c(l3_counts[i, "no_genes_core"], l3_counts[i, "no_genes_unique"],
                                l3_counts[i, "not_core"], l3_counts[i, "not_unique"]),
                              nrow = 2)
  
  # Apply Fisher exact test
  fisher_result <- fisher.test(contingency_table)
  
  # Store the p-value in the vector
  p_values_coreunique[i] <- fisher_result$p.value
  lor_coreunique[i] <- fisher_result$estimate
}

# Loop through each row and apply Fisher exact test to compare core/shell
for (i in 1:nrow(l3_counts)) {
  contingency_table2 <- matrix(c(l3_counts[i, "no_genes_core"], l3_counts[i, "no_genes_shell"],
                                l3_counts[i, "not_core"], l3_counts[i, "not_shell"]),
                              nrow = 2)
  
  # Apply Fisher exact test
  fisher_result2 <- fisher.test(contingency_table2)
  
  # Store the p-value in the vector
  p_values_coreshell[i] <- fisher_result2$p.value
  lor_coreshell[i] <- fisher_result2$estimate
}

# Loop through each row and apply Fisher exact test to compare unique/shell
for (i in 1:nrow(l3_counts)) {
  contingency_table3 <- matrix(c(l3_counts[i, "no_genes_shell"], l3_counts[i, "no_genes_unique"],
                                l3_counts[i, "not_shell"], l3_counts[i, "not_unique"]),
                              nrow = 2)
  
  # Apply Fisher exact test
  fisher_result3 <- fisher.test(contingency_table3)
  
  # Store the p-value in the vector
  p_values_shellunique[i] <- fisher_result3$p.value
  lor_shellunique[i] <- fisher_result3$estimate
}

# Add the p-values/lor as a new column the original L3 data frame
l3_counts$p_value_core_unique <- p_values_coreunique
l3_counts$p_value_core_shell <- p_values_coreshell
l3_counts$p_value_shell_unique <- p_values_shellunique

l3_counts$lor_core_unique <- lor_coreunique
l3_counts$lor_core_shell <- lor_coreshell
l3_counts$lor_shell_unique <- lor_shellunique

#correct p-values with Benjamini-Hochberg correction
l3_counts$p_value_core_unique <- p.adjust(l3_counts$p_value_core_unique, method = "BH")
l3_counts$p_value_core_shell <- p.adjust(l3_counts$p_value_core_shell, method = "BH")
l3_counts$p_value_shell_unique <- p.adjust(l3_counts$p_value_shell_unique, method = "BH")

#write data to a file
write.table(l3_counts, 'Pseudo_fluor/fisher_test_results.txt', sep='\t', quote=F, row.names=F)

#plot log odds ratio
library(ggplot2)
ggplot(l3_counts, aes(Level3, lor_core_unique, color=p_value_core_unique))+
  geom_point()+
  ylab("Log Odds Ratio")+
  xlab("")+
  theme_bw()+
  coord_flip()

ggplot(l3_counts, aes(Level3, lor_core_shell, color=p_value_core_shell))+
  geom_point()+
  ylab("Log Odds Ratio")+
  xlab("")+
  theme_bw()+
  coord_flip()

ggplot(l3_counts, aes(Level3, lor_shell_unique, color=p_value_shell_unique))+
  geom_point()+
  ylab("Log Odds Ratio")+
  xlab("")+
  theme_bw()+
  coord_flip()


##########################################
#read in break down of gene shells for all genomes
library(ggplot2)
summary_statistics <- read.delim("~/GitHub/Pseudo_fluor/PFC_panaroo_results/summary_statistics.txt", header=FALSE)
ggplot(summary_statistics[-5,], aes(x="", y=V3, fill=V2))+
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  ggtitle("All isolates")

#read in shell break down of Sixty lake gene shells
sixty_sum_stats<-read.delim("~/GitHub/Pseudo_fluor/PFC_panaroo_results_sixty/summary_statistics.txt", header=F)
ggplot(sixty_sum_stats[-5,], aes(x="", y=V3, fill=V2))+
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  ggtitle("Sixty Lake Isolates")

#read in shell break down of Coness Pond gene shells
coness_sum_stats<-read.delim("~/GitHub/Pseudo_fluor/PFC_panaroo_results_coness/summary_statistics.txt", header=F)
ggplot(coness_sum_stats[-5,], aes(x="", y=V3, fill=V2))+
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  ggtitle("Coness Pond Isolates")


###AntiSmash
anti_out<-read.delim("~/Github/Pseudo_fluor/AntiSmash_output/antismash_out.txt")
library(dplyr)
library(plyr)

#summarize
anti_sum<-ddply(anti_out, c("Genome", "Pond", "Type"), summarize, no_genes=length(Type))
ggplot(anti_sum, aes(Genome, no_genes, fill=Type))+
  geom_bar(stat='identity')+
  theme_bw()+
  ylab("Number of Genes")+
  theme(axis.text.x = element_text( angle = 45))+
  #coord_flip()+
  scale_fill_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
  facet_wrap(~Pond, scales='free_x')

#load in Bd-inhibition data
bd_data<-read.delim("~/Github/Pseudo_fluor/bd_inhibition_data.txt")

#remove controls from data
bd_data<-bd_data[-which(bd_data$Group =='Tryptone' | bd_data$Group =='Heat killed' | bd_data$Group =='Bd_Tryptone' | bd_data$Group =='Bd_Prova'),]

#plot data
library(ggplot2)
ggplot(bd_data, aes(Group, Prop_Inhibition))+
  geom_boxplot()+
  theme_bw()+
  facet_wrap(~Bd_strain)+
  ylab("Percent Bd Inhibition")+
  xlab("")+
  coord_flip()


##Metabolic modeling data based on KGGG annotation in ANVIO (v7.1)

# Set the directory path
directory_path <- "~/Github/Pseudo_fluor/metabolic_reconstruction/"

# Get a list of files in the directory
files <- list.files(directory_path, full.names = TRUE)

# Initialize an empty data frame to store the combined data
metabolic_data <- data.frame()

# Loop through each file
for (file in files) {
  # Read the data from the file
  current_data <- read.table(file, header = TRUE, sep = "\t")  # Adjust parameters based on your file format
  
  # Add a new column with the file name
  current_data$genome <- basename(file)
  
  # Combine the current data with the existing data
  metabolic_data <- rbind(metabolic_data, current_data)
}

#strip extra junk off genome names
metabolic_data$genome<-gsub("_scaffolds.fasta.db_modules.txt", "",metabolic_data$genome)
metabolic_data$genome<-gsub("meta_", "",metabolic_data$genome)

#write catenated metabolisms to a file
write.table(metabolic_data, '~/GitHub/Pseudo_fluor/catenated_metabolic_data.txt', sep='\t', row.names=F, quote=F)

#summarize data
library(plyr)
metabolic_sum<-ddply(metabolic_data, c("genome", "module_name"), summarize, presence=length(module_name))

library(ggplot2)
ggplot(metabolic_sum, aes(genome,module_name, fill=presence))+
  geom_tile()


