##Look at ANI
#read in dataframe from fastANI
PF_cali_ani_output <- read.delim("Pseudo_fluor/PFC_ANI.txt")
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
################Fisher exact test to determine genes that vary between core/shell/cloud
###Read in panaroo pan matrix
library(readr)
library(ggplot2)
gene_presence_absence <- read_delim("Pseudo_fluor/PFC_panaroo_results/gene_presence_absence.Rtab", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
dim(gene_presence_absence)
#[1] 7779   20

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
l3_counts<-read.delim('Pseudo_fluor/fisher_test_results.txt')

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
  ylab("P value")+
  xlab("")+
  theme_bw()+
  coord_flip()

ggplot(l3_counts, aes(Level3, lor_shell_unique, color=p_value_shell_unique))+
  geom_point()+
  ylab("Log Odds Ratio")+
  xlab("")+
  theme_bw()+
  coord_flip()

#reshape data to plot
library(reshape2)
l3_melt<-melt(l3_counts[,c(1,3,8)])


ggplot(l3_melt, aes(Level3, value, fill=variable))+
  geom_bar(stat='identity')+
  theme_light()+
  ylab("Number of genes (Log10)")+
  xlab("KEGG Orthology- Level 3")+
  coord_flip()+
  scale_y_log10()

#remove pathways not significantly different
l3_sig<-l3_counts[which(l3_counts$p_value_core_shell<0.05),]

#get percentage of genes in these categories
l3_sig$per_cor<-100*(l3_sig$no_genes_core/sum(l3_counts$total))
l3_sig$per_shell<-100*(l3_sig$no_genes_core/sum(l3_counts$shellunique))

####################################################
########Fisher exact test to determine differences between sites
####################################################
##########################
###Read in panaroo pan matrix
library(readr)
library(ggplot2)
library(plyr)
library(dplyr)
gene_presence_absence <- read_delim("Pseudo_fluor/PFC_panaroo_results/gene_presence_absence.Rtab", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
dim(gene_presence_absence)
#[1] 7779   20

#read in KEGG annotations of genes
kegg_annotations <- read.delim("~/GitHub/Pseudo_fluor/kegg_annotations.txt")
kegg_annotations$KO<-trimws(kegg_annotations$KO)

#read in KEGG DB
ko_db<-read.delim("Pseudo_fluor/full_kegg.txt", header=T)
ko_db$KO<-trimws(ko_db$KO)

#create table of Sixty and Conness genes
gene_presence_absence$sum<-rowSums(gene_presence_absence[,-1])
conness_genes<-select(gene_presence_absence, contains(c("CP", "Gene")))
sixty_genes<-select(gene_presence_absence, -contains("CP"))
sixty_genes$Type<-"Sixty Lake"
conness_genes$Type<-'Conness Pond'

#see how many genes each has
dim(sixty_genes)
#[1] 7779   12

dim(conness_genes)
#[1] 7779   12

#remove genes that are not annotated in KEGG
sixty_genes<-sixty_genes[-grep("group_", sixty_genes$Gene),]
conness_genes<-conness_genes[-grep("group_", conness_genes$Gene),]

#see how many genes each has after removing non annotated genes
dim(sixty_genes)
#[1] 3427   12

dim(conness_genes)
#[1] 3427   12

#reshape data
library(reshape2)
conness_m<-melt(conness_genes[,-12])
sixty_m<-melt(sixty_genes[,-12])

#append KO information
conness_m<-merge(conness_m, kegg_annotations, by='Gene', all.y=F)
sixty_m<-merge(sixty_m, kegg_annotations, by='Gene', all.y=F)

#append KO DB annotations
conness_m<-merge(conness_m, unique(ko_db), by.x='KO', by.y='KO', all.y=F, all.x=T)
sixty_m<-merge(sixty_m, unique(ko_db), by.x='KO', by.y='KO', all.y=F, all.x=T)

#calculate L3 number of genes for each
conness_m_count_l3<-ddply(conness_m, c("Level3"), summarize, no_conness=length(unique(Gene.x)))
sixty_m_count_l3<-ddply(sixty_m, c("Level3"), summarize, no_sixty=length(unique(Gene.x)))

#merge datasets
l3_counts<-merge(sixty_m_count_l3, conness_m_count_l3, by='Level3')

#reshape data to plot
l3_m<-melt(l3_counts)
ggplot(l3_m, aes(Level3, value, fill=variable))+
  geom_bar(stat='identity')+
  theme_bw()+
  ylab("Number of Genes")+
  coord_flip()+
  xlab("")

#calculate total genes 
l3_counts$total<-l3_counts$no_conness+l3_counts$no_sixty

#calculate the number of genes not found in each
l3_counts$not_sixty<-sum(l3_counts$no_sixty)-l3_counts$no_sixty
l3_counts$not_conness<-sum(l3_counts$no_conness)-l3_counts$no_conness

#for each L3 KEGG category, calculate Fisher Exact Test for comparisons between shell, core, and unique categories
# Create an empty vector to store the p-values and LOR
p_values_ponds <- numeric(nrow(l3_counts))
lor_ponds <- numeric(nrow(l3_counts))

# Loop through each row and apply Fisher exact test to compare ponds
set.seed(505)
for (i in 1:nrow(l3_counts)) {
  contingency_table <- matrix(c(l3_counts[i, "no_sixty" ], l3_counts[i, "no_conness" ],
                                l3_counts[i, "not_sixty"], l3_counts[i, "not_conness"]),
                              nrow = 2)
  
  
  # Apply Fisher exact test
  fisher_result <- fisher.test(contingency_table)
  
  # Store the p-value in the vector
  p_values_ponds[i] <- fisher_result$p.value
  lor_ponds[i] <- fisher_result$estimate
}

# Add the p-values/lor as a new column the original L3 data frame
l3_counts$p_value <- p_values_ponds
l3_counts$lor <- lor_ponds


#correct p-values with Benjamini-Hochberg correction
l3_counts$p_value_corrected <- p.adjust(l3_counts$p_value, method = "BH")

#see what's different
which(l3_counts$p_value_corrected<0.05)
which(l3_counts$p_value<0.05)
#none significant at any KEGG level :(

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

#read in shell break down of Conness Pond gene shells
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


#append Bd inhibition data to antismash data
bd_data<-read.delim("~/Github/Pseudo_fluor/bd_inhibition_data.txt")

#remove controls from data (blanks have been subtracted out already)
bd_data<-bd_data[-which(bd_data$Group =='Heat killed' | bd_data$Group =='Tryptone' | bd_data$Group =='Bd_Prova' | bd_data$Group =='Bd_Tryptone'),]

#calculate total number of BCG per genome
anti_sum2<-ddply(anti_sum, c("Genome"), summarize, n_clust=length(no_genes))

#calculate summary stats for BD inhibition data
bd_sum<-ddply(bd_data, c("Bd_strain", "Group"), summarize, mean=mean(Percent_inhibition), sd=sd(Percent_inhibition), n=length(Percent_inhibition), se=sd/n)
anti_sum_bd<-merge(anti_sum2, bd_sum, by.x='Genome', by.y='Group')

#plot clusters and Bd data
library(ggpubr)
ggplot(anti_sum_bd, aes(n_clust, mean, color=Bd_strain))+
  geom_point()+
  ylab("Percent Bd inhibition")+
  xlab("Number of BGC")+
  geom_smooth(method='lm')+
  stat_cor(method='spearman')

##################analayze Bd inhibition data
#load in Bd-inhibition data
bd_data<-read.delim("~/Github/Pseudo_fluor/bd_inhibition_data.txt")

#remove controls from data (blanks have been subtracted out already)
bd_data<-bd_data[-which(bd_data$Group =='Heat killed' | bd_data$Group =='Tryptone' | bd_data$Group =='Bd_Prova' | bd_data$Group =='Bd_Tryptone'),]

#correct p-values
bd_data$p.adj<-p.adjust(bd_data$P_value, method = 'hochberg')

#write table with corrected p-values
write.table(bd_data, 'Pseudo_fluor/bd_inhib_data_corrcted_p.txt', quote=F, row.names = F, sep='\t')
bd_data<-read.delim('Pseudo_fluor/bd_inhib_data_corrcted_p.txt', header=T)

#plot data
library(ggplot2)
ggplot(bd_data, aes(Group, Percent_inhibition, fill=Bd_strain))+
  geom_boxplot()+
  coord_flip()+
  xlab("PF strain")+
  ylab("Percent inhibition of Bd")

ggplot(bd_data, aes(Percent_inhibition, fill=Pond))+
  geom_histogram()+
  scale_fill_manual(values=c('grey', 'red'))+
  facet_wrap(~Bd_strain)+
  ylab("Frequency")+
  xlab("Percent Bd inhibition")


##Metabolic modeling data based on KGGG annotation in ANVIO (v7.1)

##ANVIO 7.1 analysis
#annotate with KEGG/diamond
for i in *.db
do
anvi-run-kegg-kofams -c $i -T 4
done

#estimate metabolism for each genome
for i in *.db
do
anvi-estimate-metabolism -c $i -O $i_metabolism --kegg-output-modes kofam_hits
done

#change to R analysis
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

#####Virus
#read in viral elements blast results
# Set the directory path
directory_path <- "~/Github/Pseudo_fluor/viral_genomes/blast_out/"

# Get a list of files in the directory
files <- list.files(directory_path, full.names = TRUE)

# Initialize an empty data frame to store the combined data
viral_data <- data.frame()

# Loop through each file and add it to a data frame, creating a column with the file name
for (file in files) {
  # Read the data from the file
  current_data <- read.table(file, header = F, sep = "\t")  # Adjust parameters based on your file format
  
  # Add a new column with the file name
  current_data$genome <- basename(file)
  
  # Combine the current data with the existing data
  viral_data <- rbind(viral_data, current_data)
}

#remove stuff that's not of interest and add relavant column names
viral_data<-viral_data[,c(2:5, 12)]

#how many viral hits?
dim(viral_data)
#[1] 387   5

#write to file
write.table(viral_data, 'Pseudo_fluor/viral_data.txt', sep='\t', quote=F, row.names=F)

#reread in data
viral_data<-read.delim("Pseudo_fluor/viral_data.txt", header=T)

#summarize data
library(plyr)
library(ggplot2)
viral_sum<-ddply(viral_data, c('Genome', 'Pond', 'Hit'), summarize, abun=length(Hit))

#plot histogram of blast percent identity hits
ggplot(viral_data, aes(Per_ident, fill=Pond))+
  geom_histogram()+
  facet_wrap(~Hit)+
  ylab("Number of BLASTn hits")+
  xlab("Percent Identity")+
  theme_bw()

#plot data
ggplot(viral_sum, aes(Genome, abun, fill=Hit))+
  facet_wrap(~Pond, scales='free_y')+
  theme_bw()+
  scale_fill_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
  geom_bar(stat='identity')+
  coord_flip()

###############Look at AmphiBac database
#read in data & library
amphibac<-read.delim("Pseudo_fluor/amphibac_database.txt", header=T)
library(ggplot2)
library(tidyverse)

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

#look at own data
bd_data<-read.delim('Pseudo_fluor/bd_inhib_data_corrcted_p.txt', header=T)

#separate to 197/423
bd_197<-bd_data[which(bd_data$Bd_strain=='Jel197'),]
bd_423<-bd_data[which(bd_data$Bd_strain=='JEL423'),]

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

length(which(bd_summary$mean <50))
#6
length(which(bd_summary$mean >50))
#16

#use Fisher exact test to determine categories of interest for testing
###Read in panaroo pan matrix
library(readr)
library(ggplot2)
library(plyr)
library(dplyr)
gene_presence_absence <- read_delim("Pseudo_fluor/Pangenome_results/PFC_panaroo_results/gene_presence_absence.Rtab", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
dim(gene_presence_absence)
#[1] 7779   20

#read in KEGG annotations of genes
kegg_annotations <- read.delim("~/GitHub/Pseudo_fluor/kegg_annotations.txt")
kegg_annotations$KO<-trimws(kegg_annotations$KO)

#read in KEGG DB
ko_db<-read.delim("Pseudo_fluor/full_kegg.txt", header=T)
ko_db$KO<-trimws(ko_db$KO)

#create table of weakly and strongly inhibitory
weak_genes<-gene_presence_absence[,c(1,4,5,8,10,12,15,17,18)]
strong_genes<-gene_presence_absence[,c(1,2,3,6,7,9,11,13,14,16,19,20)]
weak_genes$Type<-'Weakly inhibitory'
strong_genes$Type<-'Inhibitory'

#remove genes that are not annotated in KEGG
#weak_genes<-weak_genes[-grep("group_", weak_genes$Gene),]
#strong_genes<-strong_genes[-grep("group_", strong_genes$Gene),]

#reshape data
library(reshape2)
weak_m<-melt(weak_genes)
strong_m<-melt(strong_genes)

#append KO information
weak_m<-merge(weak_m, kegg_annotations, by='Gene', all.y=F)
strong_m<-merge(strong_m, kegg_annotations, by='Gene', all.y=F)

#append KO DB annotations
weak_m<-merge(weak_m, unique(ko_db), by.x='KO', by.y='KO', all.y=F, all.x=T)
strong_m<-merge(strong_m, unique(ko_db), by.x='KO', by.y='KO', all.y=F, all.x=T)

#calculate L3 number of genes for each
weak_m_count_l3<-ddply(weak_m, c("Level3"), summarize, no_weak=length(Gene.x))
strong_m_count_l3<-ddply(strong_m, c("Level3"), summarize, no_strong=length(Gene.x))

#merge datasets
l3_counts<-merge(strong_m_count_l3, weak_m_count_l3, by='Level3')

#calculate total genes 
l3_counts$total<-l3_counts$no_strong+l3_counts$no_weak

#calculate the number of genes not found in each
l3_counts$not_strong<-abs(l3_counts$no_weak-l3_counts$no_strong)
l3_counts$not_weak<-abs(l3_counts$no_strong-l3_counts$no_weak)

#for each L3 KEGG category, calculate Fisher Exact Test for comparisons between shell, core, and unique categories
# Create an empty vector to store the p-values and LOR
p_values_ponds <- numeric(nrow(l3_counts))
lor_ponds <- numeric(nrow(l3_counts))

# Loop through each row and apply Fisher exact test to compare ponds
set.seed(505)
for (i in 1:nrow(l3_counts)) {
  contingency_table <- matrix(c(l3_counts[i, "no_weak" ], l3_counts[i, "no_strong" ],
                                l3_counts[i, "not_weak"], l3_counts[i, "not_strong"]),
                              nrow = 2)

  
  # Apply Fisher exact test
  fisher_result <- fisher.test(contingency_table)
  
  # Store the p-value in the vector
  p_values_ponds[i] <- fisher_result$p.value
  lor_ponds[i] <- fisher_result$estimate
}

# Add the p-values/lor as a new column the original L3 data frame
l3_counts$p_value <- p_values_ponds
l3_counts$lor <- lor_ponds

#correct p-values with Benjamini-Hochberg correction
l3_counts$p_value_corrected <- p.adjust(l3_counts$p_value, method = "BH")

#see what's different
length(which(l3_counts$p_value_corrected<0.05))
#17

#create significant table
l3_sig<-l3_counts[which(l3_counts$p_value_corrected<0.05),]

ggplot(l3_sig, aes(Level3, ))
