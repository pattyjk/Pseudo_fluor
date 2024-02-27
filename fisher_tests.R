###########Fisher tests
##########################
################Fisher exact test to determine genes that vary between core/shell/cloud
###Read in panaroo pan matrix
library(readr)
library(ggplot2)
gene_presence_absence <- read_delim("Pseudo_fluor/Pangenome_results/PFC_panaroo_results/struct_presence_absence.Rtab", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
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
set.seed(515)
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
gene_presence_absence <- read_delim("Pseudo_fluor/Pangenome_results/PFC_panaroo_results/gene_presence_absence.Rtab", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
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

#remove genes that are not annotated in KEGG
#sixty_genes<-sixty_genes[-grep("group_", sixty_genes$Gene),]
#conness_genes<-conness_genes[-grep("group_", conness_genes$Gene),]

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
l3_counts$not_sixty<-abs(l3_counts$no_sixty-l3_counts$no_conness)
l3_counts$not_conness<-abs(l3_counts$no_conness-l3_counts$no_sixty)

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

#write table
write.table(l3_sig, "Pseudo_fluor/l3_sig_bd_inhibition.txt", sep='\t', row.names=F, quote=F)