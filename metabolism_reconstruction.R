################Metabolic reconstruction based on KEGG
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
