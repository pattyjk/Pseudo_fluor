###############Viral analysis
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
viral_sum<-ddply(viral_data, c('Genome', 'Pond', 'Hit2'), summarize, abun=length(Hit))

#plot histogram of blast percent identity hits
ggplot(viral_data, aes(Per_ident, fill=Pond))+
  geom_histogram()+
  facet_wrap(~Hit2)+
  scale_fill_manual(values=c("grey", 'black'))+
  ylab("Number of BLASTn hits")+
  xlab("Percent Identity")+
  theme_bw()

#plot data
ggplot(viral_sum, aes(Genome, abun, fill=Hit2))+
  facet_wrap(~Pond, scales='free_y')+
  theme_bw()+
  ylab("Number of viral fragments")+
  xlab("")+
  scale_fill_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
  geom_bar(stat='identity')+
  coord_flip()

#calculate average umber of viral hits per genome
library(plyr)
viral_sum2<-ddply(viral_sum, c("Genome", "Pond"), summarize, total=sum(abun), mean=mean(abun), n=length(abun))
viral_sum3<-ddply(viral_sum, c("Pond"), summarize, total=sum(abun), mean=mean(abun), n=length(abun), sd=sd(abun))

#calcualte average identity for each genome
mean(viral_data$Per_ident)
#[1] 86.30184
sd(viral_data$Per_ident)
#7.635666

viral_ident<-ddply(viral_data, c("Genome", "Pond"), summarize, mean=mean(Per_ident), sd=sd(Per_ident))
