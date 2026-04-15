## Code for processing, annotating, etc.
### Assembly
```
## Assemble genomes with SPAdes
mkdir -p ./ass

# Loop through all forward read files
for R1 in ./fastq/*-1.fq; do
    
    # Identify the reverse read by replacing -1.fq with -2.fq
    R2="${R1%-1.fq}-2.fq"
    
    # Extract the genome name (e.g., ./fastq/M1492 from ./fastq/M1492-1.fq)
    # Then use 'basename' to get just 'M1492'
    genome_name=$(basename "${R1%-1.fq}")
    
    # Define the output directory for this specific genome
    out_dir="./ass/${genome_name}"
    
    echo "Starting assembly for: $genome_name"
    echo "Using R1: $R1"
    echo "Using R2: $R2"

    # Run SPAdes
    $HOME/SPAdes-3.15.5-Linux/bin/spades.py \
        -t 64 \
        --isolate \
        -1 "$R1" \
        -2 "$R2" \
        -o "$out_dir"

    echo "Finished $genome_name. Results in $out_dir"
    echo "------------------------------------------"
done

## Aggregate and rename scaffolds by genome name

#rename the bins by folder name, assumes your in the folder where all the genomes are located
## Specify the extension 
target_extension=".fasta"

# Function to process files
process_file() {
    local filepath="$1"
    local foldername="$2"
    local filename=$(basename "$filepath")
    
    # Append folder name to the file
    new_filename="${foldername}_${filename}"
    
    # Rename the file
    mv "$filepath" "$(dirname $filepath)/$new_filename"
}

# Main loop to find and process files
find . -type f -name "*$target_extension" | while read filepath; do
    foldername=$(basename "$(dirname "$filepath")")
    process_file "$filepath" "$foldername"
done

#######find and copy all files into a new folder, assumes your are in the folder where all the genomes are located
#tell it the extension for each genome
mkdir aggregated_genomes
target_extension="_scaffolds.fasta"

# Specify the destination folder where the files will be copied
destination_folder="./aggregated_genomes"

# Create the destination folder if it doesn't exist
mkdir -p "$destination_folder"

# Main loop to find and copy files
find . -type f -name "*$target_extension" -exec cp {} "$destination_folder" \;

#remove extra files
cd aggregated_genomes/
rm *final*
rm K77*
rm  misc_broken_scaffolds.fasta
```

### Check quality with checkM
```
source activate checkm
cd /hpcstor6/scratch01/p/patrick.kearns/genomes/Assemblies_fastq_checkM/fastq/PF_Cali
mkdir checkM_results
checkm lineage_wf ./aggregated_genomes ./checkM_results -x fasta -t 48 --pplacer_threads 24 --tab_table 
conda deactivate
```
Tossed a few genomes due to bad quality (647D) and bad taxonomic alignment. 

### Annotate genomes with Prokka
```
conda activate prokka
mkdir annotations
cd aggregated_genomes
for i in *.fasta
do
prokka $i --cpus 24 --outdir $i --force --outdir /hpcstor6/scratch01/p/patrick.kearns/genomes/Assemblies_fastq_checkM/fastq/PF_Cali/annotations/$i
done
```

### Pangenome with panaroo
```
## Aggregate gff files and protein files for Pangenome analysis

## Specify the extension 
target_extension=".gff"

# Function to process files
process_file() {
    local filepath="$1"
    local foldername="$2"
    local filename=$(basename "$filepath")
    
    # Append folder name to the file
    new_filename="${foldername}_${filename}"
    
    # Rename the file
    mv "$filepath" "$(dirname $filepath)/$new_filename"
}

# Main loop to find and process files
find . -type f -name "*$target_extension" | while read filepath; do
    foldername=$(basename "$(dirname "$filepath")")
    process_file "$filepath" "$foldername"
done

#######find and copy all files into a new folder, assumes your are in the folder where all the genomes are located
#tell it the extension for each genome
target_extension=".gff"

# Specify the destination folder where the files will be copied
destination_folder="./gff_files"

# Create the destination folder if it doesn't exist
mkdir -p "$destination_folder"

# Main loop to find and copy files
find . -type f -name "*$target_extension" -exec cp {} "$destination_folder" \;


## Run Panaroo

conda activate panaroo
cd /hpcstor6/scratch01/p/patrick.kearns/genomes/Assemblies_fastq_checkM/fastq/Micro/gff_files
panaroo -i *.gff -o micro_panaroo_results --clean-mode strict --remove-invalid-genes -t 48
```
