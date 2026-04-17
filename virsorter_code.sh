#!/bin/bash

# Sample slurm submission script for the Chimera compute cluster
# Lines beginning with # are comments, and will be ignored by
# the interpreter.  Lines beginning with #SBATCH are directives
# to the scheduler.  These in turn can be commented out by
# adding a second # (e.g. ##SBATCH lines will not be processed
# by the scheduler).
#
#
# set name of job
#SBATCH --job-name=pf-vs
#
# set the number of processors/tasks needed
#SBATCH -n 64

#set an account to use
#if not used then default will be used
# for scavenger users, use this format:
#BATCH --account=patrick.kearns
# for contributing users, use this format:
##SBATCH --account=

# set max wallclock time  DD-HH:MM:SS

# the default time will be 1 hour if not set
#SBATCH --time=00-24:00:00

# set a memory request
#SBATCH --mem=128gb

# Set filenames for stdout and stderr.  %j can be used for the jobid.
# see "filename patterns" section of the sbatch man page for
# additional options
#SBATCH --error=%x-%j.err
#SBATCH --output=%x-%j.out
#

# set the partition where the job will run.  Multiple partitions can
# be specified as a comma separated list
# Use command "sinfo" to get the list of partitions
#SBATCH --partition=DGXA100
##SBATCH --partition=Intel6240,Intel6248,DGXA100

#When submitting to the GPU node, these following three lines are needed:
##SBATCH --gres=gpu:1
##SBATCH --export=NONE
#source /etc/profile
 

#Optional
# mail alert at start, end and/or failure of execution
# see the sbatch man page for other options
#SBATCH --mail-type=ALL
# send mail to this address
#SBATCH --mail-user=patrick.kearns@umb.edu

# Put your job commands here, including loading any needed
# modules or diagnostic echos.
## Assess quality with checkM
#cd /hpcstor6/scratch01/p/patrick.kearns/genomes/Assemblies_fastq_checkM/fastq/PFCALI/aggregated_genomes
#conda activate vs2
#source activate vs2
#configure vs database if necesary, only needs to be done once if you don't delete the databse like I tend to do :)
#virsorter setup -d ./virsorter2

#make new directory to store output
#mkdir virsorter_out

#run a for loop to do the thing
#for i in *.fasta;
#do
#virsorter run -i $i -w virsorter_out/$i --high-confidence-only 
#done

cd /hpcstor6/scratch01/p/patrick.kearns/genomes/Assemblies_fastq_checkM/fastq/PFCALI/aggregated_genomes

#rename files in each
mkdir -p viral_fragments

#rename fasta files for blast
for dir in */; do
    # Strip the trailing slash from the directory name for cleaner file naming
    dirname="${dir%/}"

    # Define the exact path to the fasta file
    target_file="${dirname}/final-viral-combined.fa"

    # Check if the file actually exists to avoid errors
    if [[ -f "$target_file" ]]; then
        # Move and rename the file into the new folder
        mv "$target_file" "viral_fragments/${dirname}_final-viral-combined.fa"
        
        # Optional: Print out what is happening so you can track progress
        echo "Moved $dirname to viral_fragments/${dirname}_final-viral-combined.fa"
    fi
done

#move the boundries file for analysis
mkdir -p viral_boundries

# 2. Loop through all directories in your current location
for dir in */; do
    # Strip the trailing slash from the directory name
    dirname="${dir%/}"

    # Define the exact path to the TSV file
    target_file="${dirname}/final-viral-boundary.tsv"

    # Check if the file actually exists
    if [[ -f "$target_file" ]]; then
        # Move and rename the file into the new folder
        mv "$target_file" "viral_boundries/${dirname}_final-viral-boundary.tsv"
        
        # Track progress in the terminal
        echo "Moved $dirname to viral_boundries/${dirname}_final-viral-boundary.tsv"
    fi
done

#get NCBI virus database for annotation
$HOME/ncbi-blast-2.14.1+/bin/update_blastdb.pl --decompress ref_viruses_rep_genome

#blastn the fragments
for i in *.fa; do
    $HOME/ncbi-blast-2.14.1+/bin/blastn \
    -db /hpcstor6/scratch01/p/patrick.kearns/genomes/Assemblies_fastq_checkM/fastq/PFCALI/ref_viruses_rep_genomes \
    -query "$i" \
    -out "${i}_blast_out.txt" \
    -max_target_seqs 1 \
    -outfmt 6
done

#bind all files togehter, adding a column for the name of the file
for f in *.tsv; do awk -v fname="$f" 'BEGIN{FS=OFS="\t"} {print fname, $0}' "$f"; done > combined_boundry_data.txt
for f in *.txt; do awk -v fname="$f" 'BEGIN{FS=OFS="\t"} {print fname, $0}' "$f"; done > combined_blast_data.txt
