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
#SBATCH --job-name=pf_blast
#
# set the number of processors/tasks needed
#SBATCH -n 24

#set an account to use
#if not used then default will be used
# for scavenger users, use this format:
#BATCH --account=patrick.kearns
# for contributing users, use this format:
##SBATCH --account=

# set max wallclock time  DD-HH:MM:SS

# the default time will be 1 hour if not set
#SBATCH --time=00-06:00:00

# set a memory request
#SBATCH --mem=48gb

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
cd /hpcstor6/scratch01/p/patrick.kearns/pseudomonas

#get refseq prok genomes database for commandline blast
#$HOME/ncbi-blast-2.14.1+/bin/update_blastdb.pl ref_prok_rep_genomes --decompress

#blast it
$HOME/ncbi-blast-2.14.1+/bin/blastn -query pan_genome_reference.fa -db ref_prok_rep_genomes -out blastn_results.tsv -outfmt "6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore" -evalue 1e-5 -perc_identity 90 -max_target_seqs 10 -num_threads 24

