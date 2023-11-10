## Code for processing, annotating, etc.
### Assembly
```
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_628C.bam-1.fq -2 2023-08-15_256samples_628C.bam-2.fq -o 628C
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_629A2A.bam-1.fq -2 2023-08-15_256samples_629A2A.bam-2.fq -o 629A2A
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_630A.bam-1.fq -2 2023-08-15_256samples_630A.bam-2.fq -o 630A
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_631B.bam-1.fq -2 2023-08-15_256samples_631B.bam-2.fq -o 631B
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_632B.bam-1.fq -2 2023-08-15_256samples_632B.bam-2.fq -o 632B
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_633A.bam-1.fq -2 2023-08-15_256samples_633A.bam-2.fq -o 633A
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_638A1.bam-1.fq -2 2023-08-15_256samples_638A1.bam-2.fq -o 638A1
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_642A.bam-1.fq -2 2023-08-15_256samples_642A.bam-2.fq -o 642A
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_643E.bam-1.fq -2 2023-08-15_256samples_643E.bam-2.fq -o 643E
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_647D.bam-1.fq -2 2023-08-15_256samples_647D.bam-2.fq -o 647D
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_649D1.bam-1.fq -2 2023-08-15_256samples_649D1.bam-2.fq -o 649D1
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_CP10A.bam-1.fq -2 2023-08-15_256samples_CP10A.bam-2.fq -o CP10A
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_CP12E.bam-1.fq -2 2023-08-15_256samples_CP12E.bam-2.fq -o CP12E
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_CP13A.bam-1.fq -2 2023-08-15_256samples_CP13A.bam-2.fq -o CP13A
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_CP15C.bam-1.fq -2 2023-08-15_256samples_CP15C.bam-2.fq -o CP15C
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_CP17D.bam-1.fq -2 2023-08-15_256samples_CP17D.bam-2.fq -o CP17D
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_CP19F2.bam-1.fq -2 2023-08-15_256samples_CP19F2.bam-2.fq -o CP19F2
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_CP20G1.bam-1.fq -2 2023-08-15_256samples_CP20G1.bam-2.fq -o CP20G1
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_CP2A.bam-1.fq -2 2023-08-15_256samples_CP2A.bam-2.fq -o CP2A
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_CP3F2.bam-1.fq -2 2023-08-15_256samples_CP3F2.bam-2.fq -o CP3F2
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_CP5A.bam-1.fq -2 2023-08-15_256samples_CP5A.bam-2.fq -o CP5A
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_CP6A.bam-1.fq -2 2023-08-15_256samples_CP6A.bam-2.fq -o CP6A
$HOME/SPAdes-3.15.5-Linux/bin/spades.py -t 48 --isolate -1 2023-08-15_256samples_CP8C.bam-1.fq -2 2023-08-15_256samples_CP8C.bam-2.fq -o CP8C

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
source activate panaroo
cd /hpcstor6/scratch01/p/patrick.kearns/genomes/Assemblies_fastq_checkM/fastq/PF_Cali/annotations/gff_files
panaroo -i *.gff -o results --clean-mode strict --remove-invalid-genes -t 48
```
