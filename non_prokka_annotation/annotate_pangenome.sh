#!/usr/bin/env bash
# =============================================================================
# Functional annotation of a Panaroo pangenome (pan_genome_reference.fa)
# -----------------------------------------------------------------------------
# Sources (as run Oct 2026):
#   UniRef90  - UniProt REST, clusters with members in Pseudomonas (taxid 286,
#               1,303,969 clusters) + Caudoviricetes (taxid 2731619, 243,806)
#   COG 2020  - NCBI FTP pub/COG/COG2020/data (cog-20.fa.gz, cog-20.cog.csv,
#               cog-20.def.tab, fun-20.tab)
#   KOfam     - genome.jp ftp/db/kofam (profiles.tar.gz, ko_list.gz); prokaryote.hal
#   Pfam-A    - EBI Pfam current_release (38.2) + GO pfam2go mapping
# Tools: DIAMOND 2.2.8, HMMER 3.4, KofamScan 1.3.0 (bioconda)
#   conda create -n annot -c conda-forge -c bioconda diamond hmmer kofamscan pandas requests
# Thresholds:
#   UniRef90 : blastp --sensitive, e<=1e-5, query cover>=50%, top 5 kept
#   COG      : blastp --sensitive, e<=1e-5, query cover>=50%, identity>=25%, best hit;
#              COG(s) of the hit kept only where the COG footprint overlaps >=30% of the alignment
#   KOfam    : family-specific adaptive thresholds (rows marked '*'); best-scoring KO per gene
#   Pfam     : hmmscan --cut_ga (gathering thresholds)
# Translation: genetic code 11 (done in Python, see translate step); genes <100 aa flagged.
# Merging / UniProt cross-reference retrieval: merge_annotations.py
# =============================================================================
set -euo pipefail
FA=${1:-pan_genome_reference.fa}
CPU=${CPU:-10}
mkdir -p db res

# ---- databases (skip if present) -------------------------------------------
C=https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2020/data
for f in cog-20.fa.gz cog-20.cog.csv cog-20.def.tab fun-20.tab; do
  [ -s db/$f ] || curl -sSL --retry 5 -o db/$f $C/$f; done
[ -s db/ko_list ]  || { curl -sSL -o db/ko_list.gz https://www.genome.jp/ftp/db/kofam/ko_list.gz; gunzip -f db/ko_list.gz; }
[ -d db/profiles ] || { curl -sSL -C - --speed-limit 2000 --speed-time 90 --retry 20 \
                          -o db/profiles.tar.gz https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
                        tar -xzf db/profiles.tar.gz -C db; }
[ -s db/Pfam-A.hmm.h3m ] || { curl -sSL -o db/Pfam-A.hmm.gz \
     https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
     gunzip -f db/Pfam-A.hmm.gz; hmmpress -f db/Pfam-A.hmm; }
[ -s db/pfam2go ] || curl -sSL -o db/pfam2go https://current.geneontology.org/ontology/external2go/pfam2go
U="https://rest.uniprot.org/uniref/stream?format=fasta&compressed=true&query="
[ -s db/uniref90_sel.dmnd ] || {
  curl -sSL -o db/uniref90_pseudomonas.fa.gz    "${U}%28taxonomy_id%3A286%29%20AND%20%28identity%3A0.9%29"
  curl -sSL -o db/uniref90_caudoviricetes.fa.gz "${U}%28taxonomy_id%3A2731619%29%20AND%20%28identity%3A0.9%29"
  gzip -dc db/uniref90_pseudomonas.fa.gz db/uniref90_caudoviricetes.fa.gz > db/uniref90_sel.fa
  diamond makedb --in db/uniref90_sel.fa -d db/uniref90_sel -p $CPU && rm -f db/uniref90_sel.fa; }
[ -s db/cog20.dmnd ] || diamond makedb --in db/cog-20.fa.gz -d db/cog20 -p $CPU

# ---- translate (code 11) ------------------------------------------------------
python merge_annotations.py translate "$FA" pangenome_proteins.faa translation_qc.tsv

# ---- searches -------------------------------------------------------------------
diamond blastp -q pangenome_proteins.faa -d db/uniref90_sel -o res/uniref90.tsv --sensitive \
  -e 1e-5 --query-cover 50 -k 5 -p $CPU \
  -f 6 qseqid sseqid pident length qlen slen evalue bitscore qcovhsp scovhsp stitle
diamond blastp -q pangenome_proteins.faa -d db/cog20 -o res/cog.tsv --sensitive \
  -e 1e-5 --query-cover 50 --id 25 -k 1 -p $CPU \
  -f 6 qseqid sseqid pident length qlen slen evalue bitscore qcovhsp sstart send
hmmscan --cut_ga --cpu $CPU --noali --domtblout res/pfam.domtbl db/Pfam-A.hmm pangenome_proteins.faa > /dev/null
exec_annotation -o res/kofam_detail.tsv -p db/profiles/prokaryote.hal -k db/ko_list \
  --cpu $CPU --tmp-dir tmp_kofam -f detail-tsv pangenome_proteins.faa

# ---- merge (+ UniProt GO/KEGG/eggNOG/Pfam/InterPro xrefs, KO -> pathways) --------
python merge_annotations.py merge --rtab gene_presence_absence2.Rtab \
  --meta isolate_inhibi_category.txt --kegg full_kegg.txt --swissprot best_hits.tsv \
  --out pangenome_annotation.tsv
