#!/bin/bash
#set -e
set -u
set -o pipefail

# Creating directories
WD=$(date +%F); WD=2016-05-24; mkdir -p ./$WD; cd $WD

# Downloading data
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/overview.txt
curl -O ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.1.rna.fna.gz
gunzip -c plasmid.1.rna.fna.gz > plasmid.1.rna.fna
DB=plasmid.1.rna.fna
ls -lh $DB
grep '^>' $DB | wc -l

# Building a BLAST database
makeblastdb -in $DB -dbtype nucl -hash_index -parse_seqids

# Extracting query from BLAST databases with blastdbcmd
NAME=Ruegeria.sp
NAME=Deinococcus.geothermalis
NAME=Bacillus.megaterium
QUERY=$DB.$NAME.fasta
blastdbcmd -db $DB -entry all -outfmt "%i %t" | grep "$NAME" | \
awk '{print $1}' | blastdbcmd -db $DB -entry_batch - > $QUERY

# Running BLAST
PROGRAM=blastn
EVALUE=1e-05
OUTPUT=${PROGRAM}-$(basename $QUERY)-$(basename $DB).out
$PROGRAM -db $DB -query $QUERY -evalue $EVALUE -max_hsps 1 -outfmt 7 > $OUTPUT # -max_target_seqs 1

# Parsing BLAST outputs
# Extracting sbjct from BLAST databases with blastdbcmd
grep -v '#' $OUTPUT | awk '{print $2}' | uniq | blastdbcmd -db $DB -entry_batch - > $OUTPUT.fasta
grep '^>' $OUTPUT.fasta > $OUTPUT.fasta.header.txt

# Running R scripts
ln -s $OUTPUT blast.out
ln -s $OUTPUT.fasta.header.txt blast.out.fasta.header.txt
Rscript --vanilla ../scripts/2016-05-24.R 

# Print operating system characteristics
uname -a

echo "[$(date)] $0 has been successfully completed."
