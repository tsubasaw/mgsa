#!/bin/bash
#set -e
set -u
set -o pipefail

# Creating directories
#mkdir -p ./{data/$(date +%F),analysis/$(date +%F)}
mkdir -p ./$(date +%F)

# Downloading data
ACCESSIONs=(NC_001621 NC_005088 NC_007337 NC_008459 NC_001735 NC_005793 NC_008357 NC_010935) # IncP plasmids
ACCESSIONs=(NC_001621 NC_001735 NC_008357) # http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3525142/table/mss210-T1/
URL=ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Plasmids
for ACCESSION in ${ACCESSIONs[@]}; do wget -nv -P data/ $URL/faa/$ACCESSION.faa $URL/fna/$ACCESSION.fna; done

cd data/

# Inspecting and Manipulating data
grep '^>' NC_*.fna
grep -c '^>' NC_*.faa

# Building a BLAST database with local sequences
#for DB in NC_*.faa; do makeblastdb -in $DB -dbtype prot -hash_index -parse_seqids; done
for DB in NC_*.fna; do makeblastdb -in $DB -dbtype nucl -hash_index -parse_seqids; done

# Running BLAST
QUERY=NC_008357.faa
PROGRAM=tblastn
EVALUE=1e-05

for DB in NC_*.fna; do 
 OUTPUT=${PROGRAM}-$(basename $QUERY)-$(basename $DB).out;
 $PROGRAM -db $DB -query $QUERY -evalue $EVALUE -max_hsps 1 -outfmt 7 > $OUTPUT; # -max_target_seqs $(grep -c '^>' $DB)
done

# Inspecting and Manipulating BLAST $OUTPUT
grep -vc '^#' tblastn-*.out
#cat tblastn-*.out | grep -v '^#' | cut -f1 | cut -d'|' -f4

# Running R scripts
cd ..
Rscript --vanilla scripts/2016-05-16.R 

# Print operating system characteristics
uname -a

echo "[$(date)] $0 has been successfully completed."





