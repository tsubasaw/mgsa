#!/bin/bash
set -e
set -u
set -o pipefail

wget ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/RELEASE.metalink \
     ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz \
     ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz

for FILE in uniprot_sprot.dat.gz uniprot_sprot.fasta.gz; do
    echo; echo $FILE
    md5sum $FILE
    grep -A 3 "file name=\"$FILE\"" RELEASE.metalink
    gunzip -c $FILE > $(basename $FILE .gz)
done

makeblastdb -in uniprot_sprot.fasta -dbtype prot -hash_index -parse_seqids

# Print operating system characteristics
uname -a

echo "[$(date)] $0 has been successfully completed."

: <<'#__COMMENT_OUT__'

#__COMMENT_OUT__


