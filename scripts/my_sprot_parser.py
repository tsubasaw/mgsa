
from Bio import SwissProt
for record in SwissProt.parse(open("uniprot_sprot.dat")):
    print (record.entry_name),"\t",
    print (record.description),"\t",
    print (record.organism),"\t",
    #print (record.organism_classification)
    for i in range(len(record.organism_classification)):
        print (record.organism_classification[i]),";",
    print

# wget ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/uniprot_sprot.dat.gz

# python my_sprot_parser.py > parse_sprot.txt

# Parsing Swissprot With Biopython https://www.biostars.org/p/59648/
# Biopython Tutorial and Cookbook | 10.1  Parsing Swiss-Prot files http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc136


