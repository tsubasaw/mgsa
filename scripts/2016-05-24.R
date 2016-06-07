cat("\n  This R script adds taxonomic information to BLAST output.\n\n")

# Set Working Directory
#setwd("~/projects/mgsa/2016-05-24/")

# Loading Data into R
## ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/overview.txt
ov <- read.delim('overview.txt', stringsAsFactors=FALSE, na.strings="-", check.names=FALSE); dim(ov)

colnames(ov)[1] <- "Organism"
genus = apply(as.matrix(ov[,1]), MARGIN=c(1,2), function(x) unlist(strsplit(x," "))[1])[,1]
taxon <- unique(cbind(genus, apply(ov[,2:4], 1, paste, collapse=';'))) # Kingdom;Group;SubGroup
TF = regexpr("-", taxon[,1], perl=TRUE) > 0; sum(TF); taxon = taxon[!TF,]
TF = regexpr("^[A-Z]", taxon[,1], perl=TRUE) > 0; sum(TF); taxon = taxon[TF,]
TF = regexpr("^(Viruses)", taxon[,2], perl=TRUE) > 0; sum(TF); taxon = taxon[!TF,]

## BLAST output file
blast <- read.delim('blast.out', header = FALSE, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "#", col.names=c("query id", "subject id", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score"))

## FASTA headers
annot <- read.delim('blast.out.fasta.header.txt', header = FALSE, stringsAsFactors = FALSE, quote = "")
annot <- unique(annot[,1])

# Convert BLAST sbjct_id Into organism/taxon names
blast$definition <- apply(as.matrix(blast[,2]), MARGIN=c(1,2), function(x) annot[ grep(pattern=unlist(strsplit(x, split="\\|"))[4], x=annot) ] )

#x <- blast$definition[3]; grep(pattern=unlist(strsplit(x, split=" "))[2], x=taxon[,1]); #taxon[c(800,801), 1]
blast$taxonomy <- apply(as.matrix(blast$definition), MARGIN=c(1,2), function(x) taxon[ grep(pattern=unlist(strsplit(x, split=" "))[2], x=taxon[,1]), 2 ] )[,1]

# Exporting Data
write.table(blast, file="output.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Print R version and packages
sessionInfo()


