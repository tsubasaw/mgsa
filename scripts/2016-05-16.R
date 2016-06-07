cat("\n  This R script creates protein phylogenetic profiles, and identifies taxon-specific genes.\n\n")

# Set Working Directory
#setwd("~/projects/mgsa/")

# List files in a directory
files <- list.files(path="data", pattern="tblastn-NC_008357.faa-", full.names=TRUE)
#files <- files[-grep(pattern="data/tblastn-NC_008357.faa-NC_008357.fna.out", x=files)]

# Loading Data into R
ld <- lapply(files, read.delim, header = FALSE, comment.char = "#", col.names=c("query id", "subject id", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score"))

# name each list item with the filename (sans full path)
names(ld) <- basename(files)

#head(ld[[1]])

# Using do.call() and rbind to merge this data together
d.f <- do.call(rbind, ld)

# extract ref from query/subject id
d.f$query.id <- sub("gi\\|(.+)\\|ref\\|(.+)\\|", "\\2", d.f$query.id)
d.f$subject.id <- sub("gi\\|(.+)\\|ref\\|(.+)\\|", "\\2", d.f$subject.id)

# reshape data: reshape2::dcast turns long data into wide data.
library(reshape2)
mx <- dcast(d.f, query.id ~ subject.id, value.var = "evalue")
rownames(mx) <- mx[,1]; mx <- mx[,-1]
mx[!is.na(mx)] <- 1
mx[ is.na(mx)] <- 0

# Exploring Data Visually
#pdf(file="Rplots.pdf")
#par(mfrow=c(2,2), cex = 0.7, oma = c(0, 1.5, 0, 0)) # oma = c(bottom, left, top, right)
par(mfrow=c(1,1), cex = 0.8, oma = c(2, 2, 0, 0)) # oma = c(bottom, left, top, right)
lab <- "Numbers of QUERY genes found in each genome"
barplot(apply(mx, 2, sum), las = 2, ylab = lab)
#barplot(apply(mx, 2, sum), horiz = TRUE, las = 1, xlab = lab)
#dev.off()

# clustering
heatmap(as.matrix(mx), scale="none", col=c("white","black"))
#plot(hclust(dist(mx)))

## Distance
d <- dist(mx, method="manhattan")

## Hierarchical Clustering
hc <- hclust(d)
plot(hc, hang = -1)
group <- cutree(hc, h = 0:1)
colnames(group) <- paste0("cutree.h", colnames(group))

## Classical multidimensional scaling (a.k.a. principal coordinates analysis)
loc <- cmdscale(d)
plot(loc[,1], loc[,2], type="n")
text(loc[,1], loc[,2], labels=rownames(loc))

# annotation
system("grep '^>' data/NC_008357.faa > data/NC_008357.faa.header.txt")
annot <- read.delim("data/NC_008357.faa.header.txt", header = FALSE, quote = "")
annot <- annot[,1]
annotation <- apply(as.matrix(rownames(mx)), MARGIN=c(1,2), function(x) annot[ grep(pattern=x, x=annot) ] )

# taxon-specific genes
num <- which(regexpr(pattern="NC_008357", text=colnames(mx)) > 0)
#rownames(mx)[apply(mx[,-num], 1, sum) == 0]
specific <- as.numeric(apply(mx[,-num], 1, sum) == 0)
cat("# taxon-specific genes\n")
annotation[apply(mx[,-num], 1, sum) == 0]

# Exporting Data
write.csv(mx, file="phylogenetic_profiles.csv")
write.csv(data.frame(mx, group, specific, annotation), file="phylogenetic_profiles_supp.csv")

# Print R version and packages
sessionInfo()


