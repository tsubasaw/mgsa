cat("\n  This R script create tests, small datasets for gene content analysis.\n\n")

# gene content (1,0) for strain 1, 2, and 3
s1 <- c(1,1,1,1,1,1,1,0,0,0)
s2 <- c(1,1,1,1,1,1,0,1,0,0)
s3 <- c(1,1,1,1,0,0,0,0,1,1)
mx <- cbind(s1,s2,s3) # gene content matrix
rownames(mx) <- paste0("p",1:nrow(mx)) # protein IDs
write.csv(mx, file="table.csv")

cat("# Pan-genome size\n")
nrow(mx)

cat("# No. of protein families for each genome\n")
apply(mx,2,sum)

cat("# Core genome size\n")
sum(apply(mx,1,sum) == ncol(mx))

cat("# Singleton (strain-specific genes)\n")
apply(mx,1,sum) == 1

# vennDiagram
# https://cell-innovation.nig.ac.jp/wiki/tiki-index.php?page=limma
#source("http://bioconductor.org/biocLite.R"); biocLite("limma")
library(limma)
vennDiagram(mx)

# clustering
heatmap(t(mx), scale="none", col=c("white","black"))
plot(hclust(dist(t(mx))))

## Distance
d <- dist(t(mx), method="binary") # Jaccard Distance

## Hierarchical Clustering
hc <- hclust(d, "average") # UPGMA
plot(hc)

## Classical multidimensional scaling (a.k.a. principal coordinates analysis)
loc <- cmdscale(d)
plot(loc[,1], loc[,2], type="n")
text(loc[,1], loc[,2], labels=rownames(loc))

# Pan-genome plot
#install.packages("vegan")
library(vegan)
sac <- specaccum(t(mx)); plot(sac, xlab="No. of genomes", ylab="No. of protein families")

# Print R version and packages # sessionInfo()
