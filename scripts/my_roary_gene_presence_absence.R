cat("\n  This R script analyses Roary output (gene_presence_absence.csv).\n\n")

# Set Working Directory
#setwd("~/projects/mgsa/roary/2016-06-07") # i 95
#setwd("~/projects/mgsa/roary/2016-06-21") # i 50

# Loading Data into R
d.f <- read.csv('analysis/gene_presence_absence.csv', stringsAsFactors=FALSE, check.names=FALSE)
	
# Exporting Data
#for(rgxp in c('ribosomal.protein', 'elongation.factor')){ # rgxp <- 'ribosomal'
for(rgxp in c('Trb')){ # rgxp <- 'transfer'
 i <- grep(pattern=rgxp, x=d.f$Annotation, ignore.case=FALSE); #d.f[i, c(3,4)]
 write.csv(d.f[i, c(3,4,15:ncol(d.f))], file=paste0('gene_content.',rgxp,'.csv'))
 write.csv(as.matrix(sort(table(d.f$Annotation[i]), decreasing=TRUE)), file=paste0('gene_count.',rgxp,'.csv'))
}

# Print R version and packages
sessionInfo()
