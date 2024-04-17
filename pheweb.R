library(dplyr)
library(tidyr)

outputFile <- "parsed output table from CSQ parser"

d <- read.table(outputFile, sep='\t', header = T)
d <- d[!duplicated(d$Uploaded_variation),]

dd <- read.table('sites.tsv', sep='\t', header = T)

ddd <- plyr::join(dd,d,by='Uploaded_variation')
ddd <- ddd %>% separate('Uploaded_variation',c('chrom','pos','ref','alt'), sep=':')

colnames(ddd) <- c('chrom','pos','ref','alt','rsids','nearest_genes','consequence')
write.table(ddd, "sites_annotated.tsv", sep='\t',quote = F,row.names = F)

