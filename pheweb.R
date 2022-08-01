library(dplyr)
library(tidyr)
d <- read.table('/mnt/tank/scratch/dusolcev/pheweb/ESSE.HRC.vep_annotated_37_CSQ_severe_persed.txt', sep='\t', header = T)
d <- d[!duplicated(d$Uploaded_variation),]


dd <- read.table('/mnt/tank/scratch/dusolcev/pheweb/sites1.tsv', sep='\t', header = T)

ddd <- plyr::join(dd,d,by='Uploaded_variation')
ddd <- ddd %>% separate('Uploaded_variation',c('chrom','pos','ref','alt'), sep=':')

colnames(ddd) <- c('chrom','pos','ref','alt','rsids','nearest_genes','consequence')
write.table(ddd, "/mnt/tank/scratch/dusolcev/pheweb/sites_annotated.tsv", sep='\t',quote = F,row.names = F)

