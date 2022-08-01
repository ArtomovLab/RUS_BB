args = commandArgs(trailingOnly=TRUE)
i <- as.integer(args[1]) - 1
library(adamethods)
pca <- read.table(paste("/humgen/atgu1/methods/dusoltsev/biobank/HRC/pca/esse.hrc.dr08_pruned_PLINK2_",i,".eigenvec",sep=''),sep = '\t', header = F)
pca <- pca[,c(2,3:6)]
colnames(pca) <- c('NP','PC1','PC2','PC3','PC4')
row.names(pca) <- pca$NP
pca$NP <- NULL
pca <- as.matrix(pca,labels=TRUE)
write.table(row.names(pca[c(do_knno(pca, 5, 3)),1,drop=F]),'/humgen/atgu1/methods/dusoltsev/biobank/HRC/test.txt', sep='\t',quote = F,row.names = F,col.names = F)
