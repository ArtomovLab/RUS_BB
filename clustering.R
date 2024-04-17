library(SVDFunctions)
library(mclust)
library(tidyr)
library(dplyr)
library(ggplot2)

pca_without_outliers = 'PCA esse without outliers txt file'

pca <- read.table(pca_without_outliers,sep = '\t', header = F)
row.names(pca) <- pca$V2
pca <- pca[,3:12]
colnames(pca) <- c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')

plotPCA(pca)

set.seed(60)
caseCl <- estimateCaseClusters(PCA = pca[,1:4],
                               plotBIC = TRUE,
                               plotDendrogram = TRUE,
                               clusters = 6,
                               keepSamples = row.names(pca))
plotPCA(pca,caseCl)
table(as.character(caseCl$samples))

caseCl <- mergeCluster(caseCl, cluster = 5) 

pca <- read.table(pca_without_outliers,sep = '\t', header = F)
pca$CLUSTER <- caseCl$samples
pca <- pca %>% mutate(CLUSTER = case_when(CLUSTER == 5 ~ 1,
                                          CLUSTER == 1 ~ 2,
                                          CLUSTER == 4 ~ 4,
                                          CLUSTER == 2 ~ 5,
                                          CLUSTER == 3 ~ 3,
                                          CLUSTER == 6 ~ 6)) 


table(pca$CLUSTER)
pca %>%
  ggplot(aes(V3,V4,col=as.character(CLUSTER)))+
  geom_point()+
  theme_classic()



table(pca$CLUSTER)
write.table(pca,'clusters.txt',sep='\t', row.names=FALSE, quote=F)
for (i in 1:length(unique(pca$CLUSTER))) {
  p <- pca %>% filter(CLUSTER == i)
  write.table(p[,2], paste('cluster',i,'.txt',sep=''),sep='\t', row.names=FALSE, quote=F)
}

