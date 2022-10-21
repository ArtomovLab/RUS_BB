library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
#initial PCA
pca <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.hrc.dr08_pruned_PLINK2.eigenvec",sep = '\t', header = F)
pca <- pca[,c(2,3:12)]
colnames(pca) <- c('NP','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')
pca <-  pca %>%  mutate (REGION = case_when(substr(NP,1,2) == '36' ~ 'SAMARA',
                                            substr(NP,1,2) == '40' ~ 'SPB',
                                            substr(NP,1,2) == '53' ~ 'ORENBURG',
                                            TRUE ~ 'SPB'))

cols <- c("#999999", "#E69F00", "#56B4E9","#009E73","#0072B2","#CC79A7","#009E73","#F0E442")
g1 <- pca  %>% filter(!is.na(PC1)) %>%
  ggplot(aes(PC1,-PC2,col=REGION))+
  geom_point(alpha=0.4)+
  theme_classic()+
  scale_color_manual(values = cols)+
  theme(legend.position="none")+
  scale_fill_manual(values = cols)

g2 <- pca  %>% filter(!is.na(PC1)) %>%
  ggplot(aes(PC3,PC4,col=REGION))+
  geom_point(alpha=0.4)+
  theme_classic()+
  scale_color_manual(values = cols)+
  theme(legend.position="none")+
  scale_fill_manual(values = cols)

g3 <- pca  %>% filter(!is.na(PC1)) %>%
  ggplot(aes(PC5,PC6,col=REGION))+
  geom_point(alpha=0.4)+
  theme_classic()+
  scale_color_manual(values = cols)+
  theme(legend.position="none")+
  scale_fill_manual(values = cols)


#filtered PCA
#pca <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/pca/esse.hrc.dr08_pruned_PLINK2_2.eigenvec",sep = '\t', header = F)

step <- c()
dist <- c()
for (i in c(2:60)) {
  pca <- read.table(paste("/humgen/atgu1/methods/dusoltsev/biobank/HRC/pca/esse.hrc.dr08_pruned_PLINK2_",i,".eigenvec",sep=''),sep = '\t', header = F)
  pca <- pca[,c(2,3:12)]
  
  
  colnames(pca) <- c('NP','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')
  pca <-  pca %>%  mutate (REGION = case_when(substr(NP,1,2) == '36' ~ 'SAMARA',
                                              substr(NP,1,2) == '40' ~ 'SPB',
                                              substr(NP,1,2) == '53' ~ 'ORENBURG',
                                              TRUE ~ 'SPB'))
  step <- c(step,i)
  dist <- c(dist,max(dist(pca[,2:5])) )
  print(c(i,max(dist(pca[,2:5]))  ) )
}
res <- data.frame(step,dist)
res %>%
  ggplot(aes(step,dist))+
  geom_line()+
  theme_classic()+
  geom_vline(xintercept=47,col='red',linetype='dashed')

pca <- read.table(paste("/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.hrc.dr08_pruned_PLINK2_",'final',".eigenvec",sep=''),sep = '\t', header = F)
pca <- pca[,c(2,3:12)]

colnames(pca) <- c('NP','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10')
pca <-  pca %>%  mutate (REGION = case_when(substr(NP,1,2) == '36' ~ 'SAMARA',
                                            substr(NP,1,2) == '40' ~ 'SPB',
                                            substr(NP,1,2) == '53' ~ 'ORENBURG',
                                            TRUE ~ 'SPB'))
table(pca$REGION)
cols <- c("#999999", "#E69F00", "#56B4E9","#009E73","#0072B2","#CC79A7","#009E73","#F0E442")
g4 <- pca  %>% filter(!is.na(PC1)) %>%
  ggplot(aes(PC1,PC2,col=REGION))+ 
  geom_point(alpha=0.4)+
  theme_classic()+
  scale_color_manual(values = cols)+
  theme(legend.position="none")+
  scale_fill_manual(values = cols)

g5 <- pca %>% filter(!is.na(PC1)) %>%
  ggplot(aes(PC3,PC4,col=REGION))+
  geom_point(alpha=0.4)+
  theme_classic()+
  scale_color_manual(values = cols)+
  theme(legend.position="none")+
  scale_fill_manual(values = cols)

g6 <- pca  %>% filter(!is.na(PC1)) %>%  mutate(s=case_when(NP == '40-K-1061' ~ '1',
                                                              TRUE ~ '0')) %>%
  ggplot(aes(PC5,PC6,col=s))+
  geom_point(alpha=0.4)+
  theme_classic()+
  scale_color_manual(values = cols)+
  theme(legend.position="none")+
  scale_fill_manual(values = cols)
grid.arrange(g1,g2,g4,g5, ncol=2, nrow=2)


write.table(pca[,c(2:11)],paste("/humgen/atgu1/methods/dusoltsev/biobank/HRC/UMAP/esse.hrc.dr08_pruned_PLINK2_final.eigenvec",sep=''), sep='\t',quote = F,row.names = F)
##umap
umap = read.table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/UMAP/esse.hrc.dr08_pruned_PLINK2_final.eigenvec_UMAP_PC10_NC3_NN10_MD0.001_euclidean_20221001_091211.txt',sep = ' ', header = F)
umap <- cbind(umap,pca)
umap %>%
  ggplot(aes(V1,V2,col=REGION))+
  geom_point()+
  theme_classic()

#PCA with 1000G

sa_1000G = read.table('/humgen/atgu1/methods/dusoltsev/biobank/1kg_annotations.txt',sep = '\t', header = T)
colnames(sa_1000G)[1] <- 'NP'

pca <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS_pca.txt",sep = '\t', header = T)
pca$scores <- gsub('\\{scores\\:\\[','',pca$scores)
pca$scores <- gsub('\\]\\}','',pca$scores)
pca$scores <- gsub('\\[','',pca$scores)
pca$scores <- gsub('\\]','',pca$scores)
pca <-  pca %>%separate(scores, c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10'),sep=',')
for (i in 2:11) {
  pca[,i] <- as.numeric(pca[,i])
}
colnames(pca)[1] <- 'NP'

pca <- as.data.frame(merge(pca,sa_1000G,by='NP',all.x=T,all.y=T))


pca <- pca %>%  mutate (Super_REGION = case_when(substr(NP,1,2) == '36' ~ 'RUS',
                                                 substr(NP,1,2) == '40' | substr(NP,1,3) == 'KBL' ~ 'RUS',
                                                 substr(NP,1,2) == '53' ~ 'RUS',
                                                 TRUE ~ SuperPopulation)) 
pca <- pca %>%  mutate (REGION = case_when(substr(NP,1,2) == '36' ~ 'SAMARA',
                                           substr(NP,1,2) == '40' | substr(NP,1,3) == 'KBL' ~ 'SPB',
                                           substr(NP,1,2) == '53' ~ 'ORENBURG',
                                           TRUE ~ Population)) 

pca <- pca %>% filter(!is.na(Super_REGION)) %>% filter(!is.na(PC1)) 

cols <- c("#009E73", "#F0E442", "#0072B2", "#D55E00",'#000000', "#CC79A7")

table(pca$Super_REGION)
g7 <- pca  %>% filter(!is.na(PC1)) %>%

  ggplot(aes(PC1,PC2,col=Super_REGION))+
  geom_point(alpha=0.4)+
  theme_classic()+
  scale_color_manual(values = cols)+
  theme(legend.position="none")+
  scale_fill_manual(values = cols)

#PCA with 1000G RUR only

sa_1000G = read.table('/humgen/atgu1/methods/dusoltsev/biobank/1kg_annotations.txt',sep = '\t', header = T)
colnames(sa_1000G)[1] <- 'NP'

pca <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS_EUR_pca.txt",sep = '\t', header = T)
pca$scores <- gsub('\\{scores\\:\\[','',pca$scores)
pca$scores <- gsub('\\]\\}','',pca$scores)
pca$scores <- gsub('\\[','',pca$scores)
pca$scores <- gsub('\\]','',pca$scores)
pca <-  pca %>%separate(scores, c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10'),sep=',')
for (i in 2:11) {
  pca[,i] <- as.numeric(pca[,i])
}
colnames(pca)[1] <- 'NP'

pca <- as.data.frame(merge(pca,sa_1000G,by='NP',all.x=T,all.y=T))


pca <- pca %>%  mutate (Super_REGION = case_when(substr(NP,1,2) == '36' ~ 'RUS',
                                                 substr(NP,1,2) == '40' | substr(NP,1,3) == 'KBL' ~ 'RUS',
                                                 substr(NP,1,2) == '53' ~ 'RUS',
                                                 TRUE ~ SuperPopulation)) 
pca <- pca %>%  mutate (REGION = case_when(substr(NP,1,2) == '36' ~ 'SAMARA',
                                           substr(NP,1,2) == '40' | substr(NP,1,3) == 'KBL' ~ 'SPB',
                                           substr(NP,1,2) == '53' ~ 'ORENBURG',
                                           TRUE ~ Population)) 

pca <- pca %>% filter(!is.na(Super_REGION)) %>% filter(!is.na(PC1)) 

cols <- c("#56B4E9", "#E69F00",'#000000', "#D55E00")
table(pca$REGION)
g8 <- pca  %>% filter(!is.na(PC1)) %>%
  
  mutate(REGION = case_when(Super_REGION %in% c('RUS') ~ 'RUS',
                            REGION %in% c('TSI','IBS') ~ 'TSI_IBS',
                            REGION %in% c('CEU','GBR') ~ 'CEU_GBR',
                            TRUE ~ REGION)) %>%
  ggplot(aes(-PC1,PC2,col=REGION))+
  geom_point(alpha=0.4)+
  theme_classic()+
  scale_color_manual(values = cols)+
  theme(legend.position="none")+
  scale_fill_manual(values = cols)

grid.arrange(g7,g8, ncol=2, nrow=1)


#without outliers
#PCA with 1000G

sa_1000G = read.table('/humgen/atgu1/methods/dusoltsev/biobank/1kg_annotations.txt',sep = '\t', header = T)
colnames(sa_1000G)[1] <- 'NP'

pca <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS_15_pca_exclude.txt",sep = '\t', header = T)
pca$scores <- gsub('\\{scores\\:\\[','',pca$scores)
pca$scores <- gsub('\\]\\}','',pca$scores)
pca$scores <- gsub('\\[','',pca$scores)
pca$scores <- gsub('\\]','',pca$scores)
pca <-  pca %>%separate(scores, c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','PC11','PC12','PC13','PC14','PC15'),sep=',')
for (i in 2:16) {
  pca[,i] <- as.numeric(pca[,i])
}
colnames(pca)[1] <- 'NP'

pca <- as.data.frame(merge(pca,sa_1000G,by='NP',all.x=T,all.y=T))


pca <- pca %>%  mutate (Super_REGION = case_when(substr(NP,1,2) == '36' ~ 'RUS',
                                                 substr(NP,1,2) == '40' | substr(NP,1,3) == 'KBL' ~ 'RUS',
                                                 substr(NP,1,2) == '53' ~ 'RUS',
                                                 TRUE ~ SuperPopulation)) 
pca <- pca %>%  mutate (REGION = case_when(substr(NP,1,2) == '36' ~ 'SAMARA',
                                           substr(NP,1,2) == '40' | substr(NP,1,3) == 'KBL' ~ 'SPB',
                                           substr(NP,1,2) == '53' ~ 'ORENBURG',
                                           TRUE ~ Population)) 

pca <- pca %>% filter(!is.na(Super_REGION)) %>% filter(!is.na(PC1)) 

cols <- c("#009E73", "#F0E442", "#0072B2", "#D55E00",'#000000', "#CC79A7")

for (i in pca$REGION) {
  d <- pca %>% filter(REGION == i)
  write.table(d$NP,paste("/humgen/atgu1/methods/dusoltsev/biobank/HRC/fst/samples_WGS_",i,".txt",sep=''), sep='\t',quote = F,row.names = F)
}

table(pca$Super_REGION)
g7 <- pca  %>% filter(!is.na(PC1)) %>%
  
  ggplot(aes(PC1,PC2,col=Super_REGION))+
  geom_point(alpha=0.4)+
  theme_classic()+
  scale_color_manual(values = cols)+
  theme(legend.position="none")+
  scale_fill_manual(values = cols)


write.table(pca[,c(2:16)],paste("/humgen/atgu1/methods/dusoltsev/biobank/HRC/UMAP/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS_15_pca_exclude.txt",sep=''), sep='\t',quote = F,row.names = F)
pca <- pca %>% filter(Super_REGION != 'RUS')
write.table(pca[,c(2:16)],paste("/humgen/atgu1/methods/dusoltsev/biobank/HRC/UMAP/1000G_WGS_15_pca_exclude.txt",sep=''), sep='\t',quote = F,row.names = F)
#UMAP
#esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS_15_pca_exclude_UMAP_PC15_NC3_NN10_MD0.001_euclidean_20221003_052303.txt
#1000G_WGS_15_pca_exclude.txt
umap = read.table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/UMAP/1000G_WGS_15_pca_exclude_UMAP_PC15_NC3_NN10_MD0.001_euclidean_20221003_052817.txt',sep = ' ', header = F)
umap <- cbind(umap,pca)

cl1 <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/fst/esse.hrc.dr08_pruned_PLINK_CLUSTERS_61new.txt",sep = '\t', header = T)
cl2 <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/fst/esse.hrc.dr08_pruned_PLINK_CLUSTERS_62new.txt",sep = '\t', header = T)
cl3 <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/fst/esse.hrc.dr08_pruned_PLINK_CLUSTERS_63new.txt",sep = '\t', header = T)
cl4 <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/fst/esse.hrc.dr08_pruned_PLINK_CLUSTERS_64new.txt",sep = '\t', header = T)
cl5 <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/fst/esse.hrc.dr08_pruned_PLINK_CLUSTERS_65new.txt",sep = '\t', header = T)
cl6 <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/fst/esse.hrc.dr08_pruned_PLINK_CLUSTERS_66new.txt",sep = '\t', header = T)

umap <- umap %>% mutate(POP = case_when(NP %in% cl1$x ~ 'cl1',
                                      NP %in% cl2$x ~ 'cl2',
                                      NP %in% cl3$x ~ 'cl3',
                                      NP %in% cl4$x ~ 'cl4',
                                      NP %in% cl5$x ~ 'cl5',
                                      NP %in% cl6$x ~ 'cl6'
                                      ))

umap %>% 
  #filter(!is.na(POP)) %>%
  ggplot(aes(V1,V2,col=SuperPopulation))+
  geom_point()+
  theme_classic()+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)



#PCA with 1000G RUR only

sa_1000G = read.table('/humgen/atgu1/methods/dusoltsev/biobank/1kg_annotations.txt',sep = '\t', header = T)
colnames(sa_1000G)[1] <- 'NP'

pca <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS_EUR_pca_exclude.txt",sep = '\t', header = T)
pca$scores <- gsub('\\{scores\\:\\[','',pca$scores)
pca$scores <- gsub('\\]\\}','',pca$scores)
pca$scores <- gsub('\\[','',pca$scores)
pca$scores <- gsub('\\]','',pca$scores)
pca <-  pca %>%separate(scores, c('PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10'),sep=',')
for (i in 2:11) {
  pca[,i] <- as.numeric(pca[,i])
}
colnames(pca)[1] <- 'NP'

pca <- as.data.frame(merge(pca,sa_1000G,by='NP',all.x=T,all.y=T))


pca <- pca %>%  mutate (Super_REGION = case_when(substr(NP,1,2) == '36' ~ 'RUS',
                                                 substr(NP,1,2) == '40' | substr(NP,1,3) == 'KBL' ~ 'RUS',
                                                 substr(NP,1,2) == '53' ~ 'RUS',
                                                 TRUE ~ SuperPopulation)) 
pca <- pca %>%  mutate (REGION = case_when(substr(NP,1,2) == '36' ~ 'SAMARA',
                                           substr(NP,1,2) == '40' | substr(NP,1,3) == 'KBL' ~ 'SPB',
                                           substr(NP,1,2) == '53' ~ 'ORENBURG',
                                           TRUE ~ Population)) 

pca <- pca %>% filter(!is.na(Super_REGION)) %>% filter(!is.na(PC1)) 

cols <- c("#56B4E9", "#E69F00",'#000000', "#D55E00")
table(pca$REGION)
g8 <- pca  %>% filter(!is.na(PC1)) %>%
  
  mutate(REGION = case_when(Super_REGION %in% c('RUS') ~ 'RUS',
                            REGION %in% c('TSI','IBS') ~ 'TSI_IBS',
                            REGION %in% c('CEU','GBR') ~ 'CEU_GBR',
                            TRUE ~ REGION)) %>%
  ggplot(aes(PC1,PC2,col=REGION))+
  geom_point(alpha=0.4)+
  theme_classic()+
  scale_color_manual(values = cols)+
  theme(legend.position="none")+
  scale_fill_manual(values = cols)

grid.arrange(g4,g7,g8, ncol=3, nrow=1)

write.table(pca[,c(2:11)],paste("/humgen/atgu1/methods/dusoltsev/biobank/HRC/UMAP/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS_EUR_pca_exclude.txt",sep=''), sep='\t',quote = F,row.names = F)
#UMAP
umap = read.table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/UMAP/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS_EUR_pca_exclude_UMAP_PC10_NC3_NN10_MD0.001_euclidean_20221001_085953.txt',sep = ' ', header = F)
umap <- cbind(umap,pca)
umap %>%
  ggplot(aes(V1,V2,col=REGION))+
  geom_point()+
  theme_classic()