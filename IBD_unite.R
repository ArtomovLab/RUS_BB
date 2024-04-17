library(dplyr)
library(tidyr)

'%!in%' <- function(x,y)!('%in%'(x,y))

G_annotations <- 'annotations of 1000G populations txt file'
outliers <- 'IDs of esse outliers txt file'
clusterN <- 'IDs of esse samples from cluster N txt file'
IBD <- 'esse + 1000 Genomes IBD txt file'
relatives <- 'IDs of esse relatives'

annotate <-  function(IBD_people) {

	sa_1000G = read.table(G_annotations,sep = '\t', header = T)
	sa_1000G <- sa_1000G[,c(1:3)]
	colnames(sa_1000G)[1] <- 'NP'
	samples <- read.table(outliers,sep = '\t', header = F)
	colnames(samples)[1] <- 'NP'
	IBD_people <- as.data.frame(merge(IBD_people,sa_1000G,by='NP',all.x=T))

	cl1 <- read.table(cluster1,sep = '\t', header = T)
	cl2 <- read.table(cluster2,sep = '\t', header = T)
	cl3 <- read.table(cluster3,sep = '\t', header = T)
	cl4 <- read.table(cluster4,sep = '\t', header = T)
	cl5 <- read.table(cluster5,sep = '\t', header = T)
	cl6 <- read.table(cluster6,sep = '\t', header = T)

	IBD_people <- IBD_people %>% mutate(Population = case_when(NP %in% cl1$x ~ 'cl1',
                                                           NP %in% cl2$x ~ 'cl2',
                                                           NP %in% cl3$x ~ 'cl3',
                                                           NP %in% cl4$x ~ 'cl4',
                                                           NP %in% cl5$x ~ 'cl5',
                                                           NP %in% cl6$x ~ 'cl6',
                                                           TRUE ~ Population))

	colnames(IBD_people) <- c('NP1','NP','IBD','POP1','Super_POP1')
	IBD_people <- as.data.frame(merge(IBD_people,sa_1000G,by='NP',all.x=T))

	IBD_people <- IBD_people %>% mutate(Population = case_when(NP %in% cl1$x ~ 'cl1',
                                                           NP %in% cl2$x ~ 'cl2',
                                                           NP %in% cl3$x ~ 'cl3',
                                                           NP %in% cl4$x ~ 'cl4',
                                                           NP %in% cl5$x ~ 'cl5',
                                                           NP %in% cl6$x ~ 'cl6',
                                                           TRUE ~ Population))

	colnames(IBD_people) <- c('NP2','NP1','IBD','POP1','Super_POP1','POP2','Super_POP2')

	IBD_people <- IBD_people %>% filter(!is.na(POP1) & !is.na(POP2))
	return (IBD_people)
}

for (i in c(1:22)) {
  print(i)
  if (i == 1) {
    IBD_all <- read.table(paste(IBD,i,sep=''),sep = '\t', header = F)
    #IBD_people <- IBD %>% group_by(V1,V3) %>% dplyr::summarise(N=sum(V8)) 
  }
  else {
    IBD <- read.table(paste(IBD,i,sep=''),sep = '\t', header = F)
    IBD_all <- rbind(IBD_all,IBD)
    #IBD_people <- IBD_people %>% group_by(V1,V3) %>% dplyr::summarise(N=sum(N)) 
    }
}

relatives <- read.table(relatives,sep = '\t', header = F)
IBD_all <- IBD_all %>% filter(V1 %!in% relatives$V1 & V3 %!in% relatives$V1)

write.table(IBD_all,paste(IBD,'all',sep='_'),sep = '\t',quote = F,row.names = F,col.names = T)

IBD_all <- IBD_all %>% filter(V9>1)
write.table(IBD_all,paste(IBD,'all','1',sep='_'),sep = '\t',quote = F,row.names = F,col.names = F)
IBD_all <- IBD_all %>% filter(V9>2)
write.table(IBD_all,paste(IBD,'all','2',sep='_'),sep = '\t',quote = F,row.names = F,col.names = F)
IBD_all <- IBD_all %>% filter(V9>3)
write.table(IBD_all,paste(IBD,'all','3',sep='_'),sep = '\t',quote = F,row.names = F,col.names = F)


IBD_all <- read.table(paste(IBD,'all',sep='_'),sep = '\t', header = F)

IBD_people <- IBD_all %>% group_by(V1,V3) %>% dplyr::summarise(N=sum(V9)) 
colnames(IBD_people)[1] <- 'NP'
IBD_people <- annotate(IBD_people)
write.table(IBD_people,paste(IBD,'all','len',sep='_'),sep = '\t',quote = F,row.names = F,col.names = F)

IBD_people <- IBD_all %>% group_by(V1,V3) %>% dplyr::summarise(N=sum(V8))
colnames(IBD_people)[1] <- 'NP'
IBD_people <- annotate(IBD_people)
write.table(IBD_people,paste(IBD,'all','LOD',sep='_'),sep = '\t',quote = F,row.names = F,col.names = F)

IBD_people <- IBD_all %>% group_by(V1,V3) %>% dplyr::summarise(N=n())
colnames(IBD_people)[1] <- 'NP'
IBD_people <- annotate(IBD_people)
write.table(IBD_people,paste(IBD,'all','count',sep='_'),sep = '\t',quote = F,row.names = F,col.names = F)

