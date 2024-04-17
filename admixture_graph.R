library(dplyr)
library(tidyr)
library(ggplot2)
'%!in%' <- function(x,y)!('%in%'(x,y))

fam_file <- 'plink fam file of esse + 1000 Genomes'
pop_file <- 'txt file with populations'
clusterN <- 'txt file IDs of esse samples from cluster N txt file'
admixture_Q_file <- 'txt file of ADMIXTURE results'
esse_G1K_pruned_pca_exclude = 'PCA merged esse and 1000G pruned txt file, outliers excluded'

fam <- read.table(fam_file)

s <- read.table(pop_file)
fam <- cbind(fam,s)
fam <- fam[,c(2,7)]
colnames(fam) <- c('NP','POP')
table(fam$POP)


cl1 <- read.table(cluster1,sep = '\t', header = T)
cl2 <- read.table(cluster2,sep = '\t', header = T)
cl3 <- read.table(cluster3,sep = '\t', header = T)
cl4 <- read.table(cluster4,sep = '\t', header = T)
cl5 <- read.table(cluster5,sep = '\t', header = T)
cl6 <- read.table(cluster6,sep = '\t', header = T)


fam <- fam %>% mutate(POP = case_when(NP %in% cl1$x ~ 'cl1',
                                                           NP %in% cl2$x ~ 'cl2',
                                                           NP %in% cl3$x ~ 'cl3',
                                                           NP %in% cl4$x ~ 'cl4',
                                                           NP %in% cl5$x ~ 'cl5',
                                                           NP %in% cl6$x ~ 'cl6',
                                                           TRUE ~ POP))

table(fam$POP)

pop <- read.table(admixture_Q_file)
pop <- cbind(fam,pop)
colnames(pop) <- c('NP','POP','CEU_GBR','FIN','CHB_CHS_JPT','AMR','CDX_KHV','IBS_TSI','SAS','AFR')

pop <- pop %>% filter(POP %in% c('cl1','cl2','cl3','cl4','cl5','cl6'))

pca <- read.table(esse_G1K_pruned_pca_exclude,sep = '\t', header = F)

pca <- pca[,c(2,3:12)]
colnames(pca)[1] <- 'NP'
pca <-  pca %>% dplyr::arrange(V3) %>% mutate(NP = as.factor(NP))

pop$NP <- factor(pop$NP,levels = pca$NP)


pop <- pop %>% gather(key,value,3:ncol(pop))


cols <- c('#009E73','#F0E442',"#999999","#56B4E9",'#0072B2','#E69F00',"#D55E00", "#CC79A7","#000000","#66FF33")
pop <- as.data.frame(merge(pop,pca,by='NP'))

pop %>% mutate(NP = factor(NP,levels = pca$NP)) %>%  dplyr::arrange(NP) %>%
  #mutate(key = case_when(key %in% c('ASW','ACB','LWK','GWD','YRI','MSL','ESN') ~ 'AFR',
  #                             key %in% c('PUR','CLM','MXL','PEL') ~ 'AMR',
  #                            key %in% c('PJL','GIH','BEB','ITU','STU') ~ 'SAS',
  #                            TRUE ~ key)) %>%
  #mutate(key = factor(key,levels=c('CEU','GBR','FIN','IBS','TSI',  'CHB','JPT','CHS','CDX','KHV', 'PJL','GIH','BEB','ITU','STU', 'PUR','CLM','MXL','PEL', 'AFR' ))) %>%
  #mutate(key = factor(key,levels=c('CEU','GBR','FIN','IBS','TSI',  'CHB','JPT','CHS','CDX','KHV', 'SAS', 'AMR', 'AFR' ))) %>%
  mutate(Population = key) %>%
  ggplot(aes(NP,value,col=Population,fill=Population))+
  geom_bar(stat='identity',pos=position_stack(reverse = FALSE))+
  #geom_point(aes(V3,0))+
  facet_wrap(~POP,scales = "free_x",ncol=2,nrow=3)+
  theme_classic()+
  theme(axis.text.x=element_blank()) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols)+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  theme(legend.position="none")
  #theme(legend.position="bottom", legend.box = "horizontal")


d <- pop %>% group_by(POP,key) %>% dplyr::summarise(N=mean(value),NN=sd(value))
d %>%
  ggplot(aes(POP,N,col=key,fill=key))+
  geom_bar(stat='identity',position=position_dodge(0.8), width=0.5)+
  theme_classic()+
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols)+
  geom_errorbar( aes(x=POP, ymin=N, ymax=N+NN,col=key),width=0.5, alpha=0.9,position=position_dodge(0.8))

pop %>%
  ggplot(aes(POP, value, fill=key))+
  geom_boxplot(outlier.size = 0.1,size = 0.2)+
  theme_classic()+
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols)
  

