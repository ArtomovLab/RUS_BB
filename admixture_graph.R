library(dplyr)
library(tidyr)
library(ggplot2)
'%!in%' <- function(x,y)!('%in%'(x,y))

fam <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/ADMIXTURE/pop9/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS1.fam")
table(fam$POP)
s <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/ADMIXTURE/pop9/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS1.pop")
fam <- cbind(fam,s)
fam <- fam[,c(2,7)]
colnames(fam) <- c('NP','POP')


cl1 <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/fst/esse.hrc.dr08_pruned_PLINK_CLUSTERS_61new.txt",sep = '\t', header = T)
cl2 <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/fst/esse.hrc.dr08_pruned_PLINK_CLUSTERS_62new.txt",sep = '\t', header = T)
cl3 <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/fst/esse.hrc.dr08_pruned_PLINK_CLUSTERS_63new.txt",sep = '\t', header = T)
cl4 <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/fst/esse.hrc.dr08_pruned_PLINK_CLUSTERS_64new.txt",sep = '\t', header = T)
cl5 <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/fst/esse.hrc.dr08_pruned_PLINK_CLUSTERS_65new.txt",sep = '\t', header = T)
cl6 <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/fst/esse.hrc.dr08_pruned_PLINK_CLUSTERS_66new.txt",sep = '\t', header = T)

fam <- fam %>% mutate(POP = case_when(NP %in% cl1$x ~ 'cl1',
                                                           NP %in% cl2$x ~ 'cl2',
                                                           NP %in% cl3$x ~ 'cl3',
                                                           NP %in% cl4$x ~ 'cl4',
                                                           NP %in% cl5$x ~ 'cl5',
                                                           NP %in% cl6$x ~ 'cl6',
                                                           TRUE ~ POP))

table(fam$POP)

pop <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/ADMIXTURE/pop9/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS1.8.Q")
pop <- cbind(fam,pop)
colnames(pop) <- c('NP','POP','CEU_GBR','FIN','CHB_CHS_JPT','AMR','CDX_KHV','IBS_TSI','SAS','AFR')

pop <- pop %>% filter(POP %in% c('cl1','cl2','cl3','cl4','cl5','cl6'))

pca <- read.table(paste("/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.hrc.dr08_pruned_PLINK2_",'final',".eigenvec",sep=''),sep = '\t', header = F)
pca <- pca[,c(2,3:12)]
colnames(pca)[1] <- 'NP'
pca <-  pca %>% dplyr::arrange(V3) %>% mutate(NP = as.factor(NP))

pop$NP <- factor(pop$NP,levels = pca$NP)


pop <- pop %>% gather(key,value,3:ncol(pop))


cols <- c('#009E73','#F0E442',"#999999","#56B4E9",'#0072B2','#E69F00',"#D55E00", "#CC79A7","#000000","#66FF33")


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
  facet_wrap(~POP,scales = "free_x",ncol=2,nrow=3)+
  theme_classic()+
  theme(axis.text.x=element_blank()) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols)+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())+
  theme(legend.position="none")
  #theme(legend.position="bottom", legend.box = "horizontal")

