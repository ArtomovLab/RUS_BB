library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(forcats)


'%!in%' <- function(x,y)!('%in%'(x,y))


results_all_fst <- data.frame(matrix(ncol=4,nrow=0))
colnames(results_all_fst) <- c('CHROM','REGION','key','value')
for (j in c('SAMARA','ORENBURG','SPB','RUS')) {
  #for (j in c('cl1','cl2','cl3','cl4','cl5','cl6')) {
  results_fst <- data.frame(1:29)
  colnames(results_fst) <- 'CHROM'
  for (i in c('ACB','ASW','BEB','CDX','CEU','CHB','CHS','CLM','ESN','FIN','GBR','GIH','GWD','IBS','ITU','JPT','KHV','LWK','MSL','MXL','PEL','PJL','PUR','STU','TSI','YRI','ORENBURG','SPB','SAMARA')) {
 
    print(c(i,j))
    samara <- read.table(paste("/humgen/atgu1/methods/dusoltsev/biobank/HRC/fst/regions_WGS_pruned/",j,"_",i,".weir.fst",sep=''),sep = '\t', header = T)
    samara <- samara %>% filter(!is.nan(WEIR_AND_COCKERHAM_FST))
    d <- samara %>% group_by(CHROM) %>% dplyr::summarise(N=mean(WEIR_AND_COCKERHAM_FST))
    colnames(d) <- c('CHROM',i)
    results_fst <- as.data.frame(merge(results_fst,d,by='CHROM'))
  }
  results_fst$REGION <- j
  results <- results_fst %>% gather(2:29,key=key,value=value) 
  results_all_fst <- rbind(results_all_fst,results)
}
results_all_fst <- results_all_fst %>% mutate(POP = case_when(key %in% c('CEU','GBR','FIN','IBS','TSI') ~ 'EUR',
                                                              key %in% c('PJL','GIH','BEB','ITU','STU') ~ 'SAS',
                                                              key %in% c('PUR','CLM','MXL','PEL') ~ 'AMR',
                                                              key %in% c('CHB','JPT','CHS','CDX','KHV') ~ 'EAS',
                                                              key %in% c('ASW','ACB','LWK','GWD','YRI','MSL','ESN') ~ 'AFR',
                                                              key %in% c('SPB','SAMARA','ORENBURG') ~ 'RUS'))
write.table(results_all_fst,"/humgen/atgu1/methods/dusoltsev/biobank/HRC/fst/regions_WGS_pruned.txt", sep='\t',quote = F,row.names = F,col.names = T)

##

library(ggpubr)
results_all_fst <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/fst/regions_WGS_pruned.txt",sep = '\t', header = T)
cols <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

g1 <- results_all_fst %>% filter(REGION !='RUS' & key %!in% c('SAMARA','ORENBURG','SPB')) %>%  
  dplyr::arrange(value) %>%
  mutate(key = fct_reorder(key, value, .desc = FALSE)) %>%
  filter(POP== 'EUR') %>%
  ggplot(aes(key,value,fill=REGION))+
  geom_boxplot(size =0.4,outlier.size=0.4)+
  theme_classic()+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  theme(legend.position = 'None')+
  ggtitle("EUR") +
  ylab('Fst')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

g2 <- results_all_fst %>% filter(REGION !='RUS' & key %!in% c('SAMARA','ORENBURG','SPB')) %>%  
  dplyr::arrange(value) %>% 
  mutate(key = fct_reorder(key, value, .desc = FALSE)) %>%
  filter(POP== 'AMR') %>%
  ggplot(aes(key,value,fill=REGION))+
  geom_boxplot(size =0.4,outlier.size=0.4)+
  theme_classic()+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  theme(legend.position = 'None')+
  ggtitle("AMR") +
  ylab('Fst')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

g3 <- results_all_fst %>% filter(REGION !='RUS' & key %!in% c('SAMARA','ORENBURG','SPB')) %>%  
  dplyr::arrange(value) %>% 
  mutate(key = fct_reorder(key, value, .desc = FALSE)) %>%
  filter(POP== 'SAS') %>%
  ggplot(aes(key,value,fill=REGION))+
  geom_boxplot(size =0.4,outlier.size=0.4)+
  theme_classic()+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  theme(legend.position = 'None')+
  ggtitle("SAS")+
  ylab('Fst')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

g4 <- results_all_fst %>% filter(REGION !='RUS' & key %!in% c('SAMARA','ORENBURG','SPB')) %>%  
  dplyr::arrange(value) %>% 
  mutate(key = fct_reorder(key, value, .desc = FALSE)) %>%
  filter(POP== 'EAS') %>%
  ggplot(aes(key,value,fill=REGION))+
  geom_boxplot(size =0.4,outlier.size=0.4)+
  theme_classic()+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  theme(legend.position = 'None')+
  ggtitle("EAS") +
  ylab('Fst')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 

g5 <- results_all_fst %>% filter(REGION !='RUS' & key %!in% c('SAMARA','ORENBURG','SPB')) %>%  
  dplyr::arrange(value) %>% 
  mutate(key = fct_reorder(key, value, .desc = FALSE)) %>%
  filter(POP== 'AFR') %>%
  ggplot(aes(key,value,fill=REGION))+
  geom_boxplot(size =0.4)+
  theme_classic()+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  theme(legend.position = 'None')+
  ggtitle("AFR") +
  ylab('Fst')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) 


grid.arrange(g1,g2,g3,g4,g5, ncol=3, nrow=2)
##################################

fst <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/fst/clusters6_WGS_pruned.txt",sep = '\t', header = T)

dd <- fst %>%  group_by(REGION,key,POP) %>% dplyr::summarise(N=mean(value)) %>% dplyr::arrange(N) %>% filter(POP != 'RUS') %>% dplyr::arrange(N) 

dd <- dd %>% 
  mutate(key=factor(key)) %>%
  group_by(POP) %>%
  dplyr::mutate(M = case_when(POP == 'EUR' ~ 1,
                              POP == 'AMR' ~ 2,
                              POP == 'SAS' ~ 3,
                              POP == 'EAS' ~ 4,
                              POP == 'AFR' ~ 5)) %>%
  ungroup() %>%
  mutate(key = reorder(key, M))


dd %>% 
  ggplot(aes(x=REGION,N,col=POP,fill=POP))+
  geom_bar(stat='identity',position='dodge')+
  facet_wrap(~key,scales='free',drop=F)+
  theme_classic()+
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ylab("Fst") 

