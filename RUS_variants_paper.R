library(dplyr)
library(tidyr)
library(ggplot2)
'%!in%' <- function(x,y)!('%in%'(x,y))

AF_RUS <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/AF_RUS.txt",sep = '\t', header = T)
AF_RUS <- AF_RUS %>% unite('locus', c('locus','alleles'),sep='_')
AF_RUS$AF <- gsub('\\[|\\]','',AF_RUS$AF)
AF_RUS <- AF_RUS %>% separate('AF',c('AF1','AF2'),sep=',')
AF_RUS <- AF_RUS %>% dplyr::select(locus,AF2)
AF_RUS$AF2 <- as.numeric(AF_RUS$AF2)

GNOMAD_filters <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/ANNOTATION_GNOMAD_filters.txt",sep = '\t', header = T)
GNOMAD_filters <- GNOMAD_filters %>% unite('locus', c('locus','alleles'),sep='_')

AF_RUS <- AF_RUS %>% filter(locus %!in% GNOMAD_filters$locus)

FIN_NEFIN <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/ANNOTATION_FREQ.txt",sep = '\t', header = T)
FIN_NEFIN <- FIN_NEFIN %>% unite('locus', c('locus','alleles'),sep='_')

dd_f <- AF_RUS %>% filter(locus %in% FIN_NEFIN$locus) 

FIN_NEFIN <- FIN_NEFIN %>% filter(locus %in% dd_f$locus)

d_f <- plyr::join(FIN_NEFIN,dd_f,by='locus')

d_f$locus <- gsub('\\[|\\]','',d_f$locus)
d_f$locus <- gsub('\\,','_',d_f$locus)
d_f$locus <- gsub('\\:','_',d_f$locus)

write.table(d_f,"/humgen/atgu1/methods/dusoltsev/biobank/HRC/AF_RUS_FIN_NEFIN.txt", sep='\t',quote = F,row.names = F)
###############################################################################################################################

d_f <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/AF_RUS_FIN_NEFIN.txt",sep = '\t', header = T)

d_f <- d_f %>% mutate(RUS_NEFIN = log(AF2/NEFIN,base=2))
d_ff <- d_f %>% filter(AF2> 0.01 & AF2 < 0.1)

d_ff %>% 
  ggplot(aes(RUS_NEFIN))+
  geom_density(col='black',alpha=0.5,adjust = 2,size = 0.5)+
  theme_classic()+
  scale_x_continuous()+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  geom_vline(xintercept = 0,linetype = "longdash",size = 0.5)


d_f_maf <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/MAF_RUS_FIN_NEFIN.txt",sep = '\t', header = T)
d_f_maf <- d_f_maf %>% unite('locus', c('locus','alleles'),sep='_')
d_f_maf <- d_f_maf  %>% filter(locus %!in% GNOMAD_filters$locus)
d_f_maf <- d_f_maf %>% drop_na()
d_f_maf$locus <- gsub(':|,','_',d_f_maf$locus)
d_f_maf$locus <- gsub('\\[|\\]','',d_f_maf$locus)
hwe <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/AF_RUShwe10e-4.txt",sep = '\t', header = T)
variants = read.table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/variants_to_exclude.txt',sep = '\t', header = T)
d_f_maf <- d_f_maf %>% filter(locus %!in% variants$rsid)
d_f_maf <- d_f_maf %>% filter(locus %!in% hwe$rsid)


fudgeit <- function(){
  xm <- get('xm', envir = parent.frame(1))
  ym <- get('ym', envir = parent.frame(1))
  z  <- get('dens', envir = parent.frame(1))
  colramp <- get('colramp', parent.frame(1))
  fields::image.plot(xm,ym,z, col = colramp(256), legend.only = T, add =F)
}

par(mar = c(5,4,4,5) + .1)


jet.colors <- colorRampPalette(c("white", "#e1eeffff", "#a5ccffff", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
#400 450
smoothScatter(x=d_f_maf$MAF_RUS, y=d_f_maf$MAF_NEFIN, colramp = jet.colors, xlab="MAF RUS", ylab="AF NFE GNOMAD", nbin = 128, nrpoints=500, pch=1, cex=0.1)

smoothScatter(x=d_f_maf$MAF_RUS, y=d_f_maf$MAF_FIN, colramp = jet.colors, xlab="MAF RUS", ylab="AF FIN GNOMAD", nbin = 128, nrpoints=500, pch=1, cex=0.1)

g1 <- d_f_maf  %>%
  ggplot(aes(MAF_RUS))+
  geom_histogram(color="white", fill="#4e4e4e",bins = 50)+
  theme_classic()


genotyped_variants <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/genotyped_variants.txt",sep = '\t', header = T)
genotyped_variants <- genotyped_variants %>% unite('locus', c('locus','alleles'),sep='_')
genotyped_variants <- genotyped_variants %>% filter(locus %!in% GNOMAD_filters$locus)
variants = read.table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/variants_to_exclude.txt',sep = '\t', header = T)
genotyped_variants <- genotyped_variants %>% filter(rsid %!in% variants$rsid)
genotyped_variants <- genotyped_variants %>% separate('AF',c('RSID','AF_FIN','AF_NEFIN'),sep=',')
genotyped_variants <- genotyped_variants %>% separate('variant_qc',c('AC1','AC2','AF1','AF2'),sep=',')
genotyped_variants$AF2 <- gsub('\\]','',genotyped_variants$AF2)
genotyped_variants$AF_FIN <- gsub('FIN:','',genotyped_variants$AF_FIN)
genotyped_variants$AF_NEFIN <- gsub('NEFIN:|\\}','',genotyped_variants$AF_NEFIN)
genotyped_variants$AF_RUS <- as.numeric(genotyped_variants$AF2)
genotyped_variants$AF_FIN <- as.numeric(genotyped_variants$AF_FIN)
genotyped_variants$AF_NEFIN <- as.numeric(genotyped_variants$AF_NEFIN)
genotyped_variants <- genotyped_variants %>% mutate(MAF_RUS_NEFIN = log(MAF_RUS/MAF_NEFIN,base=2))
genotyped_variants$locus <- gsub(':|,','_',genotyped_variants$locus)
genotyped_variants$locus <- gsub('\\[|\\]','',genotyped_variants$locus)

d_f_maf_genotyped <- d_f_maf %>% filter(locus %in% genotyped_variants$locus)

g2 <- d_f_maf_genotyped  %>%
  ggplot(aes(MAF_RUS))+
  geom_histogram(color="white", fill="#4e4e4e",bins = 50)+
  theme_classic()

grid.arrange(g2,g1, ncol=1, nrow=2)

smoothScatter(x=d_f_maf_genotyped$MAF_RUS, y=d_f_maf_genotyped$MAF_NEFIN, colramp = jet.colors, xlab="MAF RUS", ylab="AF NFE GNOMAD", nbin = 128, nrpoints=500, pch=1, cex=0.1)

smoothScatter(x=d_f_maf_genotyped$MAF_RUS, y=d_f_maf_genotyped$MAF_FIN, colramp = jet.colors, xlab="MAF RUS", ylab="AF FIN GNOMAD", nbin = 128, nrpoints=500, pch=1, cex=0.1)



AC_RUS <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/AC_AF001.txt",sep = '\t', header = T)
AC_RUS$AC <- gsub('\\[|\\]','',AC_RUS$AC)
AC_RUS <- AC_RUS %>% separate('AC',c('AC1','AC2'),sep=',')
AC_RUS$AC1 <- as.numeric(AC_RUS$AC1)
AC_RUS$AC2 <- as.numeric(AC_RUS$AC2)
AC_RUS <- AC_RUS %>% mutate(AC2 = case_when(AC2>5000 ~ AC1,
                                            TRUE ~ AC2))
AC_RUS <- AC_RUS %>% mutate(AC2_group = case_when(AC2 == 0 ~ '0',
                                                  AC2 <= 20 ~ '1-20',
                                                  AC2 <= 40 & AC2 > 20 ~ '21-40',
                                                  AC2 <= 60 & AC2 > 40 ~ '41-60',
                                                  AC2 <= 80 & AC2 > 60 ~ '61-80',
                                                  AC2 <= 100 & AC2 > 20 ~ '81-100'))
AC_RUS <- AC_RUS %>% unite('locus',c('locus','alleles'),sep='_')
AC_RUS$locus <- gsub(':|,','_',AC_RUS$locus)
AC_RUS$locus <- gsub('\\[|\\]','',AC_RUS$locus)

AC_RUS <- AC_RUS %>% filter(locus %!in% hwe$rsid)

AC_RUS <- AC_RUS %>% filter(locus %in% genotyped_variants$locus)
AC_RUS  %>% 
  ggplot(aes(AC2_group))+
  geom_bar(color="white", fill="#4e4e4e")+
  theme_classic()+
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))

##################################################################
d_f_maf <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/AF_RUS_FIN_NEFIN.txt",sep = '\t', header = T)
hwe <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/AF_RUShwe10e-4.txt",sep = '\t', header = T)
d_f_maf <- d_f_maf %>% filter(locus %!in% hwe$rsid)

d_f_maf <- d_f_maf %>% mutate(RUS_NEFIN = log(AF2/NEFIN,base=2))
d_f_maf <- d_f_maf %>% mutate(FIN_NEFIN = log(FIN/NEFIN,base=2))

d_ff_maf <- d_f_maf %>% filter(FIN > 0.01 & FIN < 0.1 & AF2 > 0.01)

d_ff_maf %>% 
  ggplot(aes(FIN_NEFIN))+
  geom_histogram(bins = 50,col='black',fill='#0072B2')+
  theme_classic()+
  scale_x_continuous()+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  geom_vline(xintercept = 1.5,linetype = "longdash",size = 0.5)

d_ff_maf %>% gather(key,value,5:6) %>%
  filter(value > -3 & value < 6) %>%
  ggplot(aes(value,col=key,fill=key))+
  geom_histogram(bins = 100,position='dodge')+
  theme_classic()

cols <- c('#0072B2',"#009E73","#999999","#E69F00",'#CC79A7',"#CCCCCC", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
d_ff_maf %>% filter(FIN_NEFIN > 2) %>%
  gather(key,value,5:6) %>%
  filter(value > -3 & value < 6) %>%
  ggplot(aes(value,col=key,fill=key))+
  geom_histogram(bins = 100,position='dodge')+
  theme_classic()+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  geom_vline(xintercept = 0,linetype = "longdash",size = 0.5)

ddd <- d_ff_maf %>% filter(FIN_NEFIN > 2)

write.table(ddd$locus,"/humgen/atgu1/methods/dusoltsev/biobank/HRC/finnish_enriched.txt",sep='\t',quote = F,row.names = F)
fin <-  read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/finnish_enriched_rsids.txt",sep = '\t', header = T)
write.table(fin$RSID,"/humgen/atgu1/methods/dusoltsev/biobank/HRC/finnish_enriched_rsids.txt",sep='\t',quote = F,row.names = F)

d %>% filter(RUS_NEFIN >= 2) %>%
  ggplot(aes(MAF_FIN,MAF_NEFIN))+
  geom_point()+
  theme_classic()

cols <- c("#009E73","#999999","#E69F00",'#CC79A7',"#CCCCCC", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

d_f %>% 
  ggplot(aes(RUS_NEFIN))+
  geom_density(col='black',alpha=0.5,adjust = 2,size = 0.5)+
  theme_classic()+
  scale_x_continuous()+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  geom_vline(xintercept = 0,linetype = "longdash",size = 0.5)

#Finnish Enriched
#########################################################################################################
d_ff_maf <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/AF_RUS_FIN_NEFIN_annotated_enriched2.txt",sep = '\t', header = T)
colnames(d_ff_maf)[4] <- 'AF_RUS'
colnames(d_ff_maf)[10] <- 'AF2'
d_ff_maf <-  d_ff_maf %>% mutate(AF1 = log(AF1/NEFIN,base=2))
d_ff_maf <-  d_ff_maf %>% mutate(AF2 = log(AF2/NEFIN,base=2))
d_ff_maf <-  d_ff_maf %>% mutate(AF3 = log(AF3/NEFIN,base=2))
d_ff_maf <-  d_ff_maf %>% mutate(AF4 = log(AF4/NEFIN,base=2))
d_ff_maf <-  d_ff_maf %>% mutate(AF5 = log(AF5/NEFIN,base=2))
d_ff_maf <-  d_ff_maf %>% mutate(AF6 = log(AF6/NEFIN,base=2))

d_ff_maf <-  d_ff_maf %>% mutate(AF_SPB = log(AF_SPB/NEFIN,base=2))
d_ff_maf <-  d_ff_maf %>% mutate(AF_ORENBURG = log(AF_ORENBURG/NEFIN,base=2))
d_ff_maf <-  d_ff_maf %>% mutate(AF_SAMARA = log(AF_SAMARA/NEFIN,base=2))

d_ff_maf <-  d_ff_maf %>% mutate(JEW = log(JEW/NEFIN,base=2))
d_ff_maf <-  d_ff_maf %>% mutate(AFRICAN = log(AFRICAN/NEFIN,base=2))
d_ff_maf <-  d_ff_maf %>% mutate(EAST_ASIAN = log(EAST_ASIAN/NEFIN,base=2))
d_ff_maf <-  d_ff_maf %>% mutate(LATINO = log(LATINO/NEFIN,base=2))
d_ff_maf <-  d_ff_maf %>% mutate(EUR_OTH = log(EUR_OTH/NEFIN,base=2))
d_ff_maf[d_ff_maf == Inf] <- 10
d_ff_maf[d_ff_maf == -Inf] <- -5

##############################################################################################################

dd <-  d_ff_maf %>% filter(FIN_NEFIN > 2) %>%
  filter(EAST_ASIAN == -5) %>%
  gather(key,value,7:14,21) %>% 
  filter(value != 10 & value != -5) %>% 
  group_by(key) %>% dplyr::summarise(N=median(value,na.rm=T))%>% arrange(desc(N)) %>% ungroup() %>%
  mutate(key = case_when(key == 'FIN_NEFIN' ~ 'FIN',
                         key == 'RUS_NEFIN' ~ 'RUS',
                         key == 'AF1' ~ 'RUS_cl1',
                         key == 'AF2' ~ 'RUS_cl2',
                         key == 'AF3' ~ 'RUS_cl3',
                         key == 'AF4' ~ 'RUS_cl4',
                         key == 'AF5' ~ 'RUS_cl5',
                         key == 'AF6' ~ 'RUS_cl6',
                         key == 'EAST_ASIAN' ~ 'EAST_ASIAN')) %>%
  mutate(key=factor(key,levels = c('FIN','RUS','RUS_cl1','RUS_cl2','RUS_cl3','RUS_cl4','RUS_cl5','RUS_cl6','EAST_ASIAN')))

cols <- c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

d_ff_maf %>% filter(FIN_NEFIN > 2) %>%   
  #filter(EAST_ASIAN != -5) %>%
  gather(key,value,7:14,21) %>%
  filter(value != 10 & value != -5) %>%
  mutate(key = case_when(key == 'FIN_NEFIN' ~ 'FIN',
                         key == 'RUS_NEFIN' ~ 'RUS',
                         key == 'AF1' ~ 'RUS_cl1',
                         key == 'AF2' ~ 'RUS_cl2',
                         key == 'AF3' ~ 'RUS_cl3',
                         key == 'AF4' ~ 'RUS_cl4',
                         key == 'AF5' ~ 'RUS_cl5',
                         key == 'AF6' ~ 'RUS_cl6',
                         key == 'EAST_ASIAN' ~ 'EAST_ASIAN')) %>%
  mutate(key=factor(key,levels = c('FIN','RUS','RUS_cl1','RUS_cl2','RUS_cl3','RUS_cl4','RUS_cl5','RUS_cl6','EAST_ASIAN'))) %>%
  ggplot(aes(value,fill=key))+
  geom_density(alpha=0.5)+
  facet_wrap(~key,ncol = 1)+
  theme_classic()+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  geom_vline(xintercept = 0,linetype = "longdash",size = 0.5)+
  geom_vline(data = dd, aes(xintercept = N),linetype = "longdash",size = 0.25,col='red')+
  theme(axis.text.x = element_text(angle = 90))+
  scale_x_continuous('log2(AF/AF_NFE)')+
  scale_y_continuous('density')+
  theme(legend.position="none")


d_ff_maf %>% filter(FIN_NEFIN > 2) %>%   
  #filter(EAST_ASIAN == -5) %>%
  gather(key,value,7:14,21) %>%
  filter(value != 10 & value != -5) %>%
  mutate(key = case_when(key == 'FIN_NEFIN' ~ 'FIN',
                         key == 'RUS_NEFIN' ~ 'RUS',
                         key == 'AF1' ~ 'RUS_cl1',
                         key == 'AF2' ~ 'RUS_cl2',
                         key == 'AF3' ~ 'RUS_cl3',
                         key == 'AF4' ~ 'RUS_cl4',
                         key == 'AF5' ~ 'RUS_cl5',
                         key == 'AF6' ~ 'RUS_cl6',
                         key == 'EAST_ASIAN' ~ 'EAST_ASIAN')) %>%
  mutate(key=factor(key,levels = c('FIN','RUS','RUS_cl1','RUS_cl2','RUS_cl3','RUS_cl4','RUS_cl5','RUS_cl6','EAST_ASIAN'))) %>%
  ggplot(aes(key,value,fill=effect))+
  #geom_density(alpha=0.5)+
  geom_boxplot(outlier.size=0.1)+
  #facet_wrap(~key,ncol = 1)+
  theme_classic()+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  #geom_vline(xintercept = 0,linetype = "longdash",size = 0.5)+
  #geom_vline(data = dd, aes(xintercept = N),linetype = "longdash",size = 0.25,col='red')+
  theme(axis.text.x = element_text(angle = 90))
  #scale_x_continuous('log2(AF/AF_NFE)')+
  #scale_y_continuous('density')+
  theme(legend.position="none")
#######################################################################################################

#RUS enriched
AF_log2 <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/AF_log0_annotated.txt",sep = '\t', header = T)
AF_log2 <- AF_log2 %>% mutate(effect = case_when(Consequence %in% c('transcript_ablation',
                                                                    'splice_acceptor_variant',
                                                                    'splice_donor_variant',
                                                                    'stop_gained',
                                                                    'frameshift_variant'
) ~ 'PTV',
Consequence %in% c('missense_variant') ~ 'missense',
Consequence %in% c('synonymous_variant') ~ 'synonymous',
TRUE ~ 'other'))
AF_log2 <- AF_log2  %>% mutate(RUS_NEFIN = log(AF_RUS/NEFIN,base=2))

genotyped_variants <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/genotyped_variants.txt",sep = '\t', header = T)
genotyped_variants <- genotyped_variants %>% unite('locus', c('locus','alleles'),sep='_')

AF_cl <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/AF_RUS_cl_log0.txt",sep = '\t', header = T)
AF_reg <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/AF_RUS_REGIONS_log0.txt",sep = '\t', header = T)

AF_log2_ <- AF_log2 %>% mutate(IMP = case_when(locus %in% genotyped_variants$rsid ~ 'IMP_no',
                                               locus %!in% genotyped_variants$rsid ~ 'IMP_yes'))
AF_log2_ <- AF_log2_ %>% filter(RUS_NEFIN >2) 
AF_log2_ <- as.data.frame(merge(AF_log2_,AF_cl,by='locus'))
AF_log2_ <- as.data.frame(merge(AF_log2_,AF_reg,by='locus'))
####
AF_log2_ <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/AF_log2_.txt",sep = '\t', header = T)

freq <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/EGDP/EGDP_rusEUR.afreq",sep = '\t', header = F)
freq <- freq %>% dplyr::select(V2,V5)
colnames(freq) <- c('locus','AF_EGDP')

AF_log2_ <- as.data.frame(merge(AF_log2_,freq,by='locus',all.x=T))

freq_block <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/ECCE_BLOCK/block_af.txt",sep = '\t', header = T)
colnames(freq_block) <- c('locus','alleles','AF')
freq_block <- freq_block %>% unite('locus', c('locus','alleles'),sep='_')
freq_block$AF <- gsub('\\[|\\]','',freq_block$AF)
freq_block <- freq_block %>% separate('AF',c('AF1','AF2'),sep=',')
freq_block <- freq_block %>% dplyr::select(locus,AF2)
freq_block$AF2 <- as.numeric(freq_block$AF2)
freq_block$locus <- gsub('\\[|\\]','',freq_block$locus)
freq_block$locus <- gsub(',|:','_',freq_block$locus)
freq_block$AF1 <- NULL
colnames(freq_block)[2] <- 'AF_block'

AF_log2_ <- as.data.frame(merge(AF_log2_,freq_block,by='locus',all.x=T))
AF_log2_ <- AF_log2_  %>% filter(!is.na(AF_block))

AF_log2_ <- AF_log2_ %>% mutate(AF1 = log(AF1/NEFIN,base=2),
                                AF2 = log(AF2/NEFIN,base=2),
                                AF3 = log(AF3/NEFIN,base=2),
                                AF4 = log(AF4/NEFIN,base=2),
                                AF5 = log(AF5/NEFIN,base=2),
                                AF6 = log(AF6/NEFIN,base=2),
                                FIN = log(FIN/NEFIN,base=2),
                                EAST_ASIAN = log(EAST_ASIAN/NEFIN,base=2))
AF_log2_[AF_log2_ == Inf] <- 10
AF_log2_[AF_log2_ == -Inf] <- -5


d <- AF_log2_ %>% filter(EAST_ASIAN == -5 & AF5 == -5)
  

dd <- AF_log2_ %>% 
  filter(EAST_ASIAN == -5) %>%
  gather(c(3,6,12,14:19),key=key,value=value) %>%
  filter(value != 10 & value != -5) %>% 
  group_by(key) %>% dplyr::summarise(N=median(value,na.rm=T)) %>% drop_na() %>% arrange(desc(N)) %>% ungroup() %>%
  mutate(key = case_when(key == 'FIN' ~ 'FIN',
                         key == 'RUS_NEFIN' ~ 'RUS',
                         key == 'AF1' ~ 'RUS_cl1',
                         key == 'AF2' ~ 'RUS_cl2',
                         key == 'AF3' ~ 'RUS_cl3',
                         key == 'AF4' ~ 'RUS_cl4',
                         key == 'AF5' ~ 'RUS_cl5',
                         key == 'AF6' ~ 'RUS_cl6',
                         key == 'EAST_ASIAN' ~ 'EAST_ASIAN')) %>%
  mutate(key=factor(key,levels = c('FIN','RUS','RUS_cl1','RUS_cl2','RUS_cl3','RUS_cl4','RUS_cl5','RUS_cl6','EAST_ASIAN')))


cols <- c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7","#D55E00")
AF_log2_ %>%  
  #filter(EAST_ASIAN == -5) %>%
  gather(c(3,6,12,14:19),key=key,value=value) %>%

  mutate(key = case_when(key == 'FIN' ~ 'FIN',
                         key == 'RUS_NEFIN' ~ 'RUS',
                         key == 'AF1' ~ 'RUS_cl1',
                         key == 'AF2' ~ 'RUS_cl2',
                         key == 'AF3' ~ 'RUS_cl3',
                         key == 'AF4' ~ 'RUS_cl4',
                         key == 'AF5' ~ 'RUS_cl5',
                         key == 'AF6' ~ 'RUS_cl6',
                         key == 'EAST_ASIAN' ~ 'EAST_ASIAN')) %>%
  mutate(key=factor(key,levels = c('FIN','RUS','RUS_cl1','RUS_cl2','RUS_cl3','RUS_cl4','RUS_cl5','RUS_cl6','EAST_ASIAN'))) %>%
  ggplot(aes(value,fill=key))+
  facet_wrap(~key,ncol=1)+
  geom_density(alpha=0.5)+
  theme_classic()+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  #theme(axis.text.x = element_text(angle = 90))+
  geom_vline(xintercept = 0,linetype = "longdash",size = 0.5)+
  geom_vline(data = dd, aes(xintercept = N),linetype = "longdash",size = 0.25,col='red')+
  scale_x_continuous('log2(AF/AF_NFE)')+
  scale_y_continuous('density')+
  theme(legend.position="none")

AF_log2_ %>%  
  #filter(EAST_ASIAN == -5) %>%
  gather(c(3,6,12,14:19),key=key,value=value) %>%
  
  filter(value != 10 & value != -5) %>%
  mutate(key = case_when(key == 'FIN' ~ 'FIN',
                         key == 'RUS_NEFIN' ~ 'RUS',
                         key == 'AF1' ~ 'RUS_cl1',
                         key == 'AF2' ~ 'RUS_cl2',
                         key == 'AF3' ~ 'RUS_cl3',
                         key == 'AF4' ~ 'RUS_cl4',
                         key == 'AF5' ~ 'RUS_cl5',
                         key == 'AF6' ~ 'RUS_cl6',
                         key == 'EAST_ASIAN' ~ 'EAST_ASIAN')) %>%
  mutate(key=factor(key,levels = c('FIN','RUS','RUS_cl1','RUS_cl2','RUS_cl3','RUS_cl4','RUS_cl5','RUS_cl6','EAST_ASIAN'))) %>%
  ggplot(aes(key,value,fill=effect))+
  #geom_density(alpha=0.5)+
  geom_boxplot(outlier.size=0.1)+
  facet_wrap(~IMP,ncol = 1)+
  theme_classic()+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)+
  #geom_vline(xintercept = 0,linetype = "longdash",size = 0.5)+
  #geom_vline(data = dd, aes(xintercept = N),linetype = "longdash",size = 0.25,col='red')+
  theme(axis.text.x = element_text(angle = 90))
  #scale_x_continuous('log2(AF/AF_NFE)')+
  #scale_y_continuous('density')
  theme(legend.position="none")



