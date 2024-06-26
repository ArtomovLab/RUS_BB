library(dplyr)
library(tidyr)
library(ggplot2)
'%!in%' <- function(x,y)!('%in%'(x,y))

ibdNE_file = 'file with ibdNE results'

for (i in c('RUS','CEU_GBR','IBS_TSI','ORENBURG','SAMARA','SPB','cl1','cl2','cl3','cl4','cl5','cl6','ASW','ACB','LWK','GWD','YRI','MSL','ESN','PUR','CLM','MXL','PEL','PJL','GIH','BEB','ITU','STU','CHB','JPT','CHS','CDX','KHV','CEU','GBR','FIN','IBS','TSI')) {
  if (i == 'RUS') {
    ibdNE_all <- read.table(paste(ibdNE_file,i,".ne",sep=''),sep = '\t', header = T)
    ibdNE_all <- ibdNE_all %>% dplyr::select(GEN,NE)
  }
  else {
    ibdNE <- read.table(paste(ibdNE_file,i,".ne",sep=''),sep = '\t', header = T)
    ibdNE <- ibdNE %>% dplyr::select(GEN,NE)
    ibdNE_all <- as.data.frame(merge(ibdNE_all,ibdNE,by='GEN'))
  }
}
colnames(ibdNE_all) <- c('GEN','RUS','CEU_GBR','IBS_TSI','ORENBURG','SAMARA','SPB','cl1','cl2','cl3','cl4','cl5','cl6','ASW','ACB','LWK','GWD','YRI','MSL','ESN','PUR','CLM','MXL','PEL','PJL','GIH','BEB','ITU','STU','CHB','JPT','CHS','CDX','KHV','CEU','GBR','FIN','IBS','TSI')
val <- ibdNE_all %>% gather(key,value,2:39) 


for (i in c('RUS','CEU_GBR','IBS_TSI','ORENBURG','SAMARA','SPB','cl1','cl2','cl3','cl4','cl5','cl6','ASW','ACB','LWK','GWD','YRI','MSL','ESN','PUR','CLM','MXL','PEL','PJL','GIH','BEB','ITU','STU','CHB','JPT','CHS','CDX','KHV','CEU','GBR','FIN','IBS','TSI')) {
  if (i == 'RUS') {
    ibdNE_all <- read.table(paste(ibdNE_file,i,".ne",sep=''),sep = '\t', header = T)
    ibdNE_all <- ibdNE_all %>% dplyr::select(GEN,LWR.95.CI)
  }
  else {
    i='FIN'
    ibdNE <- read.table(paste(ibdNE_file,i,".ne",sep=''),sep = '\t', header = T)
    ibdNE <- ibdNE %>% dplyr::select(GEN,LWR.95.CI)
    ibdNE_all <- as.data.frame(merge(ibdNE_all,ibdNE,by='GEN'))
  }
}
colnames(ibdNE_all) <- c('GEN','RUS','CEU_GBR','IBS_TSI','ORENBURG','SAMARA','SPB','cl1','cl2','cl3','cl4','cl5','cl6','ASW','ACB','LWK','GWD','YRI','MSL','ESN','PUR','CLM','MXL','PEL','PJL','GIH','BEB','ITU','STU','CHB','JPT','CHS','CDX','KHV','CEU','GBR','FIN','IBS','TSI')
low <- ibdNE_all %>% gather(key,value,2:39) 

for (i in c('RUS','CEU_GBR','IBS_TSI','ORENBURG','SAMARA','SPB','cl1','cl2','cl3','cl4','cl5','cl6','ASW','ACB','LWK','GWD','YRI','MSL','ESN','PUR','CLM','MXL','PEL','PJL','GIH','BEB','ITU','STU','CHB','JPT','CHS','CDX','KHV','CEU','GBR','FIN','IBS','TSI')) {
  if (i == 'RUS') {
    ibdNE_all <- read.table(paste(ibdNE_file,i,".ne",sep=''),sep = '\t', header = T)
    ibdNE_all <- ibdNE_all %>% dplyr::select(GEN,UPR.95.CI)
  }
  else {
    i='FIN'
    ibdNE <- read.table(paste(ibdNE_file,i,".ne",sep=''),sep = '\t', header = T)
    ibdNE <- ibdNE %>% dplyr::select(GEN,UPR.95.CI)
    ibdNE_all <- as.data.frame(merge(ibdNE_all,ibdNE,by='GEN'))
  }
}
colnames(ibdNE_all) <- c('GEN','RUS','CEU_GBR','IBS_TSI','ORENBURG','SAMARA','SPB','cl1','cl2','cl3','cl4','cl5','cl6','ASW','ACB','LWK','GWD','YRI','MSL','ESN','PUR','CLM','MXL','PEL','PJL','GIH','BEB','ITU','STU','CHB','JPT','CHS','CDX','KHV','CEU','GBR','FIN','IBS','TSI')
up <- ibdNE_all %>% gather(key,value,2:39) 

ibdNE_all <- cbind(val,low[,3])
ibdNE_all <- cbind(ibdNE_all,up[,3])
colnames(ibdNE_all ) <- c('GEN','key','value','LWR.95.CI','UPR.95.CI')

cols <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ibdNE_all %>% 
  filter(GEN<150) %>%
  filter(key %in% c('FIN','RUS')) %>%
  ggplot(aes(GEN,value,col=key))+
  geom_ribbon(aes(ymin=LWR.95.CI, ymax=UPR.95.CI,col=key), linetype=2, alpha=0.1)+
  geom_line()+
  theme_classic()+
  scale_y_log10()+
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols)
