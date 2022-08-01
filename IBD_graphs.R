library(gplots)
library(ggplot2)
library(tidyr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

graph  <- function(IBD_people,x) {
        

	IBD_pop <- IBD_people %>% group_by(POP1,POP2) %>% dplyr::summarise(IBD_mean=median(IBD))

	for (i in 1:nrow(IBD_pop)) {
  	a <- IBD_pop$POP1[i]
  	b <- IBD_pop$POP2[i]
  	test <- IBD_pop %>% filter(POP2 == a & POP1 == b)
	print(test)
  	if (nrow(test) > 0) {
    	IBD_pop$POP2[i] <- a
    	IBD_pop$POP1[i] <- b
  		}
	}

	IBD_pop <- IBD_pop %>% group_by(POP1,POP2) %>% dplyr::summarise(IBD_mean=mean(IBD_mean))
	IBD_pop_dop <- IBD_pop
	colnames(IBD_pop_dop) <- c('POP2','POP1','IBD_mean')
	IBD_pop <- rbind(IBD_pop,IBD_pop_dop)

	#IBD_pop <- IBD_pop %>% filter(POP1 != POP2)

	IBD_pop <- unique(IBD_pop)
	IBD_pop_matrix <- pivot_wider(IBD_pop, names_from = POP1, values_from = IBD_mean)
	row.names(IBD_pop_matrix ) <- IBD_pop_matrix$POP2
	a <- IBD_pop_matrix$POP2
	IBD_pop_matrix$POP2 <- NULL
	row.names(IBD_pop_matrix ) <- a
	IBD_pop_matrix <- as.matrix(IBD_pop_matrix )

	IBD_pop_matrix <- log(IBD_pop_matrix,base=2)

	#png(filename=paste("/humgen/atgu1/methods/dusoltsev/biobank/populations/IBD/",x,"_matrix.png",sep=''),width = 720, height = 720)
        svg(paste("/humgen/atgu1/methods/dusoltsev/biobank/populations/IBD/",x,"_matrix.svg",sep=''),width = 10, height = 10)    
        my_palette <- colorRampPalette(c("#313695", "#FFFFBF", "#A50026"))(n = 100)
        #rev(brewer.pal(11, 'RdYlBu'))
	#heatmap(IBD_pop_matrix, scale="none")
        heatmap.2(IBD_pop_matrix, scale="none", trace="none",col=my_palette)
	#pheatmap(IBD_pop_matrix,cutree_rows = 4)
	dev.off()

}

graph1  <- function(IBD_people,x) {

        png(filename=paste("/humgen/atgu1/methods/dusoltsev/biobank/populations/IBD/",x,".png",sep=''),width = 720, height = 720)

        IBD_people %>%
        ggplot(aes(IBD))+
        geom_histogram(position='dodge')+
        theme_classic()
        #facet_wrap(~POP2,scales='free')

        dev.off()

}

x='IBD_1000G_genotype_merge'

IBD_people <- read.table(paste("/humgen/atgu1/methods/dusoltsev/biobank/populations/IBD/",x,"_annotated_len.ibd",sep=''),sep = '\t', header = F)
colnames(IBD_people) <- c('NP2','NP1','IBD','POP1','Super_POP1','POP2','Super_POP2')
IBD_people <- IBD_people %>% dplyr::arrange(desc(IBD))
print(head(IBD_people,200))
#IBD_people <- IBD_people %>% filter(IBD < 2000)
graph(IBD_people,paste(x,'_len',sep=''))

png(filename=paste("/humgen/atgu1/methods/dusoltsev/biobank/populations/IBD/",x,"_len.png",sep=''),width = 720, height = 720)
IBD_people %>%
        ggplot(aes(IBD))+
        geom_histogram(position='dodge')+
        theme_classic()
        #facet_wrap(~POP2,scales='free')

dev.off()

IBD_people <- read.table(paste("/humgen/atgu1/methods/dusoltsev/biobank/populations/IBD/",x,"_annotated_count.ibd", sep=''),sep = '\t', header = F)
colnames(IBD_people) <- c('NP2','NP1','IBD','POP1','Super_POP1','POP2','Super_POP2')
IBD_people <- IBD_people %>% dplyr::arrange(desc(IBD))
print(head(IBD_people,200))
graph(IBD_people,paste(x,'_count',sep=''))

png(filename=paste("/humgen/atgu1/methods/dusoltsev/biobank/populations/IBD/",x,"_count.png",sep=''),width = 720, height = 720)
IBD_people %>%
        ggplot(aes(IBD))+
        geom_histogram(position='dodge')+
        theme_classic()
        #facet_wrap(~POP2,scales='free')

dev.off()

IBD_people <- read.table(paste("/humgen/atgu1/methods/dusoltsev/biobank/populations/IBD/",x,"_annotated_LOD.ibd",sep=''),sep = '\t', header = F)
colnames(IBD_people) <- c('NP2','NP1','IBD','POP1','Super_POP1','POP2','Super_POP2')
IBD_people <- IBD_people %>% dplyr::arrange(desc(IBD))
print(head(IBD_people,200))
graph(IBD_people,paste(x,'_LOD',sep=''))

png(filename=paste("/humgen/atgu1/methods/dusoltsev/biobank/populations/IBD/",x,"_LOD.png",sep=''),width = 720, height = 720)
IBD_people %>%
        ggplot(aes(IBD))+
        geom_histogram(position='dodge')+
        theme_classic()
        #facet_wrap(~POP2,scales='free')

dev.off()

