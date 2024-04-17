library(gplots)
library(ggplot2)
library(tidyr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

folder <- 'working folder'


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

        svg(paste(folder,x,"_matrix.svg",sep=''),width = 10, height = 10)    
        my_palette <- colorRampPalette(c("#313695", "#FFFFBF", "#A50026"))(n = 100)

        heatmap.2(IBD_pop_matrix, scale="none", trace="none",col=my_palette)
	dev.off()

}

graph1  <- function(IBD_people,x) {

        png(filename=paste(folder,x,".png",sep=''),width = 720, height = 720)

        IBD_people %>%
        ggplot(aes(IBD))+
        geom_histogram(position='dodge')+
        theme_classic()
        #facet_wrap(~POP2,scales='free')

        dev.off()

}

x <- 'file pattern'

IBD_people <- read.table(paste(folder,x,sep=''),sep = '\t', header = F)
colnames(IBD_people) <- c('NP2','NP1','IBD','POP1','Super_POP1','POP2','Super_POP2')
IBD_people <- IBD_people %>% dplyr::arrange(desc(IBD))
print(head(IBD_people,200))
#IBD_people <- IBD_people %>% filter(IBD < 2000)
graph(IBD_people,paste(x,'_len',sep=''))

png(filename=paste(folder,x,"_len.png",sep=''),width = 720, height = 720)
IBD_people %>%
        ggplot(aes(IBD))+
        geom_histogram(position='dodge')+
        theme_classic()
        #facet_wrap(~POP2,scales='free')

dev.off()

IBD_people <- read.table(paste(folder,x, sep=''),sep = '\t', header = F)
colnames(IBD_people) <- c('NP2','NP1','IBD','POP1','Super_POP1','POP2','Super_POP2')
IBD_people <- IBD_people %>% dplyr::arrange(desc(IBD))
print(head(IBD_people,200))
graph(IBD_people,paste(x,'_count',sep=''))

png(filename=paste(folder,x,"_count.png",sep=''),width = 720, height = 720)
IBD_people %>%
        ggplot(aes(IBD))+
        geom_histogram(position='dodge')+
        theme_classic()
        #facet_wrap(~POP2,scales='free')

dev.off()

IBD_people <- read.table(paste(folder,x,sep=''),sep = '\t', header = F)
colnames(IBD_people) <- c('NP2','NP1','IBD','POP1','Super_POP1','POP2','Super_POP2')
IBD_people <- IBD_people %>% dplyr::arrange(desc(IBD))
print(head(IBD_people,200))
graph(IBD_people,paste(x,'_LOD',sep=''))

png(filename=paste(folder,x,sep=''),width = 720, height = 720)
IBD_people %>%
        ggplot(aes(IBD))+
        geom_histogram(position='dodge')+
        theme_classic()
        #facet_wrap(~POP2,scales='free')

dev.off()

