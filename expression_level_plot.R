# mRNA expression box plot code
library(dplyr)
library(ggplot2)
library(patchwork)

mRNA <- read.table('EXOSC10_expression.ped')
names(mRNA) <- c('FID','IID','PID','MID','Sex','P','A1','A2')
mRNA <- mRNA[mRNA$P>=0,]
data <- mRNA %>% select(P)
data <- data %>% mutate(A=paste(mRNA$A1,mRNA$A2))

boxplot(P~A, data = data, main='EXOSC10')
xlab("Minor Allele")
ylab('Expression Level')


origin <- read.table('LRP8_expression.ped')
names(origin) <- c('FID','IID','PID','MID','Sex','P','A1','A2')
origin <- origin[origin$P>=0,]
data2 <- origin %>% select(P)
data2 <- data2 %>% mutate(A=paste(mRNA$A1,mRNA$A2))


p1 <- ggplot(data, aes(x=A,y=P, fill=factor(A))) + geom_violin(show.legend = F) + geom_boxplot(width=0.1, show.legend = F) +
  ggtitle("EXOSC10") + scale_fill_brewer(palette = "YlOrRd", guide='none') + 
  xlab('Allele') + ylab('Expression Level') +
  guides(fill = guide_legend(title = "Minor Allele")) +
  theme(legend.title = element_text(size = 25),
        axis.title = element_text(size=30),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 40))

p2 <- ggplot(data2, aes(x=A,y=P, fill=factor(A))) + geom_violin() + geom_boxplot(width=0.1, show.legend = F) +
  ggtitle("LRP8") + scale_fill_brewer(palette = 'YlOrRd') + 
  xlab('Allele') + ylab('Expression Level') +
  guides(fill = guide_legend(title = paste("Minor","\n","Allele"))) +
  theme(legend.title = element_text(size = 25),
        axis.title = element_text(size=30),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 40),
        axis.title.y = element_blank(),
        legend.key.size = unit(3,"line"))


p1+p2






