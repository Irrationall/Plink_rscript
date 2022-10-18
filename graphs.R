library(dplyr)
library(tidyverse)
library(ggplot2)
library(viridis)


# Gene with SNPs

df <- read_table("D:/서준서/OneDrive - 숭실대학교 - Soongsil University/3rd/2nd_semester/BS_exp/Results/LRP8/LRP8.chr1.sorted.filter0.00005")
df <- df %>%
  mutate(REGULATE = case_when(BETA > 0 ~ 'UP', BETA < 0 ~ 'DOWN'))

ggplot(data=df, aes(CHR,P)) + ggtitle("LRP8 related SNPs", ) + xlab("Chromosome") + ylab("P value") +
  geom_point(aes( fill=BETA, shape=REGULATE), size=10, color='black') +
  scale_shape_manual(values = c("UP" = 24, "DOWN" = 25)) +
  scale_fill_viridis(aesthetics = "fill") +
  scale_x_discrete(limits=c(1:21)) +
  scale_y_reverse() +
  theme(plot.title = element_text(size = 50, hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size=30),
        axis.text.x = element_text(size = 20, angle = 45),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size=20),
        panel.background = element_rect(fill = 'white', color = 'grey'),
        panel.grid.major = element_line(color = 'black', linetype = 'dotted'),
        panel.grid.minor = element_line(color = 'black', linetype = 'dotted')) 


# SNP with Genes

library(ggrepel)

df <- read.csv("D:/서준서/OneDrive - 숭실대학교 - Soongsil University/3rd/2nd_semester/BS_exp/Results/SNPs/rs12136984.csv")
df <- df[,-11]

df$X <- gsub(".assoc.linear", "", df$X)
df <- df %>% separate(X, c("Gene_chr","Gene_name"), sep = "\\.")
df <- df %>% mutate(logpval = -log10(P))
df$Gene_chr <- as.numeric(gsub("chr","", df$Gene_chr))
df <- df %>% mutate(outlier = case_when(logpval > 2.5 ~ T, logpval <= 2.5 ~ NA))

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=22)

ggplot(df,aes(x=factor(Gene_chr), y=logpval, fill=factor(Gene_chr))) + geom_violin() + geom_boxplot(width=0.1) + 
  ggtitle("Association of rs12136984 and genes") +
  geom_text_repel(data = subset(df, outlier == T), aes(label=Gene_name),force_pull = 3, size=6, max.overlaps = 20) +
  scale_fill_manual(values = color_list, guide='none') +
  theme(plot.title = element_text(size = 50, hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size=30),
        axis.text.x = element_text(size = 20, angle = 45),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size=20))
