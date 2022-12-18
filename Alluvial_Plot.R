library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggalluvial)

df <- read.csv("brain.csv")
df <- df%>% mutate(LOGP = -log10(P.Value)) 

sig <- df[df$P.Value<0.05,]


## NetworkD3 #################################
links <- data.frame(
  source = c(df$Tissue, df$SNP),
  target = c(df$SNP, df$Gene.Symbol),
  value = c(df$LOGP, df$LOGP)
)

nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

sankeyNetwork(Links = links, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              sinksRight=FALSE)


## ggalluvial #################################

ggplot(data = sig2,
       aes(axis1 = Gene.Symbol, axis2 = Tissue, axis3 = SNP)) +
  scale_x_discrete(limits = c("Gene", "Tissue", "SNP"), expand = c(.05, .05)) +
  geom_alluvium(aes(fill = LOGP)) + ylab("aaa") +
  geom_stratum(size=20) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size=7) +
  theme_minimal() +
  theme(axis.title.y = element_text(colour = 'white'),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=30),
        legend.title = element_text(size=22),
        legend.text = element_text(size=18))+
  viridis::scale_fill_viridis(name='-log(P value)', direction = 1, option = "A")


                        