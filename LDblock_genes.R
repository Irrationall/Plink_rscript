library(dplyr)
library(tidyverse)
library(ggplot2)
library(viridis)
library(pheatmap)

# Load association file between LD block SNPs and a gene
# It is convenient to make association file name as 'Genename_Chromsome_blockname.assoc.linear'


#### Method 1 - NOGADA

ACSM1_block1 <- read_table("ACSM1_chr16_block1.assoc.linear")
ACSM1_block2 <- read_table("ACSM1_chr16_block2.assoc.linear")
DNAH3_block2 <- read_table("DNAH3_chr16_block2.assoc.linear")
THUMPD1_block1 <-  read_table("THUMPD1_chr16_block1.assoc.linear")
THUMPD1_block2 <-  read_table("THUMPD1_chr16_block2.assoc.linear")

dfs <- Filter(function(x) is(x, "data.frame"), mget(ls()))

mergelist(dfs)
Gene_Snp_Heatmap(LDblockdf)


#### Method 2 - if assoc.linear files in one directory

filelist = list.files(path = getwd(), pattern = "assoc.linear")
datalist = lapply(filelist, function(x)read_table(x)) 
names(datalist) <- str_remove_all(filelist, paste(c("chr16_",".assoc.linear"), collapse = "|"))

chr16_LDblock <- mergelist(datalist)

LDblockdf <- arrange(chr16_LDblock, BP)

LDblockdf <- LDblockdf %>% mutate(LOGP = -log10(P))

Gene_Snp_Heatmap(chr16_LDblock, "Chr16")


#### Function 1 - merge list to dataframe ####

mergelist <- function(dflist) {
  for (i in 1:length(dflist)) {
    dflist[[i]]$X10 <- names(dflist[i])
    dflist[[i]] <- dflist[[i]] %>% separate(X10, c("GENE","BLOCK"), sep = "_")
  }
  LDblockdf <- do.call(rbind, dflist)
  rownames(LDblockdf) <- NULL
  return(LDblockdf)
}
####--------------------------------------####


#### Function 2 - make a heatmap between genes and SNPs ####

# Heatmap with ggplot

Gene_Snp_Heatmap <- function(df, title = "Title") {
  p1 <- ggplot(df, aes(x=fct_inorder(SNP), y=GENE, fill=P, order=BP)) + geom_tile() + ggtitle(paste(title)) +
    scale_fill_viridis(values=scales::rescale(c(0,0.1,1)),direction = -1, guide = guide_colorbar(barheight = 10, barwidth = 2)) + 
    xlab('SNP') + ylab('GENE') +
    facet_grid(~BLOCK, scales="free_x", space = "free_x") +
    theme(axis.text.x = element_text(angle = 90, size = 10),
          plot.title = element_text(size = 50, hjust = 0.5),
          axis.text.y = element_text(size = 20),
          axis.title = element_text(size=30),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(0.1, units = 'cm'),
          panel.background = element_rect(fill='white'),
          strip.background.x = element_rect(fill = c("grey")),
          strip.text = element_text(size = 10, lineheight = 0.2)) 
  
  ggplotColours <- function(n = 6, h = c(0, 360) + 15){
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
  }
  
  color_list <- ggplotColours(n=length(levels(factor(df$BLOCK))))
  
  g <- ggplot_gtable(ggplot_build(p1))
  stripr <- which(grepl('strip-t', g$layout$name))
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- color_list[k]
    k <- k+1
  }
  grid::grid.draw(g)
}
####----------------------------------------------------####

# For quick test of plot

ggplot(LDblockdf, aes(x=fct_inorder(SNP), y=GENE, fill=P, order=BP)) + geom_tile() + ggtitle("test") +
  scale_fill_viridis_c(values=scales::rescale(c(0,0.1,1)),direction = -1, guide = guide_colorbar(barheight = 10, barwidth = 2)) + xlab('SNP') + ylab('GENE') +
  facet_grid(~BLOCK, scales="free_x", space = "free_x") +
  theme(axis.text.x = element_text(angle = 90, size = 10),
        plot.title = element_text(size = 50, hjust = 0.5),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size=30),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.1, units = 'cm'),
        panel.background = element_rect(fill='white'),
        strip.background.x = element_rect(fill = c("grey")),
        strip.text = element_text(size = 10, lineheight = 0.2)) 

# Make matrix

mat <- with(LDblockdf, tapply(P, list(SNP, GENE), FUN = print))
mat[is.na(mat)] <- 1