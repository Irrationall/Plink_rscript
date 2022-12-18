library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggrepel)





####### Manhattan Plot ########

# Original source code: https://r-graph-gallery.com/101_Manhattan_plot.html
# Plink result file (assoc.linear, assoc etc.) needed.

Manhattan_custom <- function(df, 
                             threshold_line = 0.00000005, 
                             label = FALSE,
                             pt.size = 2,
                             label.size = 5,
                             oddcolor = "tomato",
                             evencolor = "midnightblue") {
  library(qqman)
  
  newdf <- df %>% 
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
    select(-chr_len) %>%
    left_join(df, ., by=c("CHR"="CHR")) %>%
    arrange(CHR, BP) %>%
    mutate(BPcum=BP+tot) %>%
    mutate(is_annotate=ifelse(-log10(P)>6, "yes", "no")) 
  
  axisdf = newdf %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  if (label == TRUE) {
    ggplot(newdf, aes(x=BPcum, y=-log10(P))) +
      geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=pt.size) +
      scale_color_manual(values = rep(c(oddcolor, evencolor), 22 )) +
      scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
      scale_y_continuous(expand = c(0, 0) ) +
      theme_bw() +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
      ) +
      geom_hline(yintercept=-log10(0.000005), linetype="dashed", 
                 color = "wheat4", size=2.5) +
      geom_hline(yintercept=-log10(threshold_line), linetype="dashed", 
                 color = "wheat4", size=2.5) +
      xlab("Chromosome") + ylab("-log10(P value)") +
      geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=label.size, max.overlaps = 50)
  } else {
    ggplot(newdf, aes(x=BPcum, y=-log10(P))) +
      geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=pt.size) +
      scale_color_manual(values = rep(c(oddcolor, evencolor), 22 )) +
      scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
      scale_y_continuous(expand = c(0, 0) ) +
      theme_bw() +
      theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
      ) +
      geom_hline(yintercept=-log10(0.000005), linetype="dashed", 
                 color = "wheat4", size=2.5) +
      geom_hline(yintercept=-log10(threshold_line), linetype="dashed", 
                 color = "wheat4", size=2.5) +
      xlab("Chromosome") + ylab("-log10(P value)")
  }
}




####### Visualize Gene expression depend on minor allele ########

# Plink result file (assoc.linear, assoc etc.) needed.
# Viridis package required.
# Only for 22 chromosomes. (No sex chromosome)

Expression_Plot <- function(df, pt.size=10) {
  
  df <- df %>%
    mutate(REGULATE = case_when(BETA > 0 ~ 'UP', BETA < 0 ~ 'DOWN'), logpval = -log10(P)) 
  
  ggplot(data=df, aes(CHR,logpval)) + xlab("Chromosome") + ylab("-log10(P value)") +
    geom_point(aes( fill=BETA, shape=REGULATE), size=pt.size, color='black') +
    scale_shape_manual(values = c("UP" = 24, "DOWN" = 25)) +
    viridis::scale_fill_viridis(aesthetics = "fill",guide = guide_colorbar(barheight = 10, barwidth = 2)) +
    scale_x_discrete(limits=c(1:22))
}




####### Visualize association with one SNP and all genes ########

# Plink result file (assoc.linear, assoc etc.) needed.
# Set Chromosome as column number contains chromosome information.
# Set Genename as column number contains Gene symbols.
# Each point indicates GENE symbol.

SingleSNP_Plot <- function(df, Chromosome, Genename, label = FALSE) {
  
  df <- df %>%
    mutate(logpval = -log10(P)) %>%
    mutate(outlier = case_when(logpval > 2.5 ~ T, logpval <= 2.5 ~ NA))
  
  ggplotColours <- function(n = 6, h = c(0, 360) + 15){
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
  }
  
  color_list <- ggplotColours(n=22)
  
  arg <- match.call()
  
  if (label == TRUE) {
    ggplot(df,aes(x=factor(eval(arg$Chromosome)), y=logpval, fill=factor(eval(arg$Chromosome)))) + geom_violin() + geom_boxplot(width=0.1) + 
      geom_text_repel(data = subset(df, outlier == T), aes(label=eval(arg$Genename)),force_pull = 3, size=6, max.overlaps = 20) +
      scale_fill_manual(values = color_list, guide='none')
  } else {
    ggplot(df,aes(x=factor(eval(arg$Chromosome)), y=logpval, fill=factor(eval(arg$Chromosome)))) + geom_violin() + geom_boxplot(width=0.1) + 
      scale_fill_manual(values = color_list, guide='none')
  }
}




####### Merge Plink results ########

# You can merge many assoc.linear files, assoc files, etc.
# This function changes yout workding directory
# Set path as directory storing many assoc.linear files.
# It is convenient to set name of assoc.linear file as chrx.
# Prepare datalist with code below

#filelist = list.files(path = getwd(), pattern = "assoc.linear")
#datalist = lapply(filelist, function(x)read_table(x)) 
#names(datalist) <- str_remove_all(filelist, paste(c("chr16_",".assoc.linear"), collapse = "|")) #-- optional
#names(datalist) <- str_remove_all(filelist, paste(c("chr12.",".assoc.linear"), collapse = "|")) #--optional

mergelist <- function(dflist, 
                      div=c("BLOCK","GENE"), 
                      sep = "\\.") {
  for (i in 1:length(dflist)) {
    dflist[[i]]$X10 <- names(dflist[i])
    dflist[[i]] <- dflist[[i]] %>% separate(X10, div, sep = sep)
  }
  LDblockdf <- do.call(rbind, dflist)
  rownames(LDblockdf) <- NULL
  return(LDblockdf)
}




####### Merge Heatmap with dataframe made from mergelist function ########

Gene_Snp_Heatmap <- function(df, Pvalscale=c(0,1,2), snp_id_size=10, gene_symbol_size=20) {
  
  df <- df %>% mutate(LOGP = -log10(P))
  
  p1 <- ggplot(df, aes(x=fct_inorder(SNP), y=GENE, fill=LOGP, order=BP)) + geom_tile() + 
    viridis::scale_fill_viridis(values=scales::rescale(Pvalscale),
                                direction = 1, 
                                guide = guide_colorbar(barheight = 10, barwidth = 2)) + 
    xlab('SNP') + ylab('GENE') +
    facet_grid(~BLOCK, scales="free_x", space = "free_x") +
    theme(axis.text.x = element_text(angle = 90, size = snp_id_size),
          plot.title = element_text(size = 50, hjust = 0.5),
          axis.text.y = element_text(size = gene_symbol_size),
          axis.title = element_text(size=30),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(0.1, units = 'cm'),
          panel.background = element_rect(fill='white'),
          strip.background.x = element_rect(fill = c("grey")),
          strip.text = element_text(size = 20, lineheight = 0.2)) 
  
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



