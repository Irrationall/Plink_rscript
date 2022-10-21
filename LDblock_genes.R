library(dplyr)
library(tidyverse)
library(ggplot2)
library(viridis)

# Load association file between LD block SNPs and a gene


## Method 1 - NOGADA

ACSM1_block1 <- read_table("ACSM1_chr16_block1.assoc.linear")
ACSM1_block2 <- read_table("ACSM1_chr16_block2.assoc.linear")
DNAH3_block2 <- read_table("DNAH3_chr16_block2.assoc.linear")
THUMPD1_block1 <-  read_table("THUMPD1_chr16_block1.assoc.linear")
THUMPD1_block2 <-  read_table("THUMPD1_chr16_block2.assoc.linear")

dfs <- Filter(function(x) is(x, "data.frame"), mget(ls()))

mergelist(dfs)



## Method 2 - if assoc.linear files in one directory

filelist = list.files(path = getwd(), pattern = "assoc.linear")
datalist = lapply(filelist, function(x)read_table(x)) 
names(datalist) <- str_remove_all(filelist, paste(c("chr16_",".assoc.linear"), collapse = "|"))

mergelist(datalist)



## Function 1 - merge list to dataframe ##

mergelist <- function(dflist) {
  for (i in 1:length(dflist)) {
    dflist[[i]]$X10 <- names(dflist[i])
    dflist[[i]] <- dflist[[i]] %>% separate(X10, c("GENE","BLOCK"), sep = "_")
  }
  LDblockdf <- do.call(rbind, dflist)
  rownames(LDblockdf) <- NULL
  assign("LDblockdf", LDblockdf, envir = .GlobalEnv)
  return(LDblockdf)
}


## Function 2 - make a heatmap between genes and SNPs





