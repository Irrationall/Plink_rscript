library(tidyverse)
library(dplyr)


# Before run this script, delete all .log filies in directory

filelist = list.files(path = "C:/Users/private/BS_exp/rs2123392")

setwd("C:/Users/private/BS_exp/rs2123392")

datalist = lapply(filelist, function(x)read_table(x)) 
names(datalist) <- filelist

df <- bind_rows(datalist)
df['GENE'] <- filelist
df <- column_to_rownames(df, "GENE")

sig <- subset(df, P < 0.05)

write.csv(df, "./Results/rs2123392.csv")


