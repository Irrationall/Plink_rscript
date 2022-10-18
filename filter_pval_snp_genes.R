library(tidyverse)
library(dplyr)


# Before run this script, delete all .log filies in directory

filelist = list.files(path = "C:/Users/junes/OneDrive - 숭실대학교 - Soongsil University/3rd/2nd_semester/BS_exp/rs2123392")

setwd("C:/Users/junes/OneDrive - 숭실대학교 - Soongsil University/3rd/2nd_semester/BS_exp/rs2123392")

datalist = lapply(filelist, function(x)read_table(x)) 
names(datalist) <- filelist

df <- bind_rows(datalist)
df['GENE'] <- filelist
df <- column_to_rownames(df, "GENE")

sig <- subset(df, P < 0.05)

setwd('C:/Users/junes/OneDrive - 숭실대학교 - Soongsil University/3rd/2nd_semester/BS_exp/Results')

write.csv(df, "rs2123392.csv")


