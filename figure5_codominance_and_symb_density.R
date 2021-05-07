library(ggplot2)
library(DESeq2)
library(tidyverse)
rm(list=ls())

setwd("~/Dropbox/project/zooxtype_sanity_check/coldata_remake/Evelyns_study")

#import coldata, counts, and zdat (zoox numbers)
coldata <- read.csv(file = "bleached_Coldata.csv")
load("/home/evelyn/Dropbox/project/zooxtype_sanity_check/coldata_remake/Evelyns_study/symbiont_type_counts.Rdata")
counts <- read.table(file = "host_counts.tsv", sep = "\t", header = T, row.names='Geneid')
counts = counts[,6:ncol(counts)]

#make it so the sample names of coldata, zdat, and counts match up, then colbind them
cleaned = sapply(colnames(counts), function(x) strsplit(x, "_")[[1]][1])
print(cleaned)
colnames(counts) = cleaned

#delete this study samples because they only have C symbionts
coldata <- coldata[!(coldata$my_title == "j1_thisStudy_PRJNA559404"),]
rownames(coldata) <- coldata$Run
keeps <- rownames(coldata)
keeps  
counts = counts[,colnames(counts) %in% keeps]

#make pool column
bar = coldata
bar = bar[(bar$my_title == "L_Barshis_bleachResillience_PRJNA177515"),]
bar = bar %>%
  separate(treatDescription, c("geno", "b","pool"), "_")
bar = bar %>% 
  select(Run,my_title,geno,pool,special,bleached)

pal = coldata
pal = pal[!(pal$my_title == "L_Barshis_bleachResillience_PRJNA177515"),]
pal = pal %>%
  separate(treatDescription, c("samp", "geno","pool","b","time"), "_")
pal = pal %>% 
  select(Run,my_title,geno,pool,special,bleached)
pal$pool <- gsub('400', 'MV', pal$pool)
pal$pool <- gsub('300', 'HV', pal$pool)

coldata = rbind(pal,bar)

#colbind zdat and coldata
dim(zdat)

#make zdat$Run into factors
zdat$Run <- as.factor(zdat$Run)
rownames(zdat) <- zdat$Run
keeps <- rownames(coldata)
keeps 
zdat = zdat[rownames(zdat) %in% keeps,]

coldata <- merge(coldata, zdat, by.x = "Run")

coldata$hostsym = ((coldata$cladeC +1)+(coldata$cladeD +1))/log((coldata$nonSym +1),10)
#coldata$hostsym = (coldata$nonSym)/((coldata$cladeC)+(coldata$cladeD))
coldata$logCDratio = log((coldata$cladeC +1)/(coldata$cladeD +1),10)
coldata <- coldata[order(coldata$logCDratio),]
#coldata = coldata[-c(1,2,3,4,5,6,7,8,9,10),]

#Turn your 'treatment' column into a character vector
rownames(coldata)=coldata$Run
coldata <- coldata[order(coldata$logCDratio),]
coldata$Run <- as.character(coldata$Run)
#Then turn it back into a factor with the levels in the correct order
coldata$Run <- factor(coldata$Run, levels=unique(coldata$Run)) 
#coldata$Run = rownames(coldata)

str(coldata)
g=ggplot(coldata,aes(Run,hostsym,color=logCDratio))+
  theme_classic()+
  geom_point()+
  #geom_smooth(color="grey30")+
  xlab("Samples from least to most Cladocopium")+
  #theme(axis.text.x = )
  ylab("Log ratio of host to symbiont counts")+
   theme(axis.text.x=element_blank(),
         axis.ticks.x=element_blank())  


png("~/Dropbox/project/bleww_final/Zoox_expr/fightzone/reviewer_concerns/figure5_res.png", width=2000, height=1200, res=400)
ggplot(coldata,aes(Run,hostsym,color=logCDratio,group=1))+
  theme_classic()+
  geom_point()+
  geom_smooth()+
  #geom_smooth(color="grey30")+
  xlab("Samples from least to most Cladocopium")+
  #theme(axis.text.x = )
  ylab("Log ratio of host to symbiont counts")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())  
plot(g)
dev.off()
