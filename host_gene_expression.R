#coral host gene expression

library(DESeq2)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)
theme_set(theme_cowplot())

#deseq================================================================================================

ll=load('~/Dropbox/project/bleww_final/scratchwork/host_fightzone.Rdata') 

#remove the dudes that you removed for the zoox
coldata <- coldata[order(coldata$logCDratio),]
coldata = coldata[-c(1,2,3,4,5,6,7,8,9,10),]
rownames(coldata)=coldata$Run
keeps = row.names(coldata)
counts = counts[,keeps]

cut=3
cc=counts
means=apply(cc,1,mean)
table(means>cut) #how many genes have mean below the cutoff?
counts=cc[means>cut,]

dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~c_t*minorLR+pool+logCcount)
dds=DESeq(dds)

resultsNames(dds)
rld = vst(dds)
rld.df=assay(rld)
colnames(rld.df) = coldata$Run
vsd = rld.df


#KOG lfc output
#for KOG, we do lfc
resultsNames(dds)

res = results(dds, name ="minorLR")
summary(res)
KminorLR = data.frame("gene"=rownames(res),
                      "lfc"=res$log2FoldChange)
save(KminorLR,file="~/Dropbox/project/bleww_final/scratchwork/KOG_minorLR1.Rdata")

res = results(dds, name ="c_th.minorLR")
summary(res)
Kixn = data.frame("gene"=rownames(res),
                  "lfc"=res$log2FoldChange)
save(Kixn,file="~/Dropbox/project/bleww_final/scratchwork/KOG_minorLR_ct_ixn1.Rdata")

res = results(dds, name ="c_t_h_vs_c")
summary(res)
ct = data.frame("gene"=rownames(res),
                "lfc"=res$log2FoldChange)
save(ct,file="~/Dropbox/project/bleww_final/scratchwork/KOG_ct1.Rdata")

#------------------------------------------------------------------------------
#KOG_MWU
#------------------------------------------------------------------------------

library(KOGMWU)

setwd("~/Dropbox/project/bleww_final/scratchwork")
data(adults.3dHeat.logFoldChange)
data(larvae.longTerm)
data(larvae.shortTerm)
data(gene2kog)

agene2kog = read.table("amil_gene2kogClass.tab",header=FALSE,sep = "\t")
ggene2kog = read.table("Amil.v2.eggnogWebsite.gene2kog.tsv",header = TRUE,sep = "\t")
load("~/Dropbox/project/bleww_final/scratchwork/KOG_minorLR_ct_ixn1.Rdata")
load("~/Dropbox/project/bleww_final/scratchwork/KOG_minorLR1.Rdata")
load("~/Dropbox/project/bleww_final/scratchwork/KOG_ct1.Rdata")
clustA = read.csv("clusterAstress_For_MWU.csv")
clustB = read.csv("clusterBstress_For_MWU.csv")
redmod = read.csv("MMred_goInput.csv")

#minorLR
minorLR = kog.mwu(KminorLR,agene2kog)
minorLR

#minorct
minorct = kog.mwu(ct,agene2kog)
minorct

#minorLR ixn
minorLR_ct_ixn = kog.mwu(Kixn, agene2kog)
minorLR_ct_ixn

#clustA stress
clustA = kog.mwu(clustA,ggene2kog)
clustA

#clustB stress
clustB = kog.mwu(clustB,ggene2kog)
clustB

#red module
redmod = kog.mwu(redmod,ggene2kog)
redmod

# Analyzing adult coral response to 3-day heat stress:
alfc.lth=kog.mwu(adults.3dHeat.logFoldChange,gene2kog) 
alfc.lth 

# coral larvae response to 5-day heat stress:
l.lth=kog.mwu(larvae.longTerm,gene2kog)
l.lth

# coral larvae response to 4-hour heat stress 
l.sth=kog.mwu(larvae.shortTerm,gene2kog)
l.sth

ktable=makeDeltaRanksTable(list("clustA"=clustA,"clustB"=clustB,"heated"=minorct,"minorLR"=minorLR,"minorLR-h interaction"=minorLR_ct_ixn))
ktable = na.omit(ktable)
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation") 

# exploring correlations between datasets
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor)
#?pairs
# p-values of these correlations in the upper panel:
pairs(ktable, lower.panel = panel.smooth, upper.panel = panel.cor.pval)

panel.lm=function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                   cex = 1, col.lmline = "red", ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    abline(lm(y[ok]~x[ok]), 
           col = col.lmline, ...)
}

pairs(ktable, lower.panel = panel.lm, upper.panel = panel.cor, cex.labels = 1.5)

# plotting individual delta-rank correlations:
#corrPlot(x="adults.long",y="larvae.long",ktable)
table = makeDeltaRanksTable(list()) #error

corrPlot(x="heated",y="clustA",ktable)

corrPlot(x="minorLR-h interaction",y="clustA",ktable)
corrPlot(x="heated",y="clustA",ktable)

corrPlot(x="minorLR",y="clustA",ktable)
corrPlot(x="minorLR",y="clustB",ktable)




