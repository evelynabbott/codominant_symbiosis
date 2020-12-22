#"FIGHTZONE" symbiont gene expression
#DESeq

#-----------------------------------------------------------------------------------------
#Make coldata
#-----------------------------------------------------------------------------------------
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

#delete study with only C symbionts
coldata <- coldata[!(coldata$my_title == "j1_thisStudy_PRJNA559404"),]
rownames(coldata) <- coldata$Run
keeps <- rownames(coldata)
keeps  
counts = counts[,colnames(counts) %in% keeps]

#colbind zdat and coldata
dim(zdat)
#make zdat$Run into factors
zdat$Run <- as.factor(zdat$Run)
rownames(zdat) <- zdat$Run
keeps <- rownames(coldata)
keeps 
zdat = zdat[rownames(zdat) %in% keeps,]

coldata <- merge(coldata, zdat, by.x = "Run")

#split up columns with "_" into their own columns
barcol <- subset(coldata, my_title == "L_Barshis_bleachResillience_PRJNA177515")
barcol <- separate(barcol, treatDescription, c("geno", "c_t", "pool"), sep = "_", remove = TRUE, convert = FALSE, extra = "warn", fill = "warn")

palcol <- subset(coldata, my_title == "k1_Palumbi_lab_heat_resilience_PRJNA274410")
palcol <- separate(palcol, treatDescription, c("null", "geno", "transplanted_pool", "c_t", "time"), sep = "_", remove = T, convert = F, extra = "warn", fill = "warn")
palcol$null = NULL
palcol$transplanted_pool = NULL

coldata <- dplyr::bind_rows(palcol, barcol)

#make logCDratio column for coldata
coldata$logCDratio = log((coldata$cladeC +1)/(coldata$cladeD +1),10)

coldata <- coldata %>%
  select(Run, my_title, geno, c_t, time, special, all, nonSym, cladeC, cladeD, logCDratio)

rownames(coldata) <- coldata$Run
#make all NAs zero
coldata$time[is.na(coldata$time)] <- 0
save(coldata, counts, file = "~/Dropbox/project/bleww_final/palbar_CvD/palbar_cvd_countscoldata.Rdata")

#cladocopium counts
load("/home/evelyn/Dropbox/project/zooxtype_sanity_check/coldata_remake/Evelyns_study/cladeC_counts.Rdata")
keeps <- rownames(coldata)
keeps  
counts = counts[,colnames(counts) %in% keeps]


#durusdinium counts
load("/home/evelyn/Dropbox/project/zooxtype_sanity_check/coldata_remake/Evelyns_study/cladeD_counts.Rdata")
keeps <- rownames(coldata)
keeps  
counts = counts[,colnames(counts) %in% keeps]


#-----------------------------------------------------------------------------------------
#DESeq2
#-----------------------------------------------------------------------------------------
#symbiont gene expression
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)
theme_set(theme_cowplot())

ll=load("~/Dropbox/project/bleww_final/scratchwork/fightzone_counts.Rdata") #C counts & metadata
#ll=load('~/Dropbox/project/bleww_final/scratchwork/fightzone_countsD.Rdata') #D counts & metadata

#make minorLR column
coldata$minorLR=coldata$logCDratio
coldata$minorLR[coldata$minorLR>0]=(-1)*coldata$minorLR[coldata$minorLR>0]
save(counts, coldata, file = "~/Dropbox/project/bleww_final/scratchwork/fightzone_countsC.Rdata") #cladocopium
#save(counts, coldata, file = "~/Dropbox/project/bleww_final/scratchwork/fightzone_countsD.Rdata") #durusdinium

#remove genes with low coverage
cut=3
cc=counts
means=apply(cc,1,mean)
table(means>cut)
counts=cc[means>cut,]

dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~c_t*minorLR+logCcount+time)
dds=DESeq(dds)

resultsNames(dds)
rld = vst(dds)
vsd=assay(rld)
colnames(vsd) = coldata$Run

#detect outliers
vsdt=as.data.frame(t(vsd))
sampleTree = hclust(dist(vsdt), method = "average");
par(cex = 0.6);
par(mar = c(0,5,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

#remove outliers
coldata <- coldata[order(coldata$logCDratio),]
coldata = coldata[-c(1,2,3,4,5,6,7,8,9,10),]
keeps = row.names(coldata)
counts = counts[,keeps]

#remove batch effect of study/time point and number of counts
library(limma)
vsd2=removeBatchEffect(vsd,batch=coldata$time,covariates=coldata$logCcount)

save(coldata,counts,vsd,vsd2,dds, file="~/Dropbox/project/bleww_final/scratchwork/fightzone.Rdata") #for C symbiont expression
#save(coldata,counts,vsd,vsd2,dds, file="~/Dropbox/project/bleww_final/scratchwork/fightzoneD.Rdata") #for D symbiont expression

#-----------------------------------------------------------------------------------------------
#WGCNA input
#-----------------------------------------------------------------------------------------------
datExpr=vsd2
datTraits=coldata
save(datExpr,datTraits,file="wgcna_inputC.Rdata") #do this twice, once for cladocopium and once for durusdinium 

#-----------------------------------------------------------------------------------------------
#GO_MWU input and volcano plot
#-----------------------------------------------------------------------------------------------
resultsNames(dds)

#volcano plot
res = results(dds, name ="minorLR")
summary(res)
res
sdf = data.frame(res) %>% 
  arrange(pvalue) %>% 
  mutate(sig = !is.na(padj) & padj < 0.1) %>% 
  as_tibble()
sdf
nsig = sum(sdf$sig)
tot = nrow(sdf)
g=ggplot(data=sdf, aes(x=log2FoldChange, y=-log(pvalue, 10), color=sig)) +
  geom_point(alpha=0.1) +
  scale_color_manual(values=c('black', 'red')) + 
  labs(x=bquote(log[2]~'fold difference'),
       y=bquote("-"*log[10]~'p-value'),
       subtitle='') +
  xlim(-2, 2)+
  guides(color=guide_legend(title="FDR<0.1"))
plot(g)

#save for GO_MWU
goout=data.frame("gene"=rownames(res),
                 "logP"=-log(res$pvalue,10))


sign=res$log2FoldChange>0
goout=goout %>% 
  mutate(logP=if_else(sign,
                      logP,
                      -logP)) %>%
  write_csv(path= "~/Dropbox/project/bleww_final/scratchwork/GO_fightzoneC.csv")


#-----------------------------------------------------------------------------------------------
#Resample counts to see if low count samples skew results
#-----------------------------------------------------------------------------------------------
zooxsums=apply(counts,2,sum)

dim(counts)

targetcount = min(zooxsums) #resample with lowest, 15k and 30k
#targetcount = 15000 
#targetcount = 30000

cladeC.resampled=c()
for (s in 1:ncol(counts)) {
  probs= counts[,s]/zooxsums[s]
  cts=hist(sample(c(1:nrow(counts)), targetcount,prob=probs,replace=TRUE),breaks=c(0:nrow(counts)),plot=F)$counts
  cladeC.resampled=data.frame(cbind(cladeC.resampled, cts))
}
#reassign the colnames back to counts
colnames(cladeC.resampled) = colnames(counts)
colnames(cladeC.resampled)
rownames(counts)
rownames(cladeC.resampled) = rownames(counts)

str(cladeC.resampled)
str(counts)

zooxsums=apply(cladeC.resampled,2,sum)

counts <- cladeC.resampled
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~c_t*minorLR+logCcount+time)
dds=DESeq(dds)

#get variance stabilized counts and save them
#rld = vst(dds)
rld = varianceStabilizingTransformation(dds)
rld.df=assay(rld)
colnames(rld.df) = coldata$Run
vsd = rld.df

library(limma)
covars=cbind(model.matrix(~0+coldata$time),logc=coldata$logCcount)[,-1]
vsd2=removeBatchEffect(vsd,batch=coldata$time,covariates=coldata$logCcount)


save(coldata,counts,vsd,vsd2,dds, file="~/Dropbox/project/bleww_final/scratchwork/resampled_fightzoneC.Rdata")
#make GO_MWU input the same way as ln 138-171









