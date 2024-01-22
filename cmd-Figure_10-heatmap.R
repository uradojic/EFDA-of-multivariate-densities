#Heatmap of clr coefficients of squared norms from the decomposition
#resulting from 20x20x20 grid data
#Karel Hron, 16.5.2023

#set working directory
setwd("...")
library(gplots)
library(robCompositions)

#abbreviations of districts
abbr=c("BN","BE","BK","BV","BM","BI","BR","CH","CV","CR","CL","CB","CK","DC","DO","FM","HB","A","HO","HK","JN","JE",
       "JC","JI","JH","KV","KI","KD","KT","KO","KM","KH","LI","LT","LN","ME","MB","MO","NA","NJ","NB","OC","OP",
       "OV","PU","PE","PI","PJ","PM","PS","PT","PY","PZ","PR","PB","PV","RA","RO","RK","SM","SO","ST","SY","SU",
       "TA","TC","TP","TR","TU","UH","UL","UO","VS","VY","ZL","ZN","ZR")

#preprocessing of the table
decomp=read.table("norm_decomp.txt") #data set with information compositions
cbind(decomp[,1],abbr)
rownames(decomp)=abbr
decomp=decomp[,-1]

################################################################
################################################################
#heatmap of plain clr coefficients

pdf("heatmap.pdf")
rgb.palette <- colorRampPalette(c("blue4","turquoise","white","orange","red4"),
                                                                   space = "rgb")
quant=quantile(c(as.matrix(cenLR(decomp[,-1])$x.clr)),probs = seq(0,1,0.01)) #color breaks according to quantiles
heatmap.2(as.matrix(cenLR(decomp[,-1])$x.clr),col = rgb.palette(100), key=TRUE, symkey=FALSE, breaks=quant, trace="none",cexCol=1,cexRow=0.4, scale = "none",
          Colv=FALSE,dendrogram="row")

dev.off()
