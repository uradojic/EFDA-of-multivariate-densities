#Anomaly detection for squared norms from decomposition of multivariate PDFs
#resulting from 20x20x20 grid data
#Karel Hron, 16.5.2023

#set working directory
setwd("...")
library(gplots)
library(robCompositions)

source("lr_ddc.R")

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

LR_DDC(decomp[,-1],pOutLR = 0.3,pOutRow = 0.75,showVals = NULL,pdfCellMap = TRUE)
