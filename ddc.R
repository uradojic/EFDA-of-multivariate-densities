setwd("C:/Users/Ivana/Desktop/OneDrive - Univerzita Palackého v Olomouci/p h d/Karel/Densities_new2/Densities")

#######
####### functions + packages

library(robCompositions)
library(fda)
library(autoimage)
library(cellWise)
library(colorRamps)
library(gplots)


source("NulBaze.R")
source("SmoothingSpline0_verze2.R")

trapzc <- function(passo,y) 
{
  int<-passo*(0.5*y[1]+sum(y[2:(length(y)-1)]) +0.5*y[length(y)])
  return (int)
}

clr2density <- function(z, z_step, clr)
{
  if(is.fd(clr))
    return(exp(eval.fd(z,clr))/trapzc(z_step,exp(eval.fd(z,clr))))
  if(!is.fd(clr))
    return(exp(clr)/trapzc(z_step,exp(clr)))
}



######
###### load data

load("arithm_geom_densities_discrete_new.RData")

data=arithm_marginals_y # here data for Pb
range_PB=c(Pbmin,Pbmax)

n=nrow(data)
for(i in 1:ncol(data)){
  data[,i]=as.numeric(data[,i])
}

data.clr = cenLR((data))$x.clr #clr transformation from robComposition

par(mfrow=c(1,2))
matplot(gy,t(data),type="l",col="grey") # density
matplot(gy,t(data.clr),col=rainbow(n), type="l",pch=16)


t.fine=seq(min(gy),max(gy),length=1000) # grid
t.step=diff(t.fine)[1]

######
###### spline smoothing 

w = rep(1,ncol(data.clr))  #weights
k = 3                      #degree
der = 1                    #derivation
alfa = 0.90                
ch = 1                     # type of functional: (1-alfa)

knots=seq(Pbmin,Pbmax,length=length(gy)-2) 

t =gy 
f=as.numeric(data.clr[1,])


par(mfcol=c(2,1))
#results for one curve
Spline1 = SmoothingSpline0(knots=knots, t=t, f=f, w=w, k=k, der=der, alfa=alfa, ch=ch) #18
abline(v=knots,col="gray",lty=2)
Spline1[[1]] #functional J
Spline1[[2]] #vector of spline coefficients
Spline1[[3]] #1000*1; evaluation of the grid

J=c()
splines_coefs=matrix(nrow=dim(Spline1[[3]])[1],ncol=nrow(data.clr))
z_coefs=matrix(nrow=length(knots),ncol=nrow(data.clr))
for (i in 1:nrow(data.clr)){
  z_coefs[,i]=SmoothingSpline0(knots=knots, t=t, f=as.numeric(data.clr[i,]), w=w, k=k, der=der, alfa=alfa, ch=ch)[[2]]
  J[i]=SmoothingSpline0(knots=knots, t=t, f=as.numeric(data.clr[i,]), w=w, k=k, der=der, alfa=alfa, ch=ch)[[1]]
  splines_coefs[,i]=SmoothingSpline0(knots=knots, t=t, f=as.numeric(data.clr[i,]), w=w, k=k, der=der, alfa=alfa, ch=ch)[[3]]
}
sum(J) #total functional
splines_mean=apply(splines_coefs,1,mean) #  mean function

densities=clr2density(t.fine, t.step, splines_coefs[,i]) # 

###smoothed functional dataset + mean function 
options(scipen = 1)
par(mfcol=c(1,1))

matplot(t.fine,splines_coefs,main = "arithmetic_marginal_Pb; clr",xlab="",ylab="clr(density)",type="l",col="darkgrey")

lines(t.fine,splines_mean,xlab="",ylab="",type="l",col="darkblue",lwd=3)
abline(v=knots,col="lightpink",lty=2) 
abline(h=0,col="darkred")

######
###### Bspline + ZBspline basis representation

nknots=length(knots)
Ncoef=nknots+2-1
#n=N

Z = ZsplineBasis(knots = knots,k)$C0

data_clr = Z%*%(z_coefs)

options(scipen = 1)
# plot clr density
matplot(t.fine,data_clr, lty="solid",
        type="l",cex.lab=1.2,cex.axis=1.2,lwd=1,
        ylab="clr(density)",xlab="", col="darkgrey",
        main="Pb - arithm_marginal - clr")

abline(v=knots,col="lightpink",lty=2)
abline(h=0,col="darkred",lty=1)

# density in B2
data_b2 = NULL
for (i in 1:n){
  data_b2 = cbind(data_b2, clr2density(t.fine, t.step, data_clr[,i]))
}

matplot(t.fine,data_b2, 
        lty="solid", type="l",las=1,cex.lab=1.2,cex.axis=1.2,
        ylab="density",xlab="",col='darkgrey',lwd=1,
        main="Pb - arithm_marginal")
abline(v=knots,col="lightpink",lty=2)

###

# Expressing the compositional splines through the B-spline basis  
b_coef = t(ZsplineBasis(knots = knots,k)$D)%*%ZsplineBasis(knots = knots,k)$K%*%z_coefs
B = create.bspline.basis(range(knots), nbasis = dim(b_coef)[1],norder=k,breaks=knots)

Blog = create.bspline.basis(range(exp(knots)), nbasis = dim(b_coef)[1],norder=k,breaks=exp(knots))

plot(Blog,main="",log="x",col="darkblue",lty="solid",lwd=2,
     xlab=expression(paste('Pb concentration (mg kg'^-1,')')),cex.lab=1.2,cex.axis=1.2)

fd_data = fd(b_coef,B)

# #bbasis as a fd dataobject - check
 bbasis = fd(diag(nknots+1),B)
 matplot(t.fine,eval.fd(t.fine,bbasis),type="l",lty="solid")
 bbasis_eval=eval.fd(t.fine,bbasis)

######
###### DDC
######
 
X_b=t(b_coef)          # data matrix of b-spline coefficients
mu_b=apply(X_b,2,mean) # mean vector
Sigma_b=cov(X_b)       # variance structure
quant=0.95 

ddc_b=DDC(X_b, DDCpars = list()) #residuals

labels=abbr
cm_b=cellMap(ddc_b$stdResid,columnlabels=rep("",19),
             rowlabels=abbr,
             drawCircles=F,mTitle="Cu - arithmetic marginals",sizetitles = 0.5) #,columnlabels=col_labels,rowlabels = row_labels)

# for visualisation purposes the figure is split in 2 parts
# here: discrete representation

cm_b1=cellMap(ddc_b$stdResid,columnlabels=rep("",19),
              rowlabels=abbr,showrows=1:39,
              drawCircles=F,mTitle="Pb - arithmetic marginals",sizetitles = 0.5) #,columnlabels=col_labels,rowlabels = row_labels)

cm_b2=cellMap(ddc_b$stdResid,columnlabels=rep("",19),
              rowlabels=abbr,showrows=39:77,
              drawCircles=F,mTitle="Pb - arithmetic marginals",sizetitles = 0.5) #,columnlabels=col_labels,rowlabels = row_labels)
cm_b1
cm_b2

### not so sophisticated way to get the color scaling
gradient=matrix(cm_b$data$grad,nrow=n,ncol=Ncoef,byrow=FALSE)

basecolor=matrix(cm_b$data$CatNr,nrow=n,ncol=Ncoef,byrow=FALSE)
for (i in (1:n)){
  for ( j in (1:Ncoef)){
    basecolor[i,j]=ifelse(basecolor[i,j]==1,-1,basecolor[i,j])
    basecolor[i,j]=ifelse(basecolor[i,j]==2,1,basecolor[i,j])
    
  }
}
col_cm=basecolor*gradient 
data_colors_numbers=matrix(nrow=n,ncol=length(t.fine))

color_data=round(col_cm%*%t(bbasis_eval)*1000)

colors_fine=matlab.like(2001)
colors_blue=colors_fine[1000:1]
colors_red=colors_fine[1002:2001]
shadesOfGrey <- colorRampPalette(c("lightgrey", "black"))

colors_color_data=data_colors_numbers
for (i in 1:nrow(colors_color_data)){
  for(j in 1:ncol(colors_color_data)){
    if (color_data[i,j]>0){
      colors_color_data[i,j]=colors_red[color_data[i,j]]
      next
    }
    if (color_data[i,j]<0){
      colors_color_data[i,j]=colors_blue[abs(color_data[i,j])]
      next
    }
    else colors_color_data[i,j]=shadesOfGrey(n)[i]
  }
}

# visualisation: curves
par(mfcol=c(1,1))
matplot(t.fine,eval.fd(t.fine,fd_data),type="l",lwd=2,col=colors_fine[1001], xlab="x",ylab="clr density",
        lty="solid",main="Pb - arithmetic marginal")
abline(v=knots,col="darkgrey",lty="dashed")


for(i in 1:n){
  for (j in 1:1000){
    if(colors_color_data[i,j]!=shadesOfGrey(n)[i]){
      matpoints(t.fine[j],eval.fd(t.fine[j],fd_data[i]),
                col=colors_color_data[i,j],pch=20)
      
    }
  }
}




#visualisation: heatmap (FIGURE 7)
data=color_data
rownames(data)=abbr

rgb.palette <- colorRampPalette(matlab.like(2001)[seq(1,2001,by=200)],space = "rgb")
quant=seq(-1000,1000,length=1001)
heatmap.2(data,col = rgb.palette(1000), key=FALSE, symkey=FALSE, breaks=quant, trace="none",cexCol=0.6,cexRow=0.4, scale = "none",
          Colv=FALSE,dendrogram="row",main="Pb - arithm marginals - anomalies",labCol = round(exp(t.fine),3),
          xlab="concentration")


######
### heatmap 2 - functional values (FIGURE 6) 
data1=t(eval.fd(t.fine,fd_data))
rgb.palette <- colorRampPalette(c("blue4","turquoise","white","orange","red4"),
                                space = "rgb")
quant=quantile(c(data1),probs = seq(0,1,0.01)) #color breaks according to quantiles
heatmap.2(data1,col = rgb.palette(100), key=FALSE,keysize=0.8,density.info="none", symkey=FALSE, breaks=quant, trace="none",cexCol=0.8,cexRow=0.6, scale = "none",
          Colv=FALSE,dendrogram="row",main="Pb - arithmetical marginals",
          labCol=round(exp(t.fine),1),xlab="concentration",
          labRow=abbr)

