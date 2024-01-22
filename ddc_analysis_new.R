  setwd("C:/Users/Ivana/Desktop/OneDrive - Univerzita Palackého v Olomouci/p h d/Karel/Densities_new2/Densities")
  
  load("arithm_geom_densities_discrete_new.RData")
  library(robCompositions)
  library(fda)
  library(autoimage)
  library(cellWise)
  library(colorRamps)
  
  #data=arithm_marginals_x
  data=arithm_marginals_y
  #data=arithm_marginals_z
  #data=geom_marginals_x
 #data=geom_marginals_y
  #data=geom_marginals_z
  range_PB=c(Pbmin,Pbmax)
  #head(data)
  #View(data)
  
  #visualisation
  matplot(gy,t(data),type="l",col="grey") #density
  
  ### ### ### ### ### ###
  
  
  
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
  
  ### ### ### ### ### ###
  n=nrow(data)
  class(data)
  for(i in 1:ncol(data)){
    data[,i]=as.numeric(data[,i])
  }
  
  data.clr = cenLR((data))$x.clr #VYSTUP
  #data.clr = data #VYSTUP
  apply(data.clr,1,sum)
  matplot(gy,t(data.clr),col=rainbow(n), type="l",pch=16)
  t.fine=seq(min(gy),max(gy),length=1000)
  t.step=diff(t.fine)[1]
  ### ###
  ###### uzly, vyhlazovac? splajny ######
  
  w = rep(1,ncol(data.clr)) 
  k = 3
  der = 1
  alfa = 0.90
  ch = 1     # funkcional s (1-alfa)
  
  knots=seq(Pbmin,Pbmax,length=length(gy)-2) # v?b?r podle funkcion?lu J
  #knots=c(seq(log(0.08),log(50),length=5),seq(log(50),log(300),length=8),seq(log(300),log(2000),length=3))# v?b?r podle funkcion?lu J
  
  t =gy #102x1 #st?edy interval?
  #f=as.numeric(zrnitost.clr[1,]) #1x102
  f=as.numeric(data.clr[1,])
  
  
  par(mfcol=c(2,1))
  ###### funkce SmoothingSpline0 ######
  #generov?n? koeficient? (z) pro splajnov? funkce [[1]] + y sou?adnice [[3]]
  #v?stup pro jedinou k?ivku, tj. jedin? vzorek
  Spline1 = SmoothingSpline0(knots=knots, t=t, f=f, w=w, k=k, der=der, alfa=alfa, ch=ch) #18
  abline(v=knots,col="gray",lty=2)
  Spline1[[1]] #funkcional J
  Spline1[[2]]
  Spline1[[3]]#1000*1; dimenze pro vyps?n? n??e
  
  hust=clr2density(exp(t.fine), t.step, Spline1[[3]])
  
  
  J=c()
  spliny_koef=matrix(nrow=dim(Spline1[[3]])[1],ncol=nrow(data.clr))
  z_koef=matrix(nrow=length(knots),ncol=nrow(data.clr))
  for (i in 1:nrow(data.clr)){
    z_koef[,i]=SmoothingSpline0(knots=knots, t=t, f=as.numeric(data.clr[i,]), w=w, k=k, der=der, alfa=alfa, ch=ch)[[2]]
    J[i]=SmoothingSpline0(knots=knots, t=t, f=as.numeric(data.clr[i,]), w=w, k=k, der=der, alfa=alfa, ch=ch)[[1]]
    spliny_koef[,i]=SmoothingSpline0(knots=knots, t=t, f=as.numeric(data.clr[i,]), w=w, k=k, der=der, alfa=alfa, ch=ch)[[3]]
  }
  sum(J) #puvodni: 2178.385
  #text(locator(1), names(response)[z], col = "darkblue")
  
  spliny_prumer=apply(spliny_koef,1,mean) # hodnoty "pr?m?rn? k?ivky - ze v?ech splajnov?ch k?ivek
  
  hust=clr2density(t.fine, t.step, spliny_koef[,i])
  
  ###### grafick? v?stup II ######
  ###vyhlazen? k?ivky v jednom grafu + "pr?m?rn?" k?ivka
  options(scipen = 1)
  par(mfcol=c(1,2))
  
  #matplot(t.fine,spliny_koef,main = "arithmetic_marginal_Cu; smoothed clr",xlab="",ylab="clr(density)",type="l",col="darkgrey")
  matplot(t.fine,spliny_koef,main = "arithmetic_marginal_Pb; smoothed clr",xlab="",ylab="clr(density)",type="l",col="darkgrey")
  #matplot(t.fine,spliny_koef,main = "arithmetic_marginal_Zn; smoothed clr",xlab="",ylab="clr(density)",type="l",col="darkgrey")
  #matplot(t.fine,spliny_koef,main = "geometric_marginal_Cu; smoothed clr",xlab="",ylab="clr(density)",type="l",col="darkgrey")
  #matplot(t.fine,spliny_koef,main = "geometric_marginal_Pb; smoothed clr",xlab="",ylab="clr(density)",type="l",col="darkgrey")
  #matplot(t.fine,spliny_koef,main = "geometric_marginal_Zn; smoothed clr",xlab="",ylab="clr(density)",type="l",col="darkgrey")
  
  #lines(t.fine,spliny_prumer,xlab="",ylab="",type="l",col="darkblue",lwd=3)
  #abline(v=exp(knots),col="lightgray",lty=2) #po nagenerov?n? uzl?
  abline(h=0,col="darkred")
  
  #matplot(t.fine,t(hustoty),log="x",lty=1:length(stredy), type="l",xlab = expression (paste("Particle size (",mu,"m)")),ylab="density",col="darkgrey")
  #lines(exp(t.fine),hust,xlab="",ylab="",type="l",col="darkblue",lwd=3)
  #abline(v=exp(knots),col="lightgray",lty=2) #po nagenerov?n? uzl?
  
  dev.off()
  
  #abline(h=0,col="darkred")
  #dev.off() 
  
  ###
  ###
  ###
  ###########
  nknots=length(knots)
  Ncoef=nknots+2
  n=77
  Ncoef=19
  
  Z = ZsplineBasis(knots = knots,k)$C0
  
  data_clr = Z%*%(z_koef) #1000x77
  
  options(scipen = 1)
  # plot clr density
  matplot(t.fine,data_clr, lty="solid",
          type="l",cex.lab=1.2,cex.axis=1.2,lwd=1,
          ylab="clr(density)",xlab="", col="darkgrey",
          #        main="Cu - arithm_marginal - clr")
          main="Pb - arithm_marginal - clr")
          #main="Zn - arithm_marginal - clr")
          #main="Cu - geom_marginal - clr")
          #main="Pb - geom_marginal - clr")
          #main="Zn - geom_marginal - clr")
  
  abline(v=knots,col="gray",lty=2)
  abline(h=0,col="red",lty=1)
  
  t.fine2=seq(range(t.fine)[1],range(t.fine)[2],length=10000)
  # density in B2
  data_b2 = NULL
  for (i in 1:n){
    data_b2 = cbind(data_b2, clr2density(t.fine, t.step, data_clr[,i]))
  }
  
  matplot(t.fine,data_b2, 
          lty="solid", type="l",las=1,cex.lab=1.2,cex.axis=1.2,
          ylab="density",xlab="",col='darkgrey',lwd=1,
          main="Pb - arithm_marginal")
  abline(v=knots,col="gray",lty=2)
  
  ###
  
  # Expressing the compositional splines through the B-spline basis  
  b_coef = t(ZsplineBasis(knots = knots,k)$D)%*%ZsplineBasis(knots = knots,k)$K%*%z_koef
  B = create.bspline.basis(range(knots), nbasis = dim(b_coef)[1],norder=k,breaks=knots)
  
  plot(B,main="",col="darkblue",lty="solid",lwd=2,
       xlab=expression(paste('Pb concentration (mg kg'^-1,')')),cex.lab=1.2,cex.axis=1.2)
  
  Blog = create.bspline.basis(range(exp(knots)), nbasis = dim(b_coef)[1],norder=k,breaks=exp(knots))
  
  plot(Blog,main="",log="x",col="darkblue",lty="solid",lwd=2,
       xlab=expression(paste('Pb concentration (mg kg'^-1,')')),cex.lab=1.2,cex.axis=1.2)
  
 range(t.fine)
  
  fd_data = fd(b_coef,B)
  #save.image("Cu_arithm_I.RData")
  #plot(fd_data_men)
  
  #bbasis as a fd dataobject
  bbasis = fd(diag(nknots+1),B)
  #plot(bbasis)
  #matplot(t.fine,eval.fd(t.fine,bbasis),type="l",lty="solid")
  
  bbasis_eval=eval.fd(t.fine,bbasis)
  apply(bbasis_eval,1,sum) #YES!!!
  
  #DDC
  X_b=t(b_coef)
  mu_b=apply(X_b,2,mean)
  Sigma_b=cov(X_b)
  quant=0.95 #0.99
  
  ddc_b=DDC(X_b, DDCpars = list()) #pro vypo?et rezidui
  
  labels=abbr
  #my_labels=function(x,y){
  #  v=c()
  #  for(i in 1:y){
  #    v=c(v,paste0(x,i+1949))
  #  }
  #  print(v)
  #}
  
  #lab_men=my_labels("M",72)
  
  # cm_b=cellMap(X_b,ddc_b$stdResid,columnlabels=rep("",19),
  #              rowlabels=abbr,
  #              drawCircles=F,mTitle="Cu - arithmetic marginals",sizetitles = 0.5) #,columnlabels=col_labels,rowlabels = row_labels)
  # 
  ### PRIDANO
   cm_b=cellMap(ddc_b$stdResid,columnlabels=rep("",19),
                rowlabels=abbr,
                drawCircles=F,mTitle="Cu - arithmetic marginals",sizetitles = 0.5) #,columnlabels=col_labels,rowlabels = row_labels)
   ### PRIDANO
  cm_b$scales
  
  #write.table(ddc_b$stdResid,"Zn_geom_residuals.csv")
  #write.table(t(b_coef),"Zn_geom_b_coeffs.csv")
  
              
  #PRIDANO
   cm_b1=cellMap(ddc_b$stdResid,columnlabels=rep("",19),
                 rowlabels=abbr,showrows=1:39,
                 drawCircles=F,mTitle="Cu - arithmetic marginals",sizetitles = 0.5) #,columnlabels=col_labels,rowlabels = row_labels)
   
   cm_b2=cellMap(ddc_b$stdResid,columnlabels=rep("",19),
                 rowlabels=abbr,showrows=39:77,
                 drawCircles=F,mTitle="Cu - arithmetic marginals",sizetitles = 0.5) #,columnlabels=col_labels,rowlabels = row_labels)
   #PRIDANO
  
  # cm_b1=cellMap(X_b,ddc_b$stdResid,columnlabels=rep("",19),
  #               rowlabels=abbr,showrows=1:39,
  #               drawCircles=F,mTitle="Cu - arithmetic marginals",sizetitles = 0.5) #,columnlabels=col_labels,rowlabels = row_labels)
  # 
  # cm_b2=cellMap(X_b,ddc_b$stdResid,columnlabels=rep("",19),
  #               rowlabels=abbr,showrows=39:77,
  #               drawCircles=F,mTitle="Cu - arithmetic marginals",sizetitles = 0.5) #,columnlabels=col_labels,rowlabels = row_labels)
  # 
  X11()
  cm_b1
  X11()
  cm_b2
  #par(cex.lab=0.5)
  #X11()
  #plot.new()
  #plot(cm_b,add=T)
  #legend(locator(1),legen=abbr[1:40],horiz=F,cex=0.5,bty="n")
  #legend(locator(1),legen=abbr[41:77],horiz=F,cex=0.5,bty="n")
  
  gradient=matrix(cm_b$data$grad,nrow=n,ncol=Ncoef,byrow=FALSE)
  
  basecolor=matrix(cm_b$data$CatNr,nrow=n,ncol=Ncoef,byrow=FALSE)
  for (i in (1:n)){
    for ( j in (1:Ncoef)){
      basecolor[i,j]=ifelse(basecolor[i,j]==1,-1,basecolor[i,j])
      basecolor[i,j]=ifelse(basecolor[i,j]==2,1,basecolor[i,j])
      
    }
  }
  col_cm=basecolor*gradient #72*20
  #bbasis_eval
  
  data_colors_numbers=matrix(nrow=n,ncol=length(t.fine))#96x1000
  #
  
  pokus1=round(col_cm%*%t(bbasis_eval)*1000)
  
  colors_fine=matlab.like(2001)
  colors_blue=colors_fine[1000:1]
  colors_red=colors_fine[1002:2001]
  shadesOfGrey <- colorRampPalette(c("lightgrey", "black"))
  
  colors_pokus1=data_colors_numbers
  for (i in 1:nrow(colors_pokus1)){
    for(j in 1:ncol(colors_pokus1)){
      if (pokus1[i,j]>0){
        colors_pokus1[i,j]=colors_red[pokus1[i,j]]
        next
      }
      if (pokus1[i,j]<0){
        colors_pokus1[i,j]=colors_blue[abs(pokus1[i,j])]
        next
      }
      else colors_pokus1[i,j]=shadesOfGrey(n)[i]
    }
  }
  X11()
  par(mfcol=c(1,1))
  #matplot(t.fine,eval.fd(t.fine,fd_data),type="l",lwd=1,col='grey',xlab="x",ylab="clr density",
  #lty="solid")
  #matplot(t.fine,eval.fd(t.fine,fd_data_men),type="l",lwd=2,col=shadesOfGrey(n), xlab="x",ylab="clr density",
  #        lty="solid")
  matplot(t.fine,eval.fd(t.fine,fd_data),type="l",lwd=2,col=colors_fine[1001], xlab="x",ylab="clr density",
                 lty="solid",main="Cu - arithmetic marginal")
          #lty="solid",main="Pb - arithmetic marginal")
          #lty="solid",main="Zn - arithmetic marginal")
          #lty="solid",main="Cu - geometric marginal")
          #lty="solid",main="Pb - geometric marginal")
          #lty="solid",main="Zn - geometric marginal")
  
  abline(v=knots,col="darkgrey",lty="dashed")
  
  
  for(i in 1:n){
    for (j in 1:1000){
      if(colors_pokus1[i,j]!=shadesOfGrey(n)[i]){
        matpoints(t.fine[j],eval.fd(t.fine[j],fd_data[i]),
                  col=colors_pokus1[i,j],pch=20)
        
      }
    }
  }
  
  #save.image("results_aritm_marginals_Cu.RData")
  #save.image("results_aritm_marginals_Pb.RData")
  #save.image("results_aritm_marginals_Zn.RData")
  #save.image("results_geom_marginals_Cu.RData")
  #save.image("results_geom_marginals_Pb.RData")
  #save.image("results_geom_marginals_Zn.RData")
  
  ################################################################################################
  ################################################################################################
  ################################################################################################
  setwd("C:/Users/Ivana/Desktop/OneDrive - Univerzita Palackého v Olomouci/p h d/Karel/Densities_new2/Densities")
  
  load("arithm_geom_densities_discrete_new.RData")
 # Cu_fine=seq(Cmin,Cmax,length=1000)
  Pb_fine=seq(Pbmin,Pbmax,length=1000)
#  Zn_fine=seq(Znmin,Znmax,length=1000)
  
#Cu_marginals  
  load("results_geom_marginals_Cu.RData")
  library(robCompositions)
  library(fda)
  library(autoimage)
  library(cellWise)
  library(colorRamps)
  library(gplots)
  #Cu_arithm_eval=t(eval.fd(t.fine,fd_data))
  data=pokus1
  rownames(data)=abbr
  
  
  #heatmap 
  pdf("PDF_heatmap_ddc_Cu_geom_new_ax.pdf")
    rgb.palette <- colorRampPalette(matlab.like(2001)[seq(1,2001,by=200)],space = "rgb")
  #quant=quantile(c(data),probs = seq(0,1,length=1000)) #color breaks according to quantiles
  quant=seq(-1000,1000,length=1001)
  heatmap.2(data,col = rgb.palette(1000), key=TRUE, symkey=FALSE, breaks=quant, trace="none",cexCol=0.6,cexRow=0.4, scale = "none",
            Colv=FALSE,dendrogram="row",main="Cu - geom marginals - anomalies",labCol = round(exp(Cu_fine),3),
            xlab="concentration")
  
  dev.off()

  
#Pb_marginals  
  load("results_geom_marginals_Pb.RData")

  data=pokus1
  rownames(data)=abbr
  
  
  #heatmap 
  
  pdf("PDF_heatmap_ddc_Pb_geom_new_ax.pdf")
  rgb.palette <- colorRampPalette(matlab.like(2001)[seq(1,2001,by=200)],space = "rgb")
  #quant=quantile(c(data),probs = seq(0,1,length=1000)) #color breaks according to quantiles
  quant=seq(-1000,1000,length=1001)
  heatmap.2(data,col = rgb.palette(1000), key=TRUE, symkey=FALSE, breaks=quant, trace="none",cexCol=0.6,cexRow=0.4, scale = "none",
            Colv=FALSE,dendrogram="row",main="Pb - geom marginals - anomalies",labCol = round(exp(Pb_fine),3),
            xlab="concentration")
  
  dev.off()
  
#Zn_marginals  
  load("results_geom_marginals_Zn.RData")

  data=pokus1
  rownames(data)=abbr
  
  
  #heatmap 
  
  pdf("PDF_heatmap_ddc_Zn_geom_new_ax.pdf")
  rgb.palette <- colorRampPalette(matlab.like(2001)[seq(1,2001,by=200)],space = "rgb")
  #quant=quantile(c(data),probs = seq(0,1,length=1000)) #color breaks according to quantiles
  quant=seq(-1000,1000,length=1001)
  heatmap.2(data,col = rgb.palette(1000), key=TRUE, symkey=FALSE, breaks=quant, trace="none",cexCol=0.6,cexRow=0.4, scale = "none",
            Colv=FALSE,dendrogram="row",main="Zn - geom marginals - anomalies",labCol = round(exp(Zn_fine),3),
            xlab="concentration")
  
  dev.off()
  
#Cu_arithm  
  load("results_aritm_marginals_Cu.RData")

  data=pokus1
  rownames(data)=abbr
  
  
  #heatmap 
  
  pdf("PDF_heatmap_ddc_Cu_arithm_new_ax.pdf")
  rgb.palette <- colorRampPalette(matlab.like(2001)[seq(1,2001,by=200)],space = "rgb")
  #quant=quantile(c(data),probs = seq(0,1,length=1000)) #color breaks according to quantiles
  quant=seq(-1000,1000,length=1001)
  heatmap.2(data,col = rgb.palette(1000), key=TRUE, symkey=FALSE, breaks=quant, trace="none",cexCol=0.6,cexRow=0.4, scale = "none",
            Colv=FALSE,dendrogram="row",main="Cu - arithmetic marginals - anomalies",labCol = round(exp(Cu_fine),3),
            xlab="concentration")
  
  dev.off()
  
  #Pb_arithm  
  load("results_aritm_marginals_Pb.RData")
  
  data=pokus1
  rownames(data)=abbr
  
  
  #heatmap 
  
  #pdf("PDF_heatmap_ddc_Pb_arithm_new_ax.pdf")
  rgb.palette <- colorRampPalette(matlab.like(2001)[seq(1,2001,by=200)],space = "rgb")
  #quant=quantile(c(data),probs = seq(0,1,length=1000)) #color breaks according to quantiles
  quant=seq(-1000,1000,length=1001)
  heatmap.2(data,col = rgb.palette(1000), key=TRUE, symkey=FALSE, breaks=quant, trace="none",cexCol=0.6,cexRow=0.4, scale = "none",
            Colv=FALSE,dendrogram="row",main="Pb - arithm marginals - anomalies",labCol = round(exp(Pb_fine),3),
            xlab="concentration")
  
  dev.off()

  labels10=rep("",1000)
  for(i in 1:10){
    labels10[(i-1)*100+1]=round(exp(Pb_fine),1)[(i-1)*100+1]
  }
  labels10[1000]=round(exp(Pb_fine),1)[1000]
  
              #new
  #pdf("PDF_heatmap_ddc_Pb_arithm_new_ax.pdf")
  rgb.palette <- colorRampPalette(matlab.like(2001)[seq(1,2001,by=200)],space = "rgb")
  #quant=quantile(c(data),probs = seq(0,1,length=1000)) #color breaks according to quantiles
  quant=seq(-1,1,length=1001)
  X11()
  heatmap.2(data/1000,col = rgb.palette(1000), key=TRUE,keysize=0.8,density.info="none", symkey=FALSE, breaks=quant, trace="none",cexCol=0.8,cexRow=0.6, scale = "none",
            Colv=FALSE,dendrogram="row",main="Pb - arithmetic marginals - anomalies",
            labCol = round(exp(Pb_fine),1),
                                                                                  #labCol = labels10,
            xlab="concentration")
  
  #Zn_arithm  
  load("results_aritm_marginals_Zn.RData")
  
  data=pokus1
  rownames(data)=abbr
  
  
  #heatmap 
  
  pdf("PDF_heatmap_ddc_Zn_arithm_new_ax.pdf")
  rgb.palette <- colorRampPalette(matlab.like(2001)[seq(1,2001,by=200)],space = "rgb")
  #quant=quantile(c(data),probs = seq(0,1,length=1000)) #color breaks according to quantiles
  quant=seq(-1000,1000,length=1001)
  heatmap.2(data,col = rgb.palette(1000), key=TRUE, symkey=FALSE, breaks=quant, trace="none",cexCol=0.6,cexRow=0.4, scale = "none",
            Colv=FALSE,dendrogram="row",main="Zn - arithm marginals - anomalies",labCol = round(exp(Zn_fine),3),
            xlab="concentration")
  
  dev.off()
  
  
  ###
  
  rgb.palette2 <- colorRampPalette(matlab.like(2001)[seq(1,2001,by=200)],space = "rgb")
  #quant=quantile(c(data),probs = seq(0,1,length=1000)) #color breaks according to quantiles
  quant=seq(-1000,1000,length=1001)
  heatmap.2(data,col = rgb.palette(1000), key=TRUE, symkey=FALSE, breaks=quant, trace="none",cexCol=0.6,cexRow=0.4, scale = "none",
            Colv=FALSE,dendrogram="row",main="Zn - arithm marginals - anomalies",labCol = round(exp(Zn_fine),3),
            xlab="concentration")
  
  dev.off()