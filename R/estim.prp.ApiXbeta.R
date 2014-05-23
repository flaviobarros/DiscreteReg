#' Estimate proportions under model Api=Xbeta
#' 
#'@param result Object with result from fucntion MR.Api.Xbeta
#'@param m.H Matrix that returns original proportions
#'@param v.F Vector that returns original proportions
#'@param gama Confidence level
#'@export
#'@import Matrix
#'@import plotrix

estim.prp.ApiXbeta<-function(result,m.H,v.F,gama)
  
{
  label<-result$label
  vpc <- result$vpc
  vbeta<-result$vbeta
  mcovbeta<-result$mcovbeta
  vpesm <- m.H%*%vbeta + v.F
  mcovsm <- m.H%*%mcovbeta%*%t(m.H)
  vepsm <- sqrt(diag(mcovsm))
  qic <- qnorm(0.5*(1+gama))
  LIIC <- c(vpesm)-c(qic*vepsm)
  LSIC <- c(vpesm)+ c(qic*vepsm)
  mIC  <- cbind(LIIC,LSIC)
  mIC[mIC[,1]<=0,1]=0
  mIC[mIC[,2]>=1,2]=1
  result<- list(vpesm=vpesm,vepsm=vepsm,mcovbeta=mcovbeta,mIC=mIC)
  plot(vpc,axes=F,ylim=c(min(vpc,mIC),max(vpc,mIC)),xlab="category",ylab="proportions",cex=1.2)
  plotCI(vpesm,ui=mIC[,2],li=mIC[,1],axes=F,pch=19,cex=1.2,cex.axis=1.2,cex.lab=1.2,add=T)
  axis(2,cex.axis=1.2) 
  axis(1,1:length(vpc),labels=label,cex.axis=1.2)
}