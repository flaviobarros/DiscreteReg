#'Estimate parameters of a table r x s from a multinomial model or a product of multinomials
#' 
#' WARNING: Each function concatenates lines. So, in the case of a product of multinomials, 
#' each line must represent a multinomial
#'@param tabela A contingency table
#'@param modelo 1 - multinomial, 2 - multinomial product
#'@param gama Confidence level to build confidence intervals 
#'@return uma list with vpc, epc, mIC e mcov
#'@export
#'@import Matrix
#'@import plotrix

estima.theta <- function(tabela,modelo,gama)
{
  auxl1 <- c(rownames(tabela))
  auxl2 <- c(colnames(tabela))
  label <- c(paste(auxl1[1],"&",auxl2[1]))
  for(i in 1:nrow(tabela))
  {
    for(j in 1:ncol(tabela))
    {
      label<-rbind(label,c(paste(auxl1[i],"&",auxl2[j])))
    }
  }
  
  label<-cbind(label[2:nrow(label),])
  
  if(modelo == 1) # multinomial
  {
    vn <- c(t(tabela))
    n  <- sum(vn)
    ncat <- length(vn)
    vp <- vn/n
    #vp[vp==0]=0.01
    vpc <- cbind(vp)
    aux <- vpc%*%t(vpc)
    mcov <- matrix(as.numeric((as.matrix(Diagonal(ncat,vpc))- aux)/n),ncat,ncat)
    ep <- sqrt(diag(mcov))
    epc<-cbind(ep)
  } 
  
  else if(modelo == 2)# produto de multinomias
  {
    vn <- as.numeric(apply(tabela,1,sum))
    vpg<- tabela/vn#c(t(tabela/vn))
    ncatr<- ncol(tabela)
    nmult<-nrow(tabela)
    #vpg[vpg==0]=0.01
    vp <- vpg[1,]
    vpc <- cbind(vp)
    aux <-vpc%*%t(vpc)
    mcov <- matrix((as.matrix(Diagonal(ncatr,vpc))- aux)/vn[1],ncatr,ncatr)
    mcovg <-mcov
    
    for(j in 2:nmult)
    {
      vp <- vpg[j,]
      vpc <- cbind(vp)
      aux <-vpc%*%t(vpc)
      mcov <- matrix((as.matrix(Diagonal(ncatr,vpc))- aux)/vn[j],ncatr,ncatr)
      mcovg<-bdiag(mcovg,mcov)
    }
    
    mcov<-as.matrix(bdiag(mcovg))
    vp<- c(t(vpg))
    vpc<-cbind(vp)
    ep<- sqrt(diag(mcov))
    epc<-cbind(ep)
  }
  
  qic <- qnorm(0.5*(1+gama))
  LIIC <- c(vpc)-c(qic*epc)
  LSIC <- c(vpc)+ c(qic*epc)
  mIC  <- cbind(LIIC,LSIC)
  mIC[mIC[,1]<=0,1]=0
  mIC[mIC[,2]>=1,2]=1
  m.result <-cbind(round(vpc,2),round(epc,2),round(mIC,2))
  rownames(m.result)<-label
  colnames(m.result)<- c("Estimates","EP","LIIC","LSIC")
  cat("Estimates:","\n")
  print(m.result)
  result<- list(vpc=vpc,epc=epc,mIC=mIC,mcov=mcov,label=label)
  return(result)
}