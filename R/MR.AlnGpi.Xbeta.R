#' Fit regression models as AlnGpi=Xbeta for a table r x s from a multinomial model or a product of multinomials
#' 
#'@param result Object with result form function estima.theta
#'@param m.A Matriz definidora do modelo de interesse
#'@param m.G Matriz definidora do modelo de interesse
#'@param m.X Matriz definidora do modelo de interesse
#'@return uma lista contendo vpc, vbeta, epbeta e mcovbeta
#'@export
#'@import Matrix
#'@import plotrix

MR.AlnGpi.Xbeta<-function(result,m.A,m.G,m.X)
{
  label<-c(result$label)
  vpc<-cbind(c(result$vpc))
  mcov<-result$mcov
  mcov[mcov==0]<-0.000000000001
  auxmG<-m.G%*%vpc
  auxmG[auxmG==0]<-0.000000000001
  m.B <- solve(diag(c(auxmG)))
  mPsi <- m.A%*%m.B%*%m.G
  mcovF <- mPsi%*%mcov%*%t(mPsi)
  auxav <- eigen(mcovF,only.values=TRUE)$values
  nlmcovF <- nrow(mcovF)
  while(min(auxav) <= 0.001)
  {
    mcovF <- mcovF + diag(0.00001,nlmcovF,nlmcovF)
    auxav <- eigen(mcovF,only.values=TRUE)$values
  }
  imcovF <- solve(mcovF)
  vF <- m.A%*%(log(auxmG))
  mcovbeta <- solve(t(m.X)%*%imcovF%*%m.X)
  vbeta<- mcovbeta%*%t(m.X)%*%imcovF%*%vF
  epbeta <- cbind(c(sqrt(diag(mcovbeta))))
  epl <- m.X%*%vbeta
  eQ <-t(vF-epl)%*%imcovF%*%(vF-epl)
  ngl<-nrow(m.X)-nrow(vbeta)
  epvalor<-1-pchisq(eQ,ngl)
  cat("Statistic from beta parameters","\n")
  print(round(cbind(vbeta,epbeta),2))
  cat("Test for goodness of fit =  ",round(eQ,2),"\n")
  cat("p-value = ",round(epvalor,4),"\n")
  cat("g.l. =",ngl,"\n")
  cat("Matrix A :","\n")
  print(m.A)
  cat("Matrix X :","\n")
  print(m.X)
  result<- list(vpc=vpc,vbeta=vbeta,epbeta=epbeta,mcovbeta=mcovbeta,label=label)
  return(result)
}