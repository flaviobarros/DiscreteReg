#' Fit regression models as Api=Xbeta for a table r x s from multinomial model ou product of multinomials
#' 
#'@param result Object with result from function estima.theta 
#'@param m.A First matrix that define the model 
#'@param m.X Second matrix that define the model
#'@return list with vpc, vbeta, epbeta e mcovbeta
#'@export
#'@import Matrix
#'@import plotrix

MR.Api.Xbeta<-function(result,m.A,m.X)
{
  label<-c(result$label)
  vpc<-cbind(c(result$vpc))
  mcov<-result$mcov
  mcov[mcov==0]<-0.000000000001
  mcovF <- m.A%*%mcov%*%t(m.A)
  imcovF <- solve(mcovF)
  vF <- m.A%*%vpc
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