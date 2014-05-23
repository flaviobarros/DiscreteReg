#' Make tests like Bpi=D for a table r x s from multinomial model or product of multinomials
#' 
#'@param result Object with results from function estima.theta
#'@param m.B Matrix that define the model  
#'@param v.D Vector that define the hypothesis
#'@export 
#'@import Matrix
#'@import plotrix   

test.Bpi.D<-function(result,m.B,v.D)
{
  vpc<-cbind(c(result$vpc))
  mcov<-result$mcov
  mcov[mcov==0]<-0.000000000001
  e.Q <-t(m.B%*%vpc-v.D)%*%solve(m.B%*%mcov%*%t(m.B))%*%(m.B%*%vpc-v.D)
  ngl<-nrow(m.B)
  e.pvalor<-1-pchisq(e.Q,ngl)
  cat("Statistics Q = ",round(e.Q,2),"\n")
  cat("p-value = ",round(e.pvalor,4),"\n")
  cat("g.l. =",ngl,"\n")
  cat("Matrix B :","\n")
  print(m.B)
  cat("Vector D :","\n")
  print(v.D)
}