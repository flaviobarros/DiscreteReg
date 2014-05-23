#' Realiza testes do tipo Bpi=D para uma única tabela r x s oriunda de um modelo multinomial ou produto de multinomiais
#' 
#'@param result Objeto com o resultado da aplicação da função estima.theta
#'@param m.B Matriz definidora da hipótese de interesse 
#'@param v.D Vetor definidor da hipótese de interesse
#'@export       

test.Bpi.D<-function(result,m.B,v.D)
{
  vpc<-cbind(c(result$vpc))
  mcov<-result$mcov
  mcov[mcov==0]<-0.000000000001
  e.Q <-t(m.B%*%vpc-v.D)%*%solve(m.B%*%mcov%*%t(m.B))%*%(m.B%*%vpc-v.D)
  ngl<-nrow(m.B)
  e.pvalor<-1-pchisq(e.Q,ngl)
  cat("Estatistica Q = ",round(e.Q,2),"\n")
  cat("pvalor = ",round(e.pvalor,4),"\n")
  cat("g.l. =",ngl,"\n")
  cat("Matriz B :","\n")
  print(m.B)
  cat("Vetor D :","\n")
  print(v.D)
}