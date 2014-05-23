#' Ajusta modelos de regressão do tipo Api=Xbeta para uma única tabela r x s oriunda de um modelo multinomial ou produto de multinomiais
#' 
#'@param result Objeto com o resultado da aplicação da função estima.theta 
#'@param m.B Primeira matriz definidora do modelo de interesse
#'@param m.X Segunda matriz definidora do modelo de interesse
#'@return uma lista contendo vpc, vbeta, epbeta e mcovbeta
#'@export

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
  cat("Est. dos par. beta","\n")
  print(round(cbind(vbeta,epbeta),2))
  cat("Teste para a qualidade do ajuste do modelo = ",round(eQ,2),"\n")
  cat("pvalor = ",round(epvalor,4),"\n")
  cat("g.l. =",ngl,"\n")
  cat("Matriz A :","\n")
  print(m.A)
  cat("Matriz X :","\n")
  print(m.X)
  result<- list(vpc=vpc,vbeta=vbeta,epbeta=epbeta,mcovbeta=mcovbeta,label=label)
  return(result)
}