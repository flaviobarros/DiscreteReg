#' Gera gráficos das porporções estimadas (sem considerar modelos)
#' 
#'@param result Objeto com o resultado da aplicação da função estima.theta
#'@param oplas Orientação dos "labels" do eixo "x"      
#'@param eixo Dimensão dos labels dos eixos     
#'@export       

plot.graf.prop <- function(result, oplas, eixo)
{
  label <- c(result$label)
  vpc <- cbind(c(result$vpc))
  mIC <- (result$mIC)
  mIC <- cbind(c(mIC[,1]),c(mIC[,2]))
  plotCI(vpc,ui=mIC[,2],li=mIC[,1],axes=F,xlab="categoria",ylab="proporções",pch=19,cex=1.2,cex.axis=1.2,cex.lab=1.2)
  axis(2,cex.axis=1.2) 
  axis(1,1:length(vpc),labels=label,cex.axis=eixo,las=oplas)
}