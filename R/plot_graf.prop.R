#' Generate plots from estimated proportions (no models)
#' 
#'@param result Object with results from function estima.theta
#'@param oplas Labels orientation on x axis     
#'@param eixo Axis labels dimension    
#'@export
#'@import Matrix
#'@import plotrix      

plot_graf.prop <- function(result, oplas, eixo, ...)
{
  label <- c(result$label)
  vpc <- cbind(c(result$vpc))
  mIC <- (result$mIC)
  mIC <- cbind(c(mIC[,1]),c(mIC[,2]))
  
  plotCI(vpc,ui=mIC[,2],li=mIC[,1],
         axes=F,xlab="category",
         ylab="proportions",
         pch=19,cex=1.2,
         cex.axis=1.2,
         cex.lab=1.2, ...)
  axis(2,cex.axis=1.2)
  
  axis(1,1:length(vpc),labels=label,cex.axis=eixo,las=oplas)
}