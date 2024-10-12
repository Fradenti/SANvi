#' Plotting the variational inference output
#' 
#' @description Plot method for objects of class \code{SANvb}. 
#' The function displays two graphs, meant to analyze the estimated distributional and observational clusters.
#' 
#' @param x object of class \code{SANvb} (the result of a call to \code{\link{variational_CAM}}, 
#' \code{\link{variational_fiSAN}}, \code{\link{variational_fSAN}}.
#' 
#' @param ... additional graphical parameters to be passed.
#' 
#' @return The function plots the path followed by the ELBO.
#' 
#' @seealso \code{\link{print.SANvb}}, \code{\link{variational_CAM}}, 
#' \code{\link{variational_fiSAN}}, \code{\link{variational_fSAN}}.
#'
#' @export
plot.SANvb <- function(x, ...){
  
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old.par))    
  
  graphics::par(mfrow= c(1,1))
  plot(x$sim$Elbo_val,
       xlab = "Iterations - log scale", ylab = "ELBO",
       main = paste(x$model, "- ELBO"),type="b",cex=.5,
       log="x", ...)
  
}





