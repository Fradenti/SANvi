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
#' @seealso \code{\link{print.SANvb}}, \code{\link{variational_CAM}}, 
#' \code{\link{variational_fiSAN}}, \code{\link{variational_fSAN}}.
#'
#' @export
plot.SANvb <- function(x, ...){
  
  oldpar <- graphics::par()$mfrow
  
  graphics::par(mfrow= c(1,2))
  plot(x$sim$Elbo_val,
       xlab = "Iterations - log scale", ylab = "ELBO",
       main = paste(x$model, "- ELBO"),type="b",cex=.5,
       log="x", ...)
  cl_col <- as.numeric(factor(apply(x$sim$RHO,1,which.max)))
  cl_row <- unique(unlist(lapply(x$sim$XI,
                                 function(x) apply(x, 1, which.max))))
  
  boxplot(x$params$y~x$params$group,col=scales::alpha(cl_col, .5), xlab="Group",ylab="y",
          main = paste(x$model, "- DC and posterior means"))
  abline(h = x$sim$theta_l[unique(cl_row),1],col=4,lty=2)
  on.exit(graphics::par(mfrow = oldpar))
  
}





