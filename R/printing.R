#' Print variational inference output
#' 
#' @description Print method for objects of class \code{SANvb}.
#'
#' @param x object of class \code{SANvb} (the result of a call to \code{\link{variational_CAM}}, 
#' \code{\link{variational_fiSAN}}, \code{\link{variational_fSAN}}.
#' @param ... further arguments passed to or from other methods.
#' 
#' @export
print.SANvb <- function(x, ... ){
  cat(paste("Variational inference results for", x$model ,"\n"))
  cat("----------------------------------------------\n")
  cat(paste("L:",x$params$L,"- K:",x$params$K,"\n"))
  cat(paste("Threshold:",x$params$epsilon,"\n"))
  cat(paste("ELBO value:", round(max(x$sim$Elbo_val),3),"\n"))
  cat(paste("Convergence reached in",length(x$sim$Elbo_val),"iterations\n"))
  cat(paste("Elapsed time:",round(as.numeric(x$time),3),units(x$time),"\n"))
  invisible(x)
}




