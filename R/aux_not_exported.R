#' Compute posterior mean of stick-breaking weights.
#'
#' @param alk,blk the posterior variational parameters for the SB betas.
#' @noRd
#'
post_sb_weight <- function(alk, blk){
  
  K  <- ncol(alk)
  L  <- nrow(alk)
  # are we sure about this?
  logomega_post <- matrix(NA,L,K)
  
  for(k in 1:K){
    
    p2 = blk[,k]/(alk[,k]+blk[,k])
    p2 = c(1,p2[-L])
    
    logomega_post[,k] <- log(alk[,k]/(alk[,k]+blk[,k])) + cumsum(log(p2))
    
  }
  
  return(exp(logomega_post))
  
}

