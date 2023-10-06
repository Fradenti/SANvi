#' Estimate the Posterior Atoms and Weights of the Discrete Mixing Distributions
#'
#' This function estimates the posterior atoms and weights characterizing the discrete 
#' mixing distributions using the variational estimates obtained from one of the model implemented in SANvi.
#'
#'
#' @param output an object of class \code{SANvb}, which is the output of one of the variational functions \code{\link{variational_CAM}}, 
#' \code{\link{variational_fiSAN}}, \code{\link{variational_fSAN}}.
#'
#' @return an object of class \code{vi_atoms_weights}, which is matrix comprising posterior means, 
#' variances, and a columns for each estimated DC containing the posterior weights.
#' 
#' @export
#'
#' @seealso \code{\link{variational_CAM}}, \code{\link{variational_fiSAN}}, 
#' \code{\link{variational_fSAN}}, \code{\link{extract_best}}.
#'
#'
#' @examples
#' # Generate example data
#' set.seed(1232)
#' y <- c(rnorm(100),rnorm(100,5))
#' g <- rep(1:2,rep(100,2))
#' 
#' # Fitting fiSAN via variational inference
#' est <- SANvi:::variational_fiSAN(y,g,verbose = FALSE)
#' 
#' # Estimate posterior atoms and weights
#' estimate_atoms_weights_vi(est)
estimate_atoms_weights_vi <- function(output){
  
  if(class(output$sim)[1] == "CAMvi"){
    exp_weights = post_sb_weight(output$sim$a_bar_lk, output$sim$b_bar_lk)
    
  }else   if(class(output$sim)[1] == "fSANvi" | 
             class(output$sim)[1] == "fiSANvi"){
    exp_weights <- apply(output$sim$beta_bar_lk,2,function(x) x/sum(x))
  }else{
    stop("Provide valid variational inference output")
  }
  
  exp_sigma2 <- output$sim$theta_l[,4]/(output$sim$theta_l[,3]-1)
  exp_mu     <- output$sim$theta_l[,1]
  ind_row  <-  unique(unlist(lapply(output$sim$XI,
                                    function(x) apply(x, 1, which.max))))
  ind_col <- unique(apply(output$sim$RHO,1,which.max))
  
  
  W <- exp_weights[ind_row,ind_col]
  
  D <- data.frame(post_mean = exp_mu[ind_row],
                  post_var  = exp_sigma2[ind_row],
                  post_weigths_DC = W)
  
  
  if(is.null(nrow(W))){
    colnames(D)[-c(1:2)] <- paste0("post_weigths_DC",1:(ncol(D)-2))
  }else{
    wmeans <- apply(W,2,function(x) sum(x*D$post_mean))
    D <- D[order(D$post_mean),c(1:2,order(wmeans)+2)]
    colnames(D)[-c(1:2)] <- paste0("post_weigths_DC",1:(ncol(D)-2))
  }
  structure(D,class = c("vi_atoms_weights",class(D)))
  
}


#' @name estimate_atoms_weights_vi
#'
#' @param x an object of class \code{variational_estimates}, which can be 
#' obtained from the function \code{\link{estimate_atoms_weights_vi}}.
#' @param DC_num an integer or a vector of integers indicating which distributional clusters to plot.
#' @param lim optional value for \code{plot} method to adjust the limits of the x-axis. The atoms are plotted on a range
#' given by \code{min(posterior means)-2, max(posterior means)+2}. Default is set to 2.
#'
#' @importFrom graphics abline points lines
#' @importFrom stats dnorm
#'
#' @export
#'
plot.vi_atoms_weights <- function(x, 
                                  DC_num = NULL, 
                                  lim = 2, ...) {
  
  if (is.null(DC_num)) {
    DC_num <- 1:(ncol(x)-2)
  }
  if(max(DC_num) > (ncol(x)-2) ){
    stop(paste0("There are less estimated DCs than requested.\n",
                "Please provide a number for DC_num between 1 and ", (ncol(x)-2),"." ))
  }
  
  atoms   <-  x[, 1:2]
  seqq    <- seq(min(atoms[, 1]) - lim,
                 max(atoms[, 1]) + lim,
                 length.out = 250)

  width <- sqrt(atoms[, 2])
  
  weights <-  x[, DC_num[1] + 2]
  dens <- sapply(seqq, function(g) {
    sum(weights * stats::dnorm(g, atoms[, 1], sqrt(atoms[, 2])))
  })
  norm_dens <- dens / max(dens)
  
  if(length(DC_num) == ncol(x)-2){
    main_title <- paste0("All distributional clusters")
  }else if(length(DC_num) > 8){
    main_title <- paste0("Distributional clusters")
  }else{
    main_title <- paste0("Distributional clusters # ", paste(DC_num, collapse = ", "))
  }
  
  plot(
    weights ~ atoms[, 1],
    type = "h",
    lwd = 2,
    lend = 1,
    xlim = c(min(seqq), max(seqq)),
    ylim = c(0, 1),
    xlab = "y",
    ylab = "Posterior weights and normalized density",
    main = main_title
  )
  graphics::points(weights ~ atoms[, 1], cex = width)
  graphics::lines(norm_dens ~ seqq)
  
  if (length(DC_num) > 1) {
    for (j in 2:length(DC_num)) {
      weights <-  x[, DC_num[j] + 2]
      dens <- sapply(seqq, function(g) {
        sum(weights * stats::dnorm(g, atoms[, 1], sqrt(atoms[, 2])))
      })
      norm_dens <- dens / max(dens)
      
      graphics::points(
        weights ~ atoms[, 1],
        type = "h",
        lwd = 2,
        lend = 1,
        xlim = c(min(seqq), max(seqq)),
        ylim = c(0, 1),
        col = j
      )
      graphics::points(weights ~ atoms[, 1], cex = width, col = j)
      graphics::lines(norm_dens ~ seqq, col = j)
    }
    
  }
  
}


#' @name estimate_atoms_weights_vi
#'
#' @param thr argument for the \code{print()} method. It should be a small positive number, 
#' representing a threshold. If the posterior weight of a shared atom is below the threshold, the
#' atom is not reported.
#' @param ... ignored.
#' 
#' @export
#' 
print.vi_atoms_weights <- function(x, thr = 1e-2, ...){
  
  rows <- list()
  for(i in 1:(ncol(x)-2)){
    rows[[i]] <- which(x[,2+i] > thr)
  }
  n_g <- length(rows)
  
  cat(paste("Atoms with posterior weight >", thr, "\n"))
  cat("----------------------------------\n")
  cat(paste("Number of detected DCs:", n_g, "\n"))
  cat("----------------------------------\n")
  
  for(j in 1:n_g){
    
    cat(paste("\nDistributional cluster #",j,"\n"))
    
    if(length(rows[[j]])==0){
      cat(paste0("No atom has weight above the selected threshold of ",thr,"\n"))
      next()
    }
    Dsubj <- round(x[rows[[j]],c(1,2,j+2)],3)
    
    if(is.null(nrow(Dsubj))){
      Dsubj <- matrix(Dsubj,1,3)
      rownames(Dsubj) <- "1"
      colnames(Dsubj) <- c("post_mean", "post_var", "post_weight")
    }else{
      colnames(Dsubj)[3] = "post_weight"
    }
    print(data.frame(Dsubj))
    
  }
}


