#' Perform variational inference using multiple starting points.
#'
#' @description \code{variational_multistart} is the main function of the package. 
#' It is used to estimate via variational inference the three models we present (CAM, fiSAN, fSAN) while adopting multiple random starting points to better explore the variational parameter space.
#' The run that provides the highest Expected Lower BOund (ELBO) is usually the one considered for inference. Note that the arguments passed to this functions are a union of the arguments 
#' of the functions \code{variational_CAM}, \code{variational_fiSAN}, and \code{variational_fSAN}.
#'  
#' @usage 
#' variational_multistart(y, group, runs, cores = 1, 
#'                        model = c("fiSAN","CAM","fSAN"), 
#'                        maxL = 30, maxK = 20, 
#'                        m0 = 0, tau0 = .01, lambda0 = 3, gamma0 = 2, 
#'                        conc_par = NULL,  conc_hyperpar = c(1,1,1,1), 
#'                        alpha_bar = 0.05, beta_bar = 0.05, 
#'                        epsilon = 1e-6, root = 1234, maxSIM = 1e5, 
#'                        warmstart = TRUE)
#'
#' @param y vector of observations.
#' @param group vector of the same length of y indicating the group membership (numeric).
#' @param runs the number of multiple runs to launch.
#' @param cores the number of cores to dedicate to the multiple runs.
#' @param model a string specifying the model to use. It can be \code{"fiSAN"},\code{"CAM"}, or \code{"fSAN"}.
#' @param maxL,maxK integers, the upper bounds for the observational and distributional clusters to fit, respectively.
#' @param m0,tau0,lambda0,gamma0 hyperparameters on \eqn{(\mu, \sigma^2) \sim NIG(m_0, \tau_0, \lambda_0,\gamma_0)}.
#' @param conc_hyperpar,conc_par,alpha_bar,beta_bar these values crucially depend on the chosen model. See \code{\link{variational_CAM}}, \code{\link{variational_fiSAN}}, \code{\link{variational_fSAN}} for proper explanations.
#' @param epsilon the tolerance that drives the convergence criterion adopted as stopping rule.
#' @param root common part of the random seeds used to control the initialization in order to provide reproducibility even in paralleled settings.
#' @param maxSIM the maximum number of CAVI iteration to perform.
#' @param warmstart logical, if \code{TRUE}, the observational means of the cluster atoms are initialized with a k-means algorithm.
#'
#' @details For the details of the single models, see their specific documentations: \code{\link{variational_CAM}}, \code{\link{variational_fiSAN}}, and \code{\link{variational_fSAN}}.
#'
#'
#' @seealso \code{\link{variational_CAM}}, \code{\link{variational_fiSAN}}, \code{\link{variational_fSAN}}, \code{\link{extract_best}}.
#'
#' @return A list with all the runs performed. Each element of the list is a fitted variational model of class \code{SANvb}.
#' @export
#'
#' @examples
#' \donttest{
#' # Generate example data
#' set.seed(123)
#' y <- c(rnorm(100),rnorm(100,5))
#' g <- rep(1:2,rep(100,2))
#' 
#' # Estimate multiple models via variational inference
#' est <- variational_multistart(y, g, runs=5)
#' }
variational_multistart <- function(y,
                            group,
                            runs,
                            cores = 1,
                            model = c("fiSAN","CAM","fSAN"),
                            maxL = 30,
                            maxK = 20,
                            m0 = 0,
                            tau0 = .01,
                            lambda0 = 3,
                            gamma0 = 2,
                            conc_par = NULL,                            
                            conc_hyperpar = c(1,1,1,1),
                            alpha_bar = 0.05,
                            beta_bar = 0.05,
                            epsilon = 1e-6,
                            root = 1234,
                            maxSIM = 1e5,
                            warmstart = TRUE) {
  if (is.null(root)) {
    root <- sample(1:1e5, 1)
  }
  
  model <- match.arg(model)
  
  if(model == "fiSAN"){
    
    tmp_fiSAN <- function(i) {
      variational_fiSAN(
        y = y,
        group = group,
        maxL = maxL,
        maxK = maxK,
        m0 = m0,
        tau0 = tau0,
        lambda0 = lambda0,
        gamma0 = gamma0,
        conc_par = conc_par,
        conc_hyperpar = conc_hyperpar,
        beta_bar = beta_bar,
        epsilon = epsilon,
        ###########################
        seed = root * i,
        ###########################
        maxSIM = maxSIM,
        warmstart = warmstart
      )
    }
    
    if (cores == 1) {
      RESULTS <- list()
      for (i in 1:runs)
        RESULTS[[i]] <- tmp_fiSAN(i)
      
    } else{
      RESULTS <- parallel::mclapply(1:runs, function(i)
        tmp_fiSAN(i), mc.cores = cores)
      
    }
    

  }else if(model =="fSAN"){
    tmp_fSAN <- function(i) {
      variational_fSAN(
        y = y,
        group = group,
        maxL = maxL,
        maxK = maxK,
        m0 = m0,
        tau0 = tau0,
        lambda0 = lambda0,
        gamma0 = gamma0,
        alpha_bar = alpha_bar,
        beta_bar = beta_bar,
        epsilon = epsilon,
        ###########################
        seed = root * i,
        ###########################
        maxSIM = maxSIM,
        warmstart = warmstart
      )
      
  }
    if (cores == 1) {
      RESULTS <- list()
      for (i in 1:runs)
        RESULTS[[i]] <- tmp_fSAN(i)
      
    } else{
      RESULTS <- parallel::mclapply(1:runs, function(i)
        tmp_fSAN(i), mc.cores = cores)
      
    }
  
  
  }else if(model == "CAM"){
  
  tmp_cam <- function(i) {
    variational_CAM(
      y = y,
      group = group,
      maxL = maxL,
      maxK = maxK,
      m0 = m0,
      tau0 = tau0,
      lambda0 = lambda0,
      gamma0 = gamma0,
      conc_par = conc_par,
      conc_hyperpar = conc_hyperpar,
      epsilon = epsilon,
      ###########################
      seed = root * i,
      ###########################
      maxSIM = maxSIM,
      warmstart = warmstart
    )
  }
  
  if (cores == 1) {
    RESULTS <- list()
    for (i in 1:runs)
      RESULTS[[i]] <- tmp_cam(i)
    
  } else{
    RESULTS <- parallel::mclapply(1:runs, function(i)
      tmp_cam(i), mc.cores = cores)
    
  }

  }
  
  structure(RESULTS, class = c("multistart",class(RESULTS)))

}


#' @name variational_multistart
#'
#' @param x an object of class \code{multistart}, obtained from the \code{variational_multistart} function.
#' @param type a string specifying the type of plot. It can be either \code{"elbo"} or \code{"time"}. The former displays 
#' the elbo trajectories, highlighting the best run. The latter provides a summary of the computational times.
#' @param log_scale_iter logical. If \code{TRUE}, when plotting the elbo trajectories, the x-axes is displayed in log-scale, enhancing the visualization of the results.
#' @param ... further arguments passed to or from other methods.
#'
#' @importFrom graphics par hist text
#' @importFrom stats lm
#'
#' @export
#'
plot.multistart <- function(x,
                            type = c("elbo", "time"),
                            log_scale_iter = TRUE,
                            ...) {
  type <- match.arg(type)
  
  
  summa <- do.call(rbind,
                  lapply(x, function(x)
                    c(
                      min = min(x$sim$Elbo_val),
                      max = max(x$sim$Elbo_val),
                      iter = length(x$sim$Elbo_val),
                      secs = as.numeric(x$time, units =
                                          "secs")
                    )))
  
  ind_max <-  which.max(summa[, 2])
  
  
  
   if (type == "elbo") {
    if (log_scale_iter) {
      plot(
        1,
        type = "n",
        xlim = c(1, max(summa[, 3])),
        ylim = c(min(summa[, 1]), max(summa[, 2])),
        log = "x",
        main = paste("Results over", nrow(summa), "runs"),
        xlab = "Iterations - log scale",
        ylab = "ELBO", ...
      )
      
      for (i in 1:nrow(summa)) {
        graphics::points(x[[i]]$sim$Elbo_val ~ c(1:summa[i, 3]),
               type = "l",
               col = "lightgray")
      }
      graphics::points(
        x[[ind_max]]$sim$Elbo_val ~ c(1:summa[ind_max, 3]),
        type = "b",
        cex = .5,
        col = 4
      )
      graphics::text(max(summa[, 3]) * 1 / 4,
           min(summa[, 1]),
           paste("Best model index:", ind_max),
           family = "sans")
    } else{
      plot(
        1,
        type = "n",
        xlim = c(0, max(summa[, 3])),
        ylim = c(min(summa[, 1]), max(summa[, 2])),
        main = paste("Results over", nrow(summa), "runs"),
        xlab = "Iterations",
        ylab = "ELBO", ...
      )
      
      for (i in 1:nrow(summa)) {
        graphics::points(x[[i]]$sim$Elbo_val ~ c(1:summa[i, 3]),
               type = "l",
               col = "lightgray")
      }
      graphics::points(
        x[[ind_max]]$sim$Elbo_val ~ c(1:summa[ind_max, 3]),
        type = "b",
        cex = .5,
        col = 4
      )
      graphics::text(max(summa[, 3]) * 3 / 4,
           min(summa[, 1]),
           paste("Best model index:", ind_max),
           family = "sans")
      
    }
  } else if (type == "time") {
    oldpar <- graphics::par()$mfrow
    graphics::par(mfrow = c(2, 1))
    graphics::hist(summa[, 4],
         xlab = "Seconds",
         main = paste("Results over", nrow(summa), "runs"))
    plot(
      summa[, 4] ~ summa[, 3],
      xlab = "Iterations",
      ylab = "Seconds",
      main = paste("Results over", nrow(summa), "runs"),
      pch = ".",
      cex = 3, ...
    )
    graphics::abline(stats::lm(summa[, 4] ~ summa[, 3]), col = 2)
    on.exit(graphics::par(mfrow=oldpar))
  }
  
}



#' @name variational_multistart
#'
#' @param x an object of class \code{multistart}, obtained from the \code{variational_multistart} function.
#' @param ... ignored
#'
#' @export
#'
print.multistart <- function(x, ...) {
  summa <- do.call(rbind,
                  lapply(x,  function(q)
                    c(
                      min = min(q$sim$Elbo_val),
                      max = max(q$sim$Elbo_val),
                      iter = length(q$sim$Elbo_val),
                      secs = as.numeric(q$time, units =
                                          "secs")
                    )))
  cat(paste("Collection of", nrow(summa), "runs using multiple starting points\n"))
  cat(paste("Model fitted:",x[[1]]$model),"\n")
  cat(paste("Avg. elapsed time:", round(mean(summa[, 4]), 4), "seconds"))
}


#' Extract the best run from multiple trials
#'
#' @description A simple function to automatically extract the best run from a collection of fitted variational models.
#'
#' @param object an object of class \code{multistart}, obtained from the \code{variational_multistart} function.
#'
#' @return the single best run, an object of class \code{SANvb}.
#' @export
#'
#' @examples
#' \donttest{
#' # Generate example dataset
#' set.seed(123)
#' y <- c(rnorm(100),rnorm(100,5))
#' g <- rep(1:2,rep(100,2))
#'
#' # Estimate multiple models via variational inference
#' est <- SANvi:::variational_multistart(y, g, runs=5,
#'                                       alpha_bar = 3, beta_bar = 3,
#'                                       root=1234, warmstart = FALSE)
#' 
#' # Obtain best run
#' extract_best(est)
#' }
extract_best <- function(object) {
  
  if(!inherits(object, "multistart")){
    warning("The passed object should be of class 'multistart'")
  }
  
  summa <- do.call(rbind,
                   lapply(object,  function(q)
                     c(max = max(q$sim$Elbo_val))))
  ind <- which.max(summa)
  
  return(object[[ind]])
}
