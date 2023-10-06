#' Mean Field Variational Bayes estimation of fiSAN 
#'
#' @description \code{variational_fiSAN} is used to perform posterior inference under the finite-infinite shared atoms nested (fiSAN) model with Gaussian likelihood. 
#' The model uses a Dirichlet process mixture prior at the distributional level, and finite Dirichlet mixture at the observational one.
#' 
#' @usage 
#' variational_fiSAN(y, group,
#'                   maxL = 30, maxK = 20,
#'                   m0 = 0, tau0 = .01, lambda0 = 3, gamma0 = 2,
#'                   conc_hyperpar = c(1,1), conc_par = NULL,
#'                   beta_bar = .005,
#'                   epsilon = 1e-6, seed = NULL,
#'                   maxSIM = 1e5, warmstart = TRUE, verbose = FALSE)
#'                   
#' @param y Numerical vector of observations (required).
#' @param group Numerical vector of the same length of \code{y}, indicating the group membership (required).
#' @param maxL,maxK integers, the upper bounds for the observational and distributional clusters to fit, respectively
#' @param m0,tau0,lambda0,gamma0 Hyperparameters on \eqn{(\mu, \sigma^2) \sim NIG(m_0, \tau_0, \lambda_0,\gamma_0)}.
#' @param conc_hyperpar,conc_par Vectors of values used the concentration parameters of  of the stick-breaking representation for the distributional and observational DPs, respectively. 
#' The following two arguments can be passed. Specifically,
#' \describe{
#'   \item{\code{conc_hyperpar}}{a vector with 2 positive entries: \eqn{(s_1^\alpha,s_2^\alpha)}. 
#'   If a random concentration parameter \eqn{\alpha} is adopted, the specification is
#'  \eqn{\alpha \sim Gamma(s_1^\alpha,s_2^\alpha)}. Default set to unitary vector.}
#'   \item{\code{conc_par}}{a vector with one positive entry \eqn{\alpha}. Default is set to \code{NULL}. If specified, the previous argument is ignored
#'   and the two concentration parameters are assumed fixed and equal to \code{alpha}.}
#'  }
#' @param beta_bar the hyperparameter of the symmetric observational Dirichlet distribution.
#' @param epsilon the tolerance that drives the convergence criterion adopted as stopping rule
#' @param seed random seed to control the initialization.
#' @param maxSIM the maximum number of CAVI iteration to perform.
#' @param warmstart logical, if \code{TRUE}, the observational means of the cluster atoms are initialized with a k-means algorithm.
#' @param verbose logical, if \code{TRUE} the iterations are printed.
#'
#'
#' @details 
#' \strong{Data structure}
#' 
#' The finite-infinite common atoms mixture model is used to perform inference in nested settings, where the data are organized into \eqn{J} groups. 
#' The data should be continuous observations \eqn{(Y_1,\dots,Y_J)}, where each \eqn{Y_j = (y_{1,j},\dots,y_{n_j,j})} 
#' contains the \eqn{n_j} observations from group \eqn{j}, for \eqn{j=1,\dots,J}. 
#' The function takes as input the data as a numeric vector \code{y} in this concatenated form. Hence \code{y} should be a vector of length
#' \eqn{n_1+\dots+n_J}. The \code{group} parameter is a numeric vector of the same size as \code{y} indicating the group membership for each
#' individual observation. 
#' Notice that with this specification the observations in the same group need not be contiguous as long as the correspondence between the variables
#' \code{y} and \code{group} is maintained.
#'
#' \strong{Model}
#' 
#' The data are modeled using a Gaussian likelihood, where both the mean and the variance are observational-cluster-specific, i.e., 
#' \deqn{y_{i,j}\mid M_{i,j} = l \sim N(\mu_l,\sigma^2_l)}
#' where \eqn{M_{i,j} \in \{1,\dots,L \}} is the observational cluster indicator of observation \eqn{i} in group \eqn{j}.
#' The prior on the model parameters is a Normal-Inverse-Gamma distribution \eqn{(\mu_l,\sigma^2_l)\sim NIG (m_0,\tau_0,\lambda_0,\gamma_0)}, 
#' i.e., \eqn{\mu_l\mid\sigma^2_l \sim N(m_0, \sigma^2_l / \tau_0)}, \eqn{1/\sigma^2_l \sim Gamma(\lambda_0, \gamma_0)} (shape, rate).
#'
#' \strong{Clustering}
#' 
#' The model performs a clustering of both observations and groups. 
#' The clustering of groups (distributional clustering) is provided by the allocation variables \eqn{S_j \in \{1,2,\dots\}}, with 
#' \deqn{Pr(S_j = k \mid \dots ) = \pi_k  \qquad \text{for } \: k = 1,2,\dots}
#' The distribution of the probabilities is \eqn{ \{\pi_k\}_{k=1}^{\infty} \sim GEM(\alpha)}, 
#' where GEM is the Griffiths-Engen-McCloskey distribution of parameter \eqn{\alpha}, 
#' which characterizes the stick-breaking construction of the DP (Sethuraman, 1994). 
#' 
#' The clustering of observations (observational clustering) is provided by the allocation variables \eqn{M_{i,j} \in \{1,\dots,L\}}, with
#' \deqn{ Pr(M_{i,j} = l \mid S_j = k, \dots ) = \omega_{l,k} \qquad \text{for } \: k = 1,2,\dots \, ; \: l = 1,\dots,L. }
#' The distribution of the probabilities is \eqn{(\omega_{1,k},\dots,\omega_{L,k})\sim Dirichlet_L(\beta/L,\dots,\beta/L)} for all \eqn{k = 1,2,\dots}. 
#' Here, the dimension \eqn{L} is fixed.
#'
#'
#'
#' @return \code{variational_fiSAN} returns a list of class \code{SANvb} containing four objects:
#' \itemize{
#'   \item \code{model}: name of the fitted model.
#'   \item \code{params}: list containing the data and the parameters used in the simulation. Details below.
#'   \item \code{sim}: list containing the simulated values (optimized variational parameters). Details below.
#'   \item \code{time}: total computation time.
#' }
#'
#'
#' \strong{Data and parameters}:
#' \code{params} is a list with the following components:
#' \describe{
#' \item{\code{y, group, Nj, J}}{Data, group labels, group frequencies, and number of groups.}
#' \item{\code{K, L}}{Number of fitted distributional and observational clusters.}
#' \item{\code{m0, tau0, lambda0, gamma0}}{Model hyperparameters.}
#' \item{\code{epsilon, seed}}{The threshold controlling the convergence criterion and the random seed adopted to replicate the run.}
#' \item{(\code{hyp_alpha1,hyp_alpha2}) or \code{alpha}}{Hyperparameters on \eqn{\alpha} (if \eqn{\alpha} random);
#'   or provided value for \eqn{\alpha} (if fixed).}
#' \item{\code{beta_bar}}{the hyperparameter governing all the finite Dirichlet distributions at the observational level.}
#' }
#'
#'
#' \strong{Simulated values}:
#' \code{sim} is a list with the following components:
#' \describe{
#' \item{\code{theta_l}}{Matrix of size (L,4).
#'    Each row is a posterior variational estimate of the four normal-inverse gamma hyperparameters.}
#' \item{\code{Elbo_val}}{Vector containing the values of the ELBO.}
#' \item{\code{XI}}{A list of length J. Each element is a matrix of size (N, L)
#'    posterior variational probability of assignment of assignment of the i-th observation in the j-th group to the l-th OC, 
#'    i.e., \eqn{\hat{\xi}_{i,j,l} = \hat{\mathbb{Q}}(M_{i,j}=l)}.}
#' \item{\code{RHO}}{Matrix of size (J, K).
#'    Each row is a posterior variational probability of assignment of the j-th group to the k-th DC, i.e., \eqn{\hat{\rho}_{j,k} = \hat{\mathbb{Q}}(S_j=k)}. }
#' \item{\code{a_tilde_k,b_tilde_k}}{Vector of updated variational parameters of the Beta distributions governing the distributional stick-breaking process.}
#' \item{\code{beta_bar_lk}}{Matrix of updated variational parameters of the Dirichlet distributions governing the observational clustering (arranged by column).}
#' \item{\code{conc_hyper}}{If the concentration parameters is chosen to be random, these object contain a vector with the two updated hyperparameters.}
#' \item{\code{alpha}}{If the concentration parameters is chosen to be fixed, this object contains the passed values.}
#'}
#'
#'@export
#'
#'
#'@examples
#' set.seed(1234)
#' y <- c( rnorm(100) ,rnorm(100,5))
#' g <- rep( 1:2, rep(100,2))
#' est <- SANvi:::variational_fiSAN( y, g, verbose = FALSE,epsilon = 1e-2)
#'
variational_fiSAN <- function(y, 
                    group,
                    maxL = 30,
                    maxK = 20,
                    m0 = 0,
                    tau0 = .01,
                    lambda0 = 3,
                    gamma0 = 2,
                    conc_hyperpar = c(1,1),
                    conc_par = NULL,
                    beta_bar = .005,
                    epsilon = 1e-6,
                    seed = NULL,
                    maxSIM = 1e5,
                    warmstart = TRUE,
                    verbose = FALSE){


  L <- maxL
  K <- maxK

  if(length(y) != length(group)){
    stop("The number of observations and groups must match")
  }

  Nj <- tapply(y,group, length)
  J  <- max(group)

  if(is.null(seed)){seed <- round(stats::runif(1,1,10000))}
  # random init
  set.seed(seed)


  params <- list(y = y,
                group = group,
                Nj = Nj,
                J = J,
                K = K,
                L = L,
                m0 = m0, tau0 = tau0,
                lambda0 = lambda0, gamma0 = gamma0,
                beta_bar = beta_bar,
                epsilon = epsilon,
                seed = seed)



  if(!is.null(conc_par)) {
    if( !(all(conc_par>0)) ){
      stop("Not all the specified concentration parameters are positive.")
    }
    if( length(conc_par) < 1 ){
      stop("Please provide one value for the concentration parameter.")
    }
    # params$alpha1 <- 1 # a tilde
    params$alpha <- conc_par[1] # b tilde
    random_conc <- FALSE
  }else if(!is.null(conc_hyperpar) ) {
    if( !(all(conc_hyperpar>0)) ){
      stop("Not all the specified concentration hyperparameters are positive.")
    }
    if( length(conc_hyperpar) < 2 ){
      stop("Please provide two values for the concentration hyper-parameters.")
    }
    params$hyp_alpha1 <- conc_hyperpar[1]
    params$hyp_alpha2 <- conc_hyperpar[2]
    random_conc <- TRUE
  }else{
    stop("Please provide a consistent specification for the concentration parameters")
  }


  Y_grouped <- list()
  for(j in 1:J){
    Y_grouped[[j]] <- y[group==j] # this is a list, each element is a vector with observations
  }



  # random init
  if(warmstart){
    ml  <- stats::kmeans(unlist(Y_grouped),centers = L, algorithm="MacQueen",
                         iter.max = 100)$centers
  }else{
    ml  <- stats::runif(L, min(unlist(Y_grouped)),max(unlist(Y_grouped)))
  }


  kl       <- stats::rgamma(L,1,10)
  lambdal  <- stats::rgamma(L,1,1)
  gammal   <- stats::rgamma(L,1,1)



  ###############################################################################################
  XI_ijl = list()
  for(j in 1:J){
    log.XI_il  <- array(stats::rbeta( Nj[j] * L, 1, 1),dim = c( Nj[j], L))
    Z           <- apply(log.XI_il, c(1), function(x) matrixStats::logSumExp(x))
    XI_ijl[[j]] <- exp(sapply(1:L, function(qq) log.XI_il[,qq]-Z,simplify = "array"))
  }
  ###############################################################################################
  if(K<J){
    log.RHO_jk <- matrix(stats::rbeta(J*K,1,1),J,K)
    Z2         <- apply(log.RHO_jk, c(1), function(x) matrixStats::logSumExp(x))
    RHO_jk     <- exp(sapply(1:K, function(qq) log.RHO_jk[,qq]-Z2,simplify = "matrix"))
  }else if(K==J){
    RHO_jk     <- diag(J)
  }else{
    RHO_jk     <- cbind(diag(J), matrix(0,J,K-J))
  }

  if(random_conc){
    
    start <- Sys.time()
    results <- main_vb_fiSAN_CP_cpp(
      Y_grouped = Y_grouped,
      XI_ijl = XI_ijl,
      L = L,
      K = K,
      J = J,
      RHO_jk = RHO_jk,
      Nj = Nj,
      m0 = m0,
      k0 = tau0,
      a0 = lambda0,
      b0 = gamma0,
      ml = ml,
      kl = kl,
      al = lambdal,
      bl = gammal,
      conc_hyper = conc_hyperpar,
      beta_bar = beta_bar,
      epsilon = epsilon,
      maxSIM = maxSIM,
      verbose = verbose
    )
    end <- Sys.time()

  }else{

    start <- Sys.time()
    results <- main_vb_fiSAN_cpp(
      Y_grouped = Y_grouped,
      XI_ijl = XI_ijl,
      L = L,
      K = K,
      J = J,
      RHO_jk = RHO_jk,
      Nj = Nj,
      m0 = m0,
      k0 = tau0,
      a0 = lambda0,
      b0 = gamma0,
      ml = ml,
      kl = kl,
      al = lambdal,
      bl = gammal,
      a_tilde = 1,
      b_tilde = params$alpha,
      beta_bar = beta_bar,
      epsilon = epsilon,
      maxSIM = maxSIM,
      verbose = verbose
    )
    end <- Sys.time()

  }
  output <- list("model" = "fiSAN",
                 "params" = params,
                 "sim"= results,
                 "time" = end - start)
  
  class(output$sim) <- "fiSANvi"
  structure(output, class = c("SANvb",class(output)))
  
}

