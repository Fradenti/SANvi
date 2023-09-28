#' Estimate Posterior Clustering Assignments
#'
#' This function estimates posterior clustering assignments based on posterior variational estimates
#' obtained from one of the model implemented in SANvi.
#' 
#' @param output an object of class \code{SANvb}, the output of one of the variational functions \code{\link{variational_CAM}}, 
#' \code{\link{variational_fiSAN}}, \code{\link{variational_fSAN}}.
#' @param ordered logical, if \code{TRUE} (default), the function sorts the distributional cluster labels reflecting the 
#' increasing values of medians of the data assigned to each DC.
#'
#' @return a list of class \code{clustering} containing 
#' \itemize{
#'   \item \code{obs_level}: a data frame containing the data values, their group indexes, the observational and distributional clustering assignments for each observation.
#'   \item \code{dis_level}: a vector with the distributional clustering assignment for each unit.
#' }
#' 
#' @export
#'
#' @seealso \code{\link{variational_CAM}}, \code{\link{variational_fiSAN}}, \code{\link{variational_fSAN}}, \code{\link{extract_best}}.
#'
#' 
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' y <- c(rnorm(100),rnorm(100,-5),rnorm(100,5),rnorm(100),
#'        rnorm(100),rnorm(100,-5),rnorm(100,5),rnorm(100))
#' g <- rep(1:4,rep(200,4))
#' 
#' # Fitting fiSAN via variational inference
#' est <- SANvi::variational_fiSAN(y,g,verbose = FALSE)
#' 
#' Estimate clustering assignments
#' estimate_clustering_vi(est)
#' }
estimate_clustering_vi  <- function(output, 
                                    ordered = TRUE){
  
  if(!inherits(output, "SANvb")){
    warning("The passed object should be of class 'SANvb'")  
  }
  
  obs_clust <- as.numeric(factor(unlist(lapply(output$sim$XI,
                                               function(x) apply(x, 1, which.max)))))
  dist_clust <- as.numeric(factor(apply(output$sim$RHO,1,which.max)))
  
  if(ordered){
    obs_clust <- unname(rank(tapply(output$params$y,
                                    INDEX = obs_clust,
                                    median))[obs_clust])
    
    dist_clust <- unname(rank(tapply(output$params$y,
                                     INDEX = rep(dist_clust,output$params$Nj),
                                     median))[dist_clust])
  }
  D = data.frame(Y  = output$params$y,
                 G  = output$params$group,
                 OC = obs_clust,
                 DC = rep(dist_clust,output$params$Nj))
  
  L <- (list(obs_level = D, 
             dis_level = dist_clust))
  structure(L,class = c("vi_clustering",class(L)))
}




#' @name estimate_clustering_vi
#'
#' @param x an object of class \code{variational_estimates}, which can be obtained from the
#' function \code{\link{estimate_atoms_weights_vi}}.
#' @param DC_num an integer or a vector of integers indicating which distributional clusters to plot.
#' @param type what type of plot should be drawn (only for the left-side plot). Possible types are "boxplot", "ecdf", and "scatter". 
#' @param palette_brewed (logical) the color palette to be used. Default is \code{R} base colors (\code{palette_brewed = FALSE}).
#' @param ... ignored.
#'
#' @importFrom graphics abline lines points boxplot
#' @importFrom stats median
#' @importFrom grDevices colorRampPalette
#' @importFrom scales alpha
#' @importFrom RColorBrewer brewer.pal
#' @export
#'
plot.vi_clustering <- function(x,
                               DC_num = NULL,
                               type = c("ecdf", "boxplot", "scatter"),
                               palette_brewed = FALSE,
                               ...) {
  
  type   <- match.arg(type)
  if (is.null(DC_num)) {
    DC_num <- 1:max(x$dis_level)
  }
  if(max(DC_num) > max(x$dis_level)){
    stop(paste0("There are less estimated DCs than requested.\n",
                "Please provide a number for DC_num between 1 and ", max(x$dis_level) ))
  }
  
  
  max_CD <- max(x$dis_level)
  if(palette_brewed){
    colpal <-
      rev(
        grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(max_CD))
  }else{
    colpal <- 1:max_CD
  }
  
  dix <- rank(tapply(x$obs_level$Y, x$obs_level$DC, stats::median))
  
  ind_ord_dis <- sapply(x$dis_level, function(g) dix[g])
  ind_ord_obs <- dix[x$obs_level$DC]
  
  
  inds_row <- which(ind_ord_obs %in% DC_num)
  inds_col <- which(ind_ord_dis %in% DC_num)
  
  suby <- (x$obs_level$Y[inds_row])
  subg <- (x$obs_level$G[inds_row])
  subDC <- x$obs_level$DC[inds_row]
  subOC <- x$obs_level$OC[inds_row]
  
  if(length(DC_num) == max(x$dis_level)){
    main_title <- paste0("All distributional clusters")
  }else if( length(DC_num) > 8 ){
    main_title <- NULL
  }else{
    main_title <- paste0("Distributional clusters # ", paste(DC_num, collapse = ", "))
  }
  
  
  
  if (type == "ecdf") {
    X <- unique(x$obs_level$G[inds_row])
    Nj <- table(x$obs_level$G)
    ysteps <- (0:(Nj[inds_col][1] - 1)) / 
      (Nj[inds_col][1] - 1)
    xsteps <- sort(suby[subg == X[1]])
    
    
    plot(
      ysteps ~ xsteps,
      type = "l",
      xlim = c(min(suby), max(suby)),
      ylim = c(0, 1),
      col = scales::alpha(colpal[ind_ord_dis[inds_col][1]], .5),
      xlab = "y",
      ylab = "eCDF",
      main = paste0("eCDFs colored by DC\n",main_title)
    )
    graphics::points(
      (ysteps ~ xsteps),
      cex = .1,
      col = scales::alpha(colpal[ind_ord_dis[inds_col][1]], .5)
    )
    graphics::abline(h = c(0, 1),
                     col = "gray",
                     lty = 3)
    
    
    for (j in 2:length(X)) {
      ysteps = (0:(Nj[X[j]] - 1)) / (Nj[X[j]] - 1)
      xsteps = sort(suby[subg == X[j]])
      
      graphics::lines((ysteps ~ xsteps),
                      col = scales::alpha(colpal[ind_ord_dis[inds_col][j]], .5))
      graphics::points(
        (ysteps ~ xsteps),
        cex = .1,
        col = scales::alpha(colpal[ind_ord_dis[inds_col][j]], .5)
      )
      
    }
    
  } else if (type=="boxplot"){
    
    graphics::boxplot(
      suby ~ subg,
      col = scales::alpha(colpal[ind_ord_dis[inds_col]], .7),
      pch = ".",
      main = paste0("Boxplots colored by DC\n",main_title),
      ylab = "y",
      xlab = "group"
    )
  }else{
    oldpar <- graphics::par()$mfrow
    
    graphics::par(mfrow=c(1,2))
    plot(suby ~ jitter(subg),
         pch=".",
         col=subDC,
         xlab = "Group index",
         ylab = "y",
         main = paste0("Observations colored by DC\n",main_title)
    )
    plot(suby ~ jitter(subg),
         pch=".",
         col=subOC,
         xlab = "Group index",
         ylab = "y",
         main = paste0("Observations colored by OC\n",main_title)
    )
    on.exit(graphics::par(mfrow=oldpar))
  }
}

#' @name estimate_clustering_vi
#'
#' @export
#'
print.vi_clustering <- function(x, ...){
  cat(paste("Number of estimated OCs:", length(unique(x$obs_level$OC))),"\n")
  cat(paste("Number of estimated DCs:",length(unique(x$dis_level))))
  invisible(x)
}
