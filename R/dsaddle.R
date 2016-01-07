#####
#' Empirical saddlepoint density 
#' @description Gives a pointwise evaluation of the empirical saddlepoint and optionally of
#'              its gradient at position y
#'
#' @param y points at which the SPA is evaluated (d dimensional vector) or an n by d matrix, each row indicating
#'          a different position.
#' @param X n by d matrix containing the data.
#' @param tol Tolerance used to assess the convergence of the rootfinding routine used to fit
#'            the saddlepoint density. Default value is 1e-6.
#' @param decay Rate at which the SPA falls back on a normal density. Should be a positive number,
#'              by default set to 0.5.
#' @param deriv If TRUE also the gradient of the log-saddlepoint density is returned.
#' @param log If TRUE the log of the saddlepoint density is returned.
#' @param control A list of control parameters with entries:
#'         \itemize{
#'         \item{ \code{method} }{The method used to calculate the normalizing constant. 
#'                                Either "LAP" (laplace) or "IS" (importance sampling).}
#'         \item{ \code{nNorm} }{If control$method == "IS", this is the number of importance samples used.}
#'         \item{ \code{tol} }{The tolerance used to assess the convergence of the solution to the saddlepoint equation.
#'                             The default is 1e-6.}
#'         }
#' @return A list with entries:
#'         \itemize{
#'         \item{ \code{llk} }{The value of the empirical saddlepoint at y;}
#'         \item{ \code{mix} }{The mixture of normal-saddlepoint used (1 means only saddlepoint);}
#'         \item{ \code{grad} }{The gradient of the log-density at y (optional);}
#'         }
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com> and Simon N. Wood.
#' @export
#'

dsaddle <- function(y, X,  decay = 0.5, deriv = FALSE, log = FALSE, 
                     normalize = FALSE, fastInit = FALSE, control = list(), 
                     multicore = !is.null(cluster), ncores = detectCores() - 1, cluster = NULL) {
  ## X[i,j] is ith rep of jth variable; y is vector of variables.
  ## evaluate saddle point approximation based on empirical CGF
  if( !is.matrix(X) ) X <- matrix(X, length(X), 1)
  
  d <- ncol(X)
  
  if( !is.matrix(y) ){ 
    if(d == 1){
      y <- as.matrix(y)
    } else {
      if(length(y) == d) y <- matrix(y, 1, d) else stop("y should be a matrix n by d.")
    }
  }
  
  if( multicore && fastInit ) stop("You can't use multicore and fastInit together")
  
  ny <- nrow( y )
  
  # Offsetting dimensionality, so decay stays pretty much at the same level for any d.
  decayI <- decay / ( d ^ 2 )
  
  # Setting up control parameter
  ctrl <- list( "method" = "LAP", 
                "nNorm" = 100 * ncol(X), 
                "tol" = 1e-6, 
                "mst" = NULL)
  
  # Checking if the control list contains unknown names
  # Entries in "control" substitute those in "ctrl"
  ctrl <- .ctrlSetup(innerCtrl = ctrl, outerCtrl = control, verbose = FALSE)
  
  if( multicore ){ 
    # Force evaluation of everything in the environment, so it will available on cluster
    .forceEval(ALL = TRUE)
    
    tmp <- .clusterSetUp(cluster = cluster, ncores = ncores, libraries = "esaddle", exportALL = TRUE)
    cluster <- tmp$cluster
    ncores <- tmp$ncores
    clusterCreated <- tmp$clusterCreated
    registerDoSNOW(cluster)
  }
  
  # Pre-calculating covariance and normalizing
  #iCov <- .robCov(t(X), alpha2 = 4, beta2 = 1.25)
  iCov <- .robCov(t(X), alpha = 10, alpha2 = 10, beta2 = 1.25)
  
  # Saving originals
  iX <- X
  iy <- y
  
  # Weighting the statistics in order to downweight outliers
  X <- iCov$weights * X
  
  # Creating normalized version
  y <- t( iCov$E %*% (t(y) - iCov$mY) )
  X <- t( iCov$E %*% (t(X) - iCov$mY) )
  
  # If there are some statistics with zero variance we remove them
  if( length(iCov$lowVar) )
  {
    stop("The columns of X indexed", iCov$lowVar, "have zero variance.")
    #y <- y[-iCov$lowVar]
    #X <- X[ , -iCov$lowVar, drop = FALSE]
  }
  
  if( fastInit && ny > 1 )
  {
    
    # Objective function
    objFun <- function(.y, .lambda)
    {
      return( .dsaddle(y = .y, X = X, tol = ctrl$tol, decay = decayI, 
                        deriv = deriv, lambda = .lambda) )
    }
    
    # Gradient function
    gradFun <- function(.y, .lambda, .extra)
    {
      return( .gradSaddle(y = .y, lambda = .lambda,  X = X, decay = decayI, onlyDlamDy = TRUE, extra = .extra) )
    }
    
    # Approximate solution of saddlepoint equation, based on Gaussian CGF: lambda = Sigma^-1 (x - mu)
    if( is.null(ctrl$mst) ) lamHat <- y else lamHat <- NULL
    
    # Evaluate density at multiple points using minimum spanning tree
    out <- mstOptim(y, lamHat, objFun, gradFun, mst = ctrl$mst)
    
  } else{
    out <- list()
    # Divide saddlepoint evaluations between cores
    withCallingHandlers({
      out <- alply(y, 1, .dsaddle, .parallel = multicore,
                   # Args for .dsaddle()
                   X = X,
                   tol = ctrl$tol,
                   decay = decayI,
                   deriv = deriv)}, warning = function(w) {
                     # There is a bug in plyr concerning a useless warning about "..."
                     if (length(grep("... may be used in an incorrect context", conditionMessage(w))))
                       invokeRestart("muffleWarning")
                   })
    
    # Close the cluster if it was opened inside this function
    if(multicore && clusterCreated) stopCluster(cluster)
    
  }
  
  # Create output list
  out <- list( "llk" = sapply(out, "[[", "llk"),
               "mix" = sapply(out, "[[", "mix"),
               "niter" = sapply(out, "[[", "niter"),
               "lambda" = t(sapply(out, "[[", "lambda")),
               "grad" = if(deriv) t( sapply(out, "[[", "DsadDy") ) else NULL)
  
  # Correct the log-likelihood for the missing normalizing constant
  if( normalize ){
    
    stopifnot( (ctrl$method %in% c("IS", "LAP")) )
    
    # Log-normalizing constant by importance sampling
    if(ctrl$method == "IS") 
    {
      aux <- rmvn(ctrl$nNorm, iCov$mY, 2*iCov$COV)
      
      logNorm <- log( .meanExpTrick( 
        dsaddle(y = aux, X = iX, decay = decay, deriv = FALSE, 
                 log = TRUE, normalize = FALSE, fastInit = fastInit, control = ctrl, 
                 multicore = multicore, ncores = ncores, cluster = cluster)$llk - dmvn(aux, iCov$mY, 2*iCov$COV, log = TRUE) )
      )
      
    }
    
    # Log-normalizing constant by Laplace approximation
    if(ctrl$method == "LAP")
    {
      tmp <- findMode(X = iX, init = iCov$mY, decay = decay, sadControl = list("tol" = ctrl$tol), hess = T)
      logNorm <- .laplApprox(tmp$logDens, tmp$hess, log = TRUE)
    }
    
    out$logNorm <- logNorm
    out$llk <- out$llk - logNorm
            
  }
  
  # Adjusting log-lik and its gradient to correct for normalization
  out$llk <- out$llk + sum(log(abs(diag(iCov$E))))
  if( deriv ) out$grad <- drop( drop(out$grad) %*% iCov$E )
  
  if( !log ) out$llk <- exp( out$llk )
  
  return( out )
  
}




##########
# INTERNAL
##########

.dsaddle <- cmpfun(function(y, X, tol, decay, 
                             deriv = FALSE, mixMethod = "mse", 
                             maxit = 100, lambda = NULL) {
  ## X[i,j] is ith rep of jth variable; y is vector of variables.
  ## evaluate saddle point approximation based on empirical CGF
  
  if( !is.vector(y) ) y <- as.vector(y)
  
  d <- length(y)
  
  if( !is.matrix(X) ){
    if(d > 1){ 
      stop("Error: simulated data must be entered in matrix form")
    }else{
      X <- matrix(X, length(X), 1)
    }
  }
  
  n <- nrow(X)
  
  # Initial guess of the root is the solution to the Gaussian case
  # the gain is one step less of Newton on average.
  if( is.null(lambda) ) lambda <- y
  
  mu <- rep(0, d)
  sig <- diag(1, d)
  
  # Choose the mixture of saddlepoint-normal, mix \in [0, 1]
  mix <- .ecgfMix(y, decay = decay, method = mixMethod, deriv = FALSE, m = d)$mix
  
  b <- .ecgf(lambda, X, kum1 = mu, kum2 = sig, mix = mix, grad = 2)
  
  ## Newton loop to minimize K(lambda) - t(lambda)%*%y or solve dK(lambda) = y wrt lambda
  kk <- jj <- 0
  # Convergence test: see [con_test] below.
  while( any( abs(b$dK-y) > tol ) && kk < maxit ) 
  { 
    kk <- kk + 1
    
    # Build scaling vector and matrix
    dd <- diag(b$d2K)^-0.5
    DD <- tcrossprod(dd, dd)
    
    d2KQR <- qr(DD * b$d2K, tol = 0)
    
    # Try solve scaled linear system fast, if that doesn't work use QR decomposition.
    d.lambda <- -  dd * drop(qr.solve(d2KQR, dd*(b$dK-y), tol = 0))
    
    lambda1 <- lambda + d.lambda ## trial lambda
    
    b1 <- .ecgf(lambda1, X, kum1 = mu, kum2 = sig, mix = mix, grad = 2)
    if ( sum( abs(b1$d2K) ) == 0 ) return(NA) ## spa breakdown (possibly too tight)
    
    jj <- 1
    c1 <- 10^-4
    alpha <- 1
    rho <- 0.5
    ## Line search checking Arminjo condition and step halving
    while( ( b1$K - crossprod(lambda1, y) ) > 
             ( b$K - crossprod(lambda, y) ) + c1 * alpha * crossprod(d.lambda, drop(b$dK-y)) && jj < 50)  
    {
      jj <- jj + 1
      alpha <- alpha * rho
      d.lambda <- d.lambda * alpha
      lambda1 <- lambda + d.lambda
      b1 <- .ecgf(lambda1, X, kum1 = mu, kum2 = sig, mix = mix, grad = 2)
    }
    
    ## now update lambda, K, dK, d2K
    lambda <- lambda1 
    b <- b1
    
  } ## end of Newton loop
  ## lambda is the SPA lambda...
  
  if(kk > 50 || jj > 20) warning(paste("The convergence of the saddlepoint root-finding is quite slow! \n",
                                       "Outer root-finding Newton-Raphson:", kk, "iter \n",
                                       "Inner line search for Arminjo condition:", jj, "iter"))
  
  # We need to recompute the QR decomposition
  # Build scaling vector and matrix
  dd <- diag(b$d2K)^-0.5
  DD <- tcrossprod(dd, dd)
  d2KQR <- qr(DD * b$d2K, tol = 0)
  
  # Determinant of b$d2K from its qr decomposition
  logDet <- sum(log(abs(diag(qr.R(d2KQR))))) - 2 * sum(log(dd))
  
  spa <- b$K - crossprod(lambda, y) - 0.5 * log( 2*pi ) * d - 0.5 * logDet
  
  # Additional stuff needed by .gradSaddle()
  b[ c("dd", "DD", "d2KQR") ] <- list(dd, DD, d2KQR)
  
  out <- list("llk" = drop(spa), "mix" = mix, "niter" = kk, "lambda" = lambda, "extra" = b)
  
  if(deriv)
  {
    
    tmp <- .gradSaddle(y = y, lambda = lambda, X = X, decay = decay, mixMethod = mixMethod, extra = b)
    
    out <- c(out, tmp)
    
  }
  
  return( out )
  
})