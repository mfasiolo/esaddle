#####
#' Simulate random variables from Extended Empirical Saddlepoint density (ESA)
#' @description Simulate random variables from the Extended Empirical Saddlepoint density (ESA), using importance 
#'              sampling and then resampling according to the importance weights.
#'
#' @param n number of simulated vectors.
#' @param X m by d matrix containing the data.
#' @param decay Rate at which the SPA falls back on a normal density. Should be a positive number,
#'              by default set to 0.5.
#' @param multicore  if TRUE each fold will be run on a different core.
#' @param ncores   number of cores to be used.
#' @param cluster an object of class \code{c("SOCKcluster", "cluster")}. This allowes the user to pass her own cluster,
#'                which will be used if \code{multicore == TRUE}. The user has to remember to stop the cluster. 
#' @param ... additional arguments to be  A list of control parameters with entries:
#'         \itemize{
#'         \item{ \code{method} }{The method used to calculate the normalizing constant. 
#'                                Either "LAP" (laplace) or "IS" (importance sampling).}
#'         \item{ \code{nNorm} }{If control$method == "IS", this is the number of importance samples used.}
#'         \item{ \code{tol} }{The tolerance used to assess the convergence of the solution to the saddlepoint equation.
#'                             The default is 1e-6.}
#'         }
#' @return An n by d matrix containing the simulated vectors.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>.
#' @examples
#' X <- matrix(rexp(2 * 1e3), 1e3, 2)
#' Z <- rsaddle(500, X, decay = 0.1)
#' hist( Z[ , 1] )
#' @export
#'


rsaddle <- function(n, X, decay,
                    multicore = !is.null(cluster), cluster = NULL, ncores = detectCores() - 1,  ...)
{
  prop <- rmvn(n, colMeans(X), 1.5 * cov(X))
  
  w <- dsaddle(prop, X = X, decay = decay, multicore = multicore, ncores = ncores, cluster = cluster, ...)$llk / 
       dmvn(prop, colMeans(X), cov(X))
    
  out <- prop[ sample(1:n, n, replace = TRUE, prob = w), ]
  
  return( out )
  
}