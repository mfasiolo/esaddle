#####
#' Empirical multivariate Gaussian density 
#' @description Gives a pointwise evaluation of the multivariate normal (MVN) fit to X at position y.
#'
#' @param y points at which the MVN is evaluated (d dimensional vector) or an n by d matrix, each row indicating
#'          a different position.
#' @param X n by d matrix containing the data.
#' @param log If TRUE the log of the log-density is returned.
#' @param verbose currently not used.
#' @return The density of the empirical MVN, evaluated at y.
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com> and Simon Wood.
#' @examples 
#' X <- matrix(rnorm(2 * 1e3), 1e3, 2)
#' demvn(rnorm(2), X, log = TRUE)
#' @export demvn
#'
#'
demvn <- function(y, X, log = FALSE, verbose = TRUE)
{
  
  if( !is.matrix(y) ) y <- matrix(y, 1, length(y))
  
  tmp <- robCov( t(X) )
  
  # If there are some statistics with zero variace we remove them
  if( length(tmp$lowVar) ) stop("The columns of X indexed ", preCov$lowVar, " have zero variance.") #y <- y[-tmp$lowVar]
  
  llk <- apply(y, 
               1, 
               function(input) .demvn(y = input, L = tmp, log = log) )
  
  return(llk)
}


.demvn <- function(y, L, log)
{
  rss <- sum( (L$E%*%as.vector(y-L$mY))^2 )
  
  llk <- -rss/2 - L$half.ldet.V - log(2 * pi) * length(y) / 2
  
  if( !log ) llk <- exp( llk )
  
  return(llk)
}