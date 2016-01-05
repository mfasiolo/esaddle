
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