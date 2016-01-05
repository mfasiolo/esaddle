########
# Simulating r.v. from a saddlepoint density
########

rsaddle <- function(n, X, decay,
                    multicore = !is.null(cluster), cluster = NULL, ncores = detectCores() - 1,  ...)
{
  prop <- rmvn(n, colMeans(X), cov(X))
  
  w <- dsaddle(prop, X = X, decay = decay, multicore = multicore, ncores = ncores, cluster = cluster, ...)$llk / 
       dmvn(prop, colMeans(X), cov(X))
    
  out <- prop[ sample(1:n, n, replace = TRUE, prob = w), ]
  
  return( out )
  
}