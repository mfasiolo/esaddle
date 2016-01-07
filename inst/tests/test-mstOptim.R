context("mstOptim")

###############################################################################
test_that("mstOptim works well with saddlepoint: Univariate case", {
  
  x <- rgamma(1000, 1, 1)
  decay <-  0.2

  n <- 1000
  xSeq <- seq(-2, 10, length.out =n)
  
  slow <- dsaddle(y = xSeq, X = x,  decay = decay, deriv = FALSE, log = TRUE)
  fast <- dsaddle(y = xSeq, X = x,  decay = decay, deriv = FALSE, log = TRUE, fastInit = TRUE)
  
  expect_less_than( max(abs(slow$llk - fast$llk) / abs(slow$llk)), 1e-4 )
  
  plot(slow$mix, slow$niter - fast$niter, xlab = "Sad-Gaus Mix", ylab = "Iterations saved", 
       main = "Univariate Case")  
  abline(h = 0, col = 2)
})

###############################################################################
test_that("mstOptim works well with saddlepoint: Multivariate Gaussian case", {
  
  library(mvnfast)
  
  ### Simulate Data
  dims <- 2
  decay <- 1
  n <- 1000
  A <- matrix(rnorm(dims^2), dims, dims)
  A <- t(A)%*%A + diag(0.5, dims)
  myMu <- seq.int(1, dims)
  SIGMA <- A 
  X <- rmvn(n = n, mu = myMu, sigma = SIGMA)
  
  #### Points at which saddlepoint is evaluated
  Y <- X[1:n, , drop = F]
  ny <- nrow(Y)
  
  # Pre-compute minimum spanning tree
  preCov <- esaddle:::.robCov(t(Y), alpha = 10, alpha2 = 10, beta2 = 1.25)
  lamHat <- t( crossprod(preCov$E, preCov$E %*% (t(Y) - preCov$mY)) )
  
  mat <- .Call("mst", X_ = t(lamHat), PACKAGE = "esaddle");
  mat <- mat[1:2, ] + 1
  
  # Convert mst to a list easier to traverse
  mst <- rep(list(NULL), ny)
  names(mst) <- 1:ny
  
  for(ii in 1:ncol(mat))
  {
    mst[[ mat[1, ii] ]] <- c(mst[[ mat[1, ii] ]], mat[2, ii])
    mst[[ mat[2, ii] ]] <- c(mst[[ mat[2, ii] ]], mat[1, ii])
  }

  ###
  # Estimate density
  slow <- dsaddle(y = Y, X = X, decay = decay, deriv = FALSE, log = TRUE)
  fast <- dsaddle(y = Y, X = X,  decay = decay, deriv = FALSE, log = TRUE, fastInit = TRUE)
  superFast <- dsaddle(y = Y, X = X,  decay = decay, deriv = FALSE, log = TRUE, 
                       fastInit = TRUE, control = list("mst" = mst))
  
  expect_less_than( max(abs(slow$llk - fast$llk) / abs(slow$llk)), 1e-4 )
  expect_less_than( max(abs(superFast$llk - fast$llk) / abs(superFast$llk)), 1e-4 )
  
  plot(slow$mix, slow$niter - fast$niter, xlab = "Sad-Gaus Mix", ylab = "Iterations saved", 
       main = "Multivariate Gaussian Case")
  abline(h = 0, col = 2)
  
})


###############################################################################
test_that("mstOptim works well with saddlepoint: Multivariate Gamma case", {
    
  ### Simulate Data
  dims <- 2
  decay <- 1
  n <- 1000
  X <- matrix(rgamma(n*dims, 4, 1),  n, dims)
  
  #### Points at which saddlepoint is evaluated
  Y <- X[1:n, , drop = F]
  ny <- nrow(Y)
  
  # Pre-compute minimum spanning tree
  preCov <- esaddle:::.robCov(t(Y), alpha = 10, alpha2 = 10, beta2 = 1.25)
  lamHat <- t( crossprod(preCov$E, preCov$E %*% (t(Y) - preCov$mY)) )
  
  mat <- .Call("mst", X_ = t(lamHat), PACKAGE = "esaddle");
  mat <- mat[1:2, ] + 1
  
  # Convert mst to a list easier to traverse
  mst <- rep(list(NULL), ny)
  names(mst) <- 1:ny
  
  for(ii in 1:ncol(mat))
  {
    mst[[ mat[1, ii] ]] <- c(mst[[ mat[1, ii] ]], mat[2, ii])
    mst[[ mat[2, ii] ]] <- c(mst[[ mat[2, ii] ]], mat[1, ii])
  }
  
  ###
  # Estimate density
  slow <- dsaddle(y = Y, X = X, decay = decay, deriv = FALSE, log = TRUE)
  fast <- dsaddle(y = Y, X = X,  decay = decay, deriv = FALSE, log = TRUE, fastInit = TRUE)
  superFast <- dsaddle(y = Y, X = X,  decay = decay, deriv = FALSE, log = TRUE, 
                       fastInit = TRUE, control = list("mst" = mst))
      
  expect_less_than( max(abs(slow$llk - fast$llk) / abs(slow$llk)), 1e-4 )
  expect_less_than( max(abs(superFast$llk - fast$llk) / abs(superFast$llk)), 1e-4 )
  
  plot(slow$mix, slow$niter - fast$niter, xlab = "Sad-Gaus Mix", ylab = "Iterations saved", 
       main = "Multivariate Gamma Case")
  abline(h = 0, col = 2)
  
})

