
This package implements the Extended Empirical Saddlepoint density approximation (ESA) described in  
[Fasiolo et al., 2016](http://arxiv.org/abs/1601.01849). The documentation is fairly sparse, but most of the information is contained in the paper.

The main functions are:
- `dsaddle` which can be used to evaluate the ESA, given some data.
- `selectDecay` which can be used to select the tuning parameter of ESA by cross-validation.
- `findMode` which maximizes ESA in order to find its mode.

Here we describe how to use them with two simple examples. In the first example we the ESA
density to approximate an univariate Gamma(2, 1) density.

```R
library(esaddle)

###################################
######### Univariate Gamma example
###################################

########## Simulating data
x <- rgamma(1000, 2, 1)

# Fixing tuning parameter of ESA
decay <-  0.05

# Evaluating ESA at several point
xSeq <- seq(-2, 8, length.out = 200)
tmp <- dsaddle(y = xSeq, X = x, decay = decay, log = TRUE)

# Plotting true density, ESA and normal approximation
plot(xSeq, exp(tmp$llk), type = 'l', ylab = "Density", xlab = "x")
lines(xSeq, dgamma(xSeq, 2, 1), col = 3)
lines(xSeq, dnorm(xSeq, mean(x), sd(x)), col = 2)
suppressWarnings( rug(x) )
legend("topright", c("ESA", "Truth", "Gaussian"), col = c(1, 3, 2), lty = 1)

# Finding mode of ESA
res <- findMode(x, init = mean(x), decay = decay)$mode
abline(v = res, lty = 2, lwd = 1.5)

# Selection tuning parameter of ESA by 2-fold cross-validation on two cores
# Complexity decreases as decay increases
tmp <- selectDecay(decay = c(5e-4, 1e-3, 5e-3, 0.01, 0.1, 0.5, 5, Inf), 
                   K = 2,
                   simulator = function() x,
                   multicore = T,
                   ncores = 2)
                  
```

In the second example we consider a banana shape or warped bivariate Gaussian
density.

```R

############################################
####### Bivariate "Banana" example 
############################################

# Function that evaluates the true density
dwarp <- function(x, alpha) {
  ## warped normal density
  d <- length(alpha) + 1
  lik <- dnorm(x[ , 1], log = TRUE)
  
  tmp <- x[ , 1]^2
  for(ii in 2:d)
    lik <- lik + dnorm(x[ , ii] - alpha[ii-1]*tmp, log = TRUE)
  
  lik
}

# Function that simulates from true distribution
rwarp <- function(n = 1, alpha) {
  d <- length(alpha) + 1
  z <- matrix(rnorm(n*d), n, d)
  
  tmp <- z[ , 1]^2
  for(ii in 2:d) z[ , ii] <- z[ , ii] + alpha[ii-1]*tmp
  z
}

# Creating 2d grid
m <- 50
expansion <- 1
x1 <- seq(-2, 3, length=m)* expansion; 
x2 <- seq(-3, 3, length=m) * expansion
x <- expand.grid(x1, x2) 

# Evaluating true density on grid
alpha <- 1
dw <- dwarp(x, alpha = alpha)

# Simulate random variables
X <- rwarp(1000, alpha = alpha)

# Evaluating ESA density
dwa <- dsaddle(as.matrix(x), X, decay = 0.1, log = FALSE)$llk

# Plotting true density
par(mfrow = c(1, 2))
plot(X, pch=".", col=1, ylim = c(min(x2), max(x2)), xlim = c(min(x1), max(x1)),
     main = "True density", xlab = expression(X[1]), ylab = expression(X[2]))
contour(x1, x2, matrix(dw, m, m), levels = quantile(as.vector(dw), seq(0.8, 0.995, length.out = 10)), col=2, add=T)

# Plotting ESA density
plot(X, pch=".",col=2, ylim = c(min(x2), max(x2)), xlim = c(min(x1), max(x1)),
     main = "ESA density", xlab = expression(X[1]), ylab = expression(X[2]))
contour(x1, x2, matrix(dwa, m, m), levels = quantile(as.vector(dwa), seq(0.8, 0.995, length.out = 10)), col=2, add=T)

# Finding mode using ESA 
init <- rnorm(2, 0, sd = c(1, 2)) # random initialization
res <- findMode(X = X, init = init, decay = decay)$mode
points(res[1], res[2], pch = 3, lwd = 2)

par(mfrow = c(1, 1))
# Selection tuning parameter of ESA by 2-fold cross-validation on 2 cores
# Complexity decreases as decay increases
tmp <- selectDecay(decay = c(0.005, 0.01, 0.1, 0.25, 0.5, 1, 5, Inf), 
                   K = 2,
                   simulator = function() X,
                   multicore = T,
                   ncores = 2)
```

References
----------------------------
  
  * Fasiolo, M., Wood, S. N., Hartig, F. and Bravington, M. V. (2016). An Extended Empirical Saddlepoint Approximation for Intractable Likelihoods. ArXiv http://arxiv.org/abs/1601.01849.
  