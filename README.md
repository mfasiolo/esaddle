This package implements the Extended Empirical Saddlepoint density approximation (ESA) described in XXX.
The documentation is fairly sparse, but most of the information is contained in the paper.

The main functions are:
- `dsaddle` which can be used to evaluate the ESA, given some data.
- `selectDecay` which can be used to select the tuning parameter of ESA by cross-validation.
- `findMode` which maximizes ESA in order to find its mode.

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
rug(x)
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
plot(X, pch=".", col=1, ylim = c(min(x2), max(x2)), xlim = c(min(x1), max(x1)),
     main = "ESA density", xlab = expression(X[1]), ylab = expression(X[2]))
contour(x1, x2, matrix(dw, m, m), levels = quantile(as.vector(dw), seq(0.8, 0.995, length.out = 10)), col=2, add=T)

# Plotting ESA density
par(mfrow = c(1, 2))
plot(X, pch=".",col=2, ylim = c(min(x2), max(x2)), xlim = c(min(x1), max(x1)),
     main = "True Density", xlab = expression(X[1]), ylab = expression(X[2]))
contour(x1, x2, matrix(dwa, m, m), levels = quantile(as.vector(dwa), seq(0.8, 0.995, length.out = 10)), col=2, add=T)

# Finding mode using ESA 
init <- rnorm(2, 0, sd = c(1, rep(2, 2-1))) # random initialization
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

The model contains the real dataset described in ([Fasiolo and Wood, 2015](http://arxiv.org/abs/1511.02644)),
which can be loaded using `data(voles_data)`. 
To fit the model with Synthetic Likelihood (Wood, 2010), we need to load the [synlik package](https://cran.r-project.org/web/packages/synlik/index.html) and
then create an object of class `synlik`. Notice that here we use the wrapper `volesWrap`, rather than `volesSimulator`, 
for compatibility with `synlik`.

```R
# This is needed because the process is observed only in 
# the months of June and September, for 45 years
obsTiming <- sort(c(6 + seq(0, 12 * 44, by = 12), 9 + seq(0, 12 * 44, by = 12)))

# Create "synlik" object
voles_sl <- new("synlik",
                simulator = volesWrap,
                summaries = volesStats,
                param = log(c(r = 4.5, e = 1, g = 0.2, h = 0.15, a = 8,
                              d = 0.06, s = 1, sigmaProc = 1.5, phi = 100)),
                extraArgs = list("nObs" = 12 * 45,  
                                 "nBurn" = 12 * 10, 
                                 "monthsToSave" = obsTiming)
)

# Put real data into object
data(voles_data)
voles_sl@data <- round(voles_data$popDensity * 10)
voles_sl@extraArgs$obsData <- round(voles_data$popDensity * 10)
```

We can than fit the model using a Metropolis-Hastings sampler.

```R
# Load MCMC proposal
data(voleFullProp_sl)

# Run MCMC on synthetic likelhood (only 100 iterations and it takes a while)
set.seed(51554)
tim <- proc.time()
voles_sl <- smcmc(voles_sl, 
                  initPar = voleFullTrueInit_sl,
                  nsim = 1000,
                  niter = 100, 
                  burn = 0,
                  priorFun = voleFullPrior_sl, 
                  propCov = 0.1 * voleFullProp_sl)
)
tim <- proc.time() - tim

plot(voles_sl)
```


References
----------------------------
  
  * Fasiolo, M and Wood, S. N. (2015). Approximate methods for dynamic ecological models, 
  ArXiv http://arxiv.org/abs/1511.02644. To appear in the Handbook of Approximate Bayesian Computation by S. Sisson, L. Fan, and M. Beaumont.

  * Turchin, P. and S. P. Ellner (2000). Living on the edge of chaos: population
dynamics of fennoscandian voles. Ecology 81 (11), 3099–3116.

  * Wood, S. N. (2010). Statistical inference for noisy nonlinear ecological dynamic systems. Nature 466 (7310), 1102–1104.
