library(purrr)
library(rstan)
library(brms)
library(glue)

set.seed(02138)

# dimensions
K <- 4L
D <- 4L
N <- 600L

# hyperparameter
alpha <- rep(1, K)


# cluster assignment
theta <- rdirichlet(1, alpha)
Z_table <- rmultinom(N, 1, theta)
Z <-  map_dbl(1:N,  ~which(Z_table[, .x] == 1))

# set parameters
mu <- list(
  `1` = rnorm(D, 0, 3),
  `2` = rnorm(D, 0, 3),
  `3` = rnorm(D, 0, 3),
  `4` = rnorm(D, 0, 3)
)

# each mu vector is K by D
stopifnot(all(map_lgl(mu, ~length(.x) == D)))


# Generate data
y <- array(NA, dim = c(N, D))
for (i in 1:N) {
  # y[i, ] <- rbinom(D, size = 1, prob = mu[[Z[i]]])
  y[i, ] <- rnorm(D, mean = mu[[Z[i]]], sd=0.05)
}

# put together data
data <- list(D = D,
             K = K,
             N = N,
             y = y,
             alpha = alpha)

# target params
params <- list(theta = theta,
               mu = mu,
               Z = Z)


# check vanilla k means
k_vanilla <- kmeans(data$y, centers = K)
kdf_vanilla <- as_tibble(k_vanilla$centers) %>%
  mutate(n = k_vanilla$size)

# Gaussian Mixture
library(mclust)
gmm <- Mclust(y, G=K)
summary(gmm)
table(Z)

# EM
X <- y  # data
init_mu <- t(matrix(apply(X, 2, mean)))
mean_ <- matrix(rep(init_mu, K), ncol=D)
cov_ <- array(rep(0,D*D*K), dim=c(D, D, K)) # Multi dimensional array
for(k in 1:K) cov_[,,k] = diag(D)

pi_ <- matrix(rep(1, K), nrow=K, ncol=1) / K

gamma_ <- matrix(rep(0, N*K), nrow=N, ncol=K)

## Iteration
iter <- 0

likelihood <- function(X, mean_, cov_, pi_){
  # log likelihood
  summation <- 0.0
  for(n in 1:dim(X)[1]){
    temp = 0.0
    for(k in 1:K){
      temp <- temp + pi_[k] * dmvnorm(t(X[n,]), t(mean_[k,]), cov_[,,k])
    }
    summation <- summation + log(temp)
  }
  return(summation)
}

like <- likelihood(X, mean_, cov_, pi_)

while(1){
  ## E-Step
  for(n in 1:N){
    denominator <- 0.0
    for(k in 1:K)
      denominator <- denominator + pi_[k] * dmvnorm(t(X[n,]), t(mean_[k,]), cov_[,,k])

    for(k in 1:K)
      gamma_[n,k] <- pi_[k] * dmvnorm(t(X[n,]), t(mean_[k,]), cov_[,,k]) / denominator
  }

  ## M-Step (intuitively, MLE with weights)
  for(k in 1:K){
    Nk <- sum(gamma_[,k])

    # Calculate mean
    mean_[k,] <- rep(0,D)
    for(n in 1:N)
      mean_[k, ] <- mean_[k, ] + gamma_[n,k] * X[n,]
    mean_[k, ] <- mean_[k, ] / Nk

    # Calculate covariance
    cov_[,,k] <- array(rep(0,D*D), dim=c(D, D))
    for(n in 1:N){
      temp <- matrix(X[n,] - mean_[k,])
      cov_[,,k] <- cov_[,,k] + gamma_[n,k] * temp %*% t(temp) 
    }
    cov_[,,k] <- cov_[,,k] / Nk

    # Calculate mixture proportion
    pi_[k] <- Nk / N
  }

  ## Judge convergence
  new_like <- likelihood(X, mean_, cov_, pi_)
  difference <- abs(new_like - like)
  if(difference < 0.001) break
  like <- new_like
  iter <- iter+1
  print(c(as.integer(iter), like))
}

pi_
table(Z)/sum(table(Z))

mean_
mu
