data {
  int<lower=1> K;            # number of clusters
  int<lower=1> D;            # number of offices
  int<lower=1> N;            # number of voters
  vector[D] y[N];            # observations
  vector<lower=0>[K] alpha;  # hyperparameter
}

transformed data {
  # following stan 9.2 but not sure if this is relevant in bernoulli case
  real<upper=0> neg_log_K;
  neg_log_K = -log(K);
}


parameters {
  simplex[K] pi;         # mixing proportions
  vector[D] mu[K];       # Bernoulli probability
}

transformed parameters {

  real<upper=0> soft_z[N, K]; # log unnormalized clusters
  for (n in 1:N)
    for (k in 1:K) {
      # I do not understand this part from user manual 9.2
      soft_z[n, k] = neg_log_K - 0.5*dot_self(mu[k] - y[n]);
    }

}



model {
  // likelihood
    for (n in 1:N) {
        z[n] ~ categorical(rep(1/K, K))
        y[n] ~ bernoulli(mu[z[n]])
  }


  // prior
  for (k in 1:K)
      mu[k] ~ beta(2, 5); # something with mode at a low pr

}