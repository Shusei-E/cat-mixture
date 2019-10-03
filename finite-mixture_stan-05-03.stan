data {
  int<lower=1> K;            # number of clusters
  int<lower=1> D;            # number of offices
  int<lower=1> N;            # number of voters
  int<lower=0, upper=1> y[N, D];            # observations
  vector<lower=0>[K] alpha;  # hyperparameter
}

transformed data {
  # following stan 9.2 but not sure if this is relevant in bernoulli case
  real<upper=0> neg_log_K;
  neg_log_K = -log(K);
}


parameters {
  simplex[K] pi;         # mixing proportions
  real<lower=0,upper=1> mu[K, D];       # Bernoulli probability
}

transformed parameters {

  real<upper=0> soft_z[N, K]; # log unnormalized clusters
  for (n in 1:N) {
    for (k in 1:K) {
      soft_z[n, k] = neg_log_K - 0.5*dot_self(mu[k] - y[n]);
    }
  }
}



model {

  // likelihood
    for (n in 1:N) {
      vector[K] lps = rep_vector(neg_log_K, K);
      for (k in 1:K) {
        for(j in 1:D) {
          lps[k] += bernoulli_lpmf(y[n, j] | mu[k, j]);
        }
      }

      #  z[n] ~ categorical(rep(1/K, K))
      #  y[n] ~ bernoulli(mu[z[n]])
  }


  // prior
  for (k in 1:K)
      mu[k] ~ beta(2, 5); # something with mode at a low pr

}