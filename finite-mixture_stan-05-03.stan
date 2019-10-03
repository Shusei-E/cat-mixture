data {
  int<lower=1> K;                // number of clusters
  int<lower=1> D;                // number of offices
  int<lower=1> N;                // number of voters
  int<lower=0, upper=1> y[N, D]; // data
  vector<lower=0>[K] alpha;      // hyperparameter for cluster prevalence
}

transformed data {
  real<upper=0> neg_log_K;
  neg_log_K = -log(K);
}


parameters {
  simplex[K] pi_mixture;          // mixing proportions
  real<lower=0,upper=1> mu[K, D]; // Bernoulli probability
}

model {
  // prior
  pi_mixture ~ dirichlet(alpha);
  for (k in 1:K) {
    for (j in 1:D) {
      // something with mode at a low pr
      mu[k, j] ~ beta(2, 5);
    }
  }

  // likelihood
  for (n in 1:N) {
    for(j in 1:D) {
      // initiate with z ~ Cat(1/K, ... 1/K)
      vector[K] lps = rep_vector(neg_log_K, K);
      for (k in 1:K) {
        // sum all possible values of
        lps[k] += bernoulli_lpmf(y[n, j] | mu[k, j]) + log(pi_mixture[k]);
      }
      target += log_sum_exp(lps);
    }
  }



}
