data {
  int<lower=1> K;                // number of clusters
  int<lower=1> D;                // number of offices
  int<lower=1> N;                // number of voters
  int<lower=0, upper=1> y[N, D]; // data
  real<lower = 0> alpha;
}

transformed data {
  real<upper=0> neg_log_K;
  neg_log_K = -log(K);
}


parameters {
  simplex[K] pi;                  // mixing proportions
  real<lower=0,upper=1> mu[K, D]; // Bernoulli probability
}


model {
  // likelihood
  for (n in 1:N) {
    for(j in 1:D) {
      vector[K] lps = rep_vector(neg_log_K, K); // assume z ~ Cat(1/K, ... 1/K)
      for (k in 1:K) {
        lps[k] += bernoulli_lpmf(y[n, j] | mu[k, j]); // sum all possible values of k
        // lps[k] += bernoulli_lpmf(y[n, j] | mu[k, j]) + log(pi[k]); // sum all possible values of k
      }
      target += log_sum_exp(lps);
    }
  }

  // prior
  for (k in 1:K) {
    for (j in 1:D) {
      // something with mode at a low pr
      mu[k, j] ~ beta(1, 1);
    }
  }
  
  pi ~ dirichlet(alpha); 
}
