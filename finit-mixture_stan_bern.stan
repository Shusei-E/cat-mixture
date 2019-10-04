data {
  int<lower=1> K;                // number of clusters
  int<lower=1> D;                // number of offices
  int<lower=1> N;                // number of voters
  int<lower=0, upper=1> y[N, D]; // data
  vector<lower=0>[K] alpha;      // hyperparameter for cluster prevalence
}


parameters {
  simplex[K] theta;          // mixing proportions
  real<lower=0,upper=1> mu[K, D]; // Bern probability
}

model {
  // prior
  theta ~ dirichlet(alpha);  // prior for theta

  // likelihood
  for (n in 1:N) {
    for(j in 1:D) {
      // initiate with z ~ Cat(theta), theta ~ Dir(alpha)
      real lps[K];
      for (k in 1:K) {
        // sum all possible values of
        lps[k] = bernoulli_lpmf(y[n, j] | mu[k, j]) + 
				          log(theta[k]) +  // prior for z
									beta_lpdf(mu[k, j] | 2.0, 5.0);  // prior for mu
      }
      target += log_sum_exp(lps);
    }
  }



}

