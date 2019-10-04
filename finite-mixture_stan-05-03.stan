data {
  int<lower=1> K;                // number of clusters
  int<lower=1> D;                // number of offices
  int<lower=1> N;                // number of voters
  real y[N, D]; // data
  vector<lower=0>[K] alpha;      // hyperparameter for cluster prevalence
}


parameters {
  simplex[K] theta;          // mixing proportions
  ordered[D] mu[K]; // Normal mean
}

model {
  // prior
  theta ~ dirichlet(alpha);  // prior for theta

	for(k in 1:K){
		mu[k] ~ normal(0, 3);
	}

  // likelihood
	real lps[K];
  for (n in 1:N) {
      // initiate with z ~ Cat(theta), theta ~ Dir(alpha)
      for (k in 1:K) {
        // sum all possible values of
        lps[k] = normal_lpdf(y[n, ] | mu[j], 0.05) + 
				          log(theta[k]);   // prior for z
      }
      target += log_sum_exp(lps);
  }



}
