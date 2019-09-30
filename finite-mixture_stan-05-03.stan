data {
  int<lower=2> M;          // number of possible values of Y
  int<lower=1> K;          // number of clusters
  int<lower=1> J;          // number of offices
  int<lower=1> N;          // number of voters
  int  Y[N, J];          // observations
  vector<lower=0>[M] alpha;  // hyperparameter
}


parameters {
  simplex[M] theta[K*J];     // abstain, straight, split
  simplex[K] psi[J];         // mixing proportions
}

model {
  for (j in 1:J) {
      psi[j] ~ dirichlet(alpha);  // prior
  }

  for (j in 1:J) {
    for (n in 1:N) {
      vector[K] lps = log(psi[j]);
      for (k in 1:K) {
        lps[k] += categorical_lpmf(Y[n, j] | theta[k*(j - 1) + k]);
      }
      target += log_sum_exp(lps);
    }
  }
}