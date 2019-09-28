data {
  int<lower=2> M;          // number of possible values of Y
  int<lower=1> K;          // number of clusters
  int<lower=1> J;          // number of offices
  int<lower=1> N;          // number of data points
  real<lower=0> alpha[M];  // hyperparameter
  matrix[N, J] Y;          // observations
}


parameters {
  simplex[M] theta[K*J];     // abstain, straight, split
  simplex[K] psi[J];         // mixing proportions
}

model {
  vector[K] log_psi = log(psi);  // log of mixing proportions
  for (j in 1:J) {
    for (n in 1:N) {
      vector[K] lps = log_psi;
      for (k in 1:K) {
        lps[k] += multinomial_lpmf(Y[n, j] | theta[k*(j - 1) + k]);
      }
      target += log_sum_exp(lps);
    }
  }
}