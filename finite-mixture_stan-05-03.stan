data {
  int<lower=2> M;          // number of possible values of Y
  int<lower=1> K;          // number of clusters
  int<lower=1> D;          // number of offices
  int<lower=1> N;          // number of voters
  vector[D] y[N];           // observations
  vector<lower=0>[D] alpha;  // hyperparameter
}

transformed data {
  real<upper=0> neg_log_K;
  neg_log_K = -log(K);
}


parameters {
  simplex[K] pi;         // mixing proportions
  vector[D] mu[K];            // Bernoulli probability
}

transformed parameters {
  real<upper=0> soft_z[N, K]; // log unnormalized clusters
  for (n in 1:N)
    for (k in 1:K)
      soft_z[n, k] = neg_log_K - 0.5*dot_self(mu[k] - y[n]);
}

model {
  vector[K] log_pi = log(pi); // cache log calculutation

  // prior
  for (k in 1:K)
      mu[k] ~ beta(2, 5); // something with mode at a small value

  // likelihood
    for (n in 1:N) {
        target += log_sum_exp(soft_z[n]);
  }
}