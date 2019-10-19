#include <Rcpp.h>
using namespace Rcpp;

// (mu^{y = l})
// [[Rcpp::export]]
double mu_ind(const double &mu, const int &y, const int &l) {
  bool ind = (y == l);
  return pow(mu, ind);
}

// for a given individual i, compute the numerators of the
// responsibility parameter zeta_ik for a given cluster k,
// using the length-D vector of votes, each of which can take
// values 0:L
// [[Rcpp::export]]
double mu_yvec(NumericMatrix mu_k, IntegerVector y, int L) {
  int D = y.size();
  double total = 1;

  for (int j = 0; j < D; ++j) {
    for (int l = 0; l <= L; ++l) {
      total *= mu_ind(mu_k(j, l), y[j], l);
    }
  }
  return total;
}
