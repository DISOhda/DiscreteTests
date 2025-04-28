#include <Rcpp.h>
using namespace Rcpp;

// Wilcoxon's sign rank probabilities
// [[Rcpp::export]]
List sign_rank_probs_int(
    const NumericVector n  // must be sorted!!
) {
  // number of distributions to be computed
  int len = n.length();
  // largest sample size
  int max_n = n[len - 1];
  // largest possible observation for largest sample size
  int max_obs = max_n * (max_n + 1) / 2;

  // list of results
  List out(len);

  // find zeros, ones and twos in sizes and store distributions
  int pos_out = 0;
  while(pos_out < len && n[pos_out] == 0) {
    out[pos_out] = NumericVector(1, 1.0);
    pos_out++;
  }
  while(pos_out < len && n[pos_out] == 1) {
    out[pos_out] = NumericVector(2, 0.5);
    pos_out++;
  }

  if(pos_out < len) {
    NumericVector dist(max_obs + 1);
    dist[0] = dist[1] = 0.5;

    // compute non-trivial distributions
    for(int k = 2; k <= max_n; k++) {
      int end_new = k * (k + 1) / 2;

      for(int j = end_new; j >= k; j--) {
        dist[j - k] /= 2;
        dist[j] += dist[j - k];
      }

      while(pos_out < len && n[pos_out] == k){
        out[pos_out] = dist[Range(0, end_new)];
        pos_out++;
      }
    }
  }

  // return results
  return out;
}
