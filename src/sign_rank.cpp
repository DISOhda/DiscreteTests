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
    // array of vectors of NumericVectors
    NumericVector dists[2] = {
      NumericVector(max_obs + 1),
      NumericVector(max_obs + 1)
    };
    dists[0][0] = dists[0][1] = 0.5;

    // which array field contains previous and new distribution?
    int olds = 1, news = 0;

    // compute non-trivial distributions
    for(int k = 2; k <= max_n; k++) {
      olds = 1 - olds;
      news = 1 - news;

      int end_olds = k * (k - 1) / 2;
      int end_news = k * (k + 1) / 2;

      for(int j = 0; j <= end_olds; j++) dists[news][j] = dists[olds][j] / 2;
      for(int j = end_news; j >= k; j--) dists[news][j] += dists[news][j - k];

      while(pos_out < len && n[pos_out] == k){
        out[pos_out] = dists[news][Range(0, end_news)];
        pos_out++;
      }
    }
  }

  // return results
  return out;
}

