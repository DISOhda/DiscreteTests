#include <Rcpp.h>
using namespace Rcpp;

// Wilcoxon's sign rank probabilities
// [[Rcpp::export]]
List sign_rank_probs_int(
    const NumericVector n  // must be sorted!!
) {
  // number of distributions to be computed
  int32_t len = n.length();
  // largest sample size
  int32_t max_n = n[len - 1];
  // largest possible observation for largest sample size
  int32_t max_obs = max_n * (max_n + 1) / 2;

  // list of results
  List out(len);

  // find zeros, ones and twos in sizes and store distributions
  int32_t pos_out = 0;
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
    for(int32_t k = 2; k <= max_n; k++) {
      int32_t end_new = k * (k + 1) / 2;

      for(int32_t j = end_new; j >= k; j--) {
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

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector sign_rank_probs_ties_int(const IntegerVector scores) {
  IntegerVector s_sorted = clone(scores);
  std::sort(s_sorted.begin(), s_sorted.end());

  uint64_t total = 0;
  for (int32_t v : s_sorted) total += v;

  std::vector<double> dp(total + 1, 0.0);
  dp[0] = 1.0;
  uint64_t last = 0;

  for(int32_t idx = 0; idx < s_sorted.size(); idx++) {
    uint64_t s = s_sorted[idx];

    for(uint64_t j = last + s; j >= s; j--) {
      dp[j - s] *= 0.5;
      dp[j]     += dp[j - s];
    }

    last += s;
  }

  NumericVector res(total + 1);
  double sum_total = 0.0;
  for(uint64_t i = 0; i <= total; i++) sum_total += dp[i];
  for(uint64_t i = 0; i <= total; i++) res[i] = dp[i] / sum_total;

  return res;
}
