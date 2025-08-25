#include <Rcpp.h>
using namespace Rcpp;

// Mann-Whitney probabilities
// [[Rcpp::export]]
List mann_whitney_probs_int(
  const IntegerVector m, // must be sorted!!
  const IntegerVector n  // each number corresponds to a value in 'm' and m >= n!!
) {
  // number of distributions to be computed
  int len = m.length();
  // largest size of first samples
  int max_m = m[len - 1];
  // largest size of second samples
  int max_n = max(n);

  // list of results
  List out(len);

  // compute distributions
  if(max_n == 0) {
    // if all sample pairs contain zeros, all distributions are trivial
    for(int i = 0; i < len; i++) out[i] = NumericVector(1, 1.0);
  } else {
    // determine maximums of second sample sizes for efficiency
    IntegerVector m_unique;
    IntegerVector n_unique;
    int len_unique = 1;
    if(len > 1) {
      m_unique = IntegerVector(len);
      n_unique = IntegerVector(len);
      m_unique[0] = m[0];
      n_unique[0] = n[0];
      for(int i = 1; i < len; i++) {
        if(m[i] > m[i - 1]) {
          len_unique++;
          m_unique[len_unique - 1] = m[i];
        }
        if(n_unique[len_unique - 1] <= n[i]) n_unique[len_unique - 1] = n[i];
      }
      m_unique = m_unique[Range(0, len_unique - 1)];
      n_unique = n_unique[Range(0, len_unique - 1)];
      n_unique = IntegerVector(rev(n_unique));
      n_unique = IntegerVector(cummax(n_unique));
      n_unique = IntegerVector(rev(n_unique));
      n_unique = IntegerVector(pmin(n_unique, m_unique));
    } else {
      m_unique = IntegerVector(1, m[0]);
      n_unique = IntegerVector(1, n[0]);
    }

    // find zeros in size pairs and store trivial distribution
    int pos_pair = 0;
    while(pos_pair < len && (m[pos_pair] == 0 || n[pos_pair] == 0)) {
      out[pos_pair] = NumericVector(1, 1.0);
      pos_pair++;
    }

    // find first unique non-zero second sample size
    int pos_n_unique = 0;
    while(pos_n_unique < len_unique && n_unique[pos_n_unique] == 0)
      pos_n_unique++;

    // array of NumericVectors for all distributions (1, ..., n) for m
    NumericVector* dist = new NumericVector[(unsigned int)max_n];
    for(int i = 0; i < max_n; i++) dist[i] = NumericVector(max_m * max_n + 1);

    // compute non-trivial distributions
    for(int i = 1; i <= max_m; i++) {
      int max_j = std::min<int>(i, n_unique[pos_n_unique]);
      for(int j = 1; j <= max_j; j++) {
        if(j == 1) {
          for(int k = 0; k <= i; k++) dist[0][k] = 1.0/(i + 1);
        } else {
          int end1 = i * (j - 1);
          int end2 = i * j;

          if(i == j) {
            for(int k = end2; k >= j; k--)
              dist[j - 1][k] = dist[j - 2][k - j] * i;
          } else {
            for(int k = end2; k >= j; k--)
              dist[j - 1][k] = dist[j - 1][k - j] * i;
          }
          for(int k = 0; k < j; k++)
            dist[j - 1][k] = dist[j - 2][k] * j / (i + j);
          for(int k = j; k <= end1; k++)
            dist[j - 1][k] += dist[j - 2][k] * j;
          for(int k = j; k <= end2; k++)
            dist[j - 1][k] /= i + j;
        }

        while(pos_pair < len && m[pos_pair] == i && n[pos_pair] == j) {
          out[pos_pair] = dist[j - 1][Range(0, i * j)];
          pos_pair++;
        }
      }

      if(pos_n_unique < len_unique && m_unique[pos_n_unique] == i && n_unique[pos_n_unique] == max_j) {
        pos_n_unique++;
      }
    }
  }

  // return results
  return out;
}
