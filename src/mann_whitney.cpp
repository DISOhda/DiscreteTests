#include <Rcpp.h>
using namespace Rcpp;

// sort order
// IntegerVector order(const NumericVector& x, bool descending = false) {
//   arma::vec y = as<arma::vec>(x);
//   IntegerVector ord = as<IntegerVector>(wrap(arma::sort_index(y)));
//
//   if(descending)
//     return rev(ord);
//   else
//     return ord;
// }

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
    while(m[pos_pair] == 0 || n[pos_pair] == 0) {
      out[pos_pair] = NumericVector(1, 1.0);
      pos_pair++;
    }

    // find first unique non-zero second sample size
    int pos_n_unique = 0;
    while(n_unique[pos_n_unique] == 0) pos_n_unique++;

    // array of vectors of NumericVectors
    std::vector<NumericVector> dists[2] = {
      std::vector<NumericVector>(max_n),
      std::vector<NumericVector>(max_n)
    };

    // which array field contains previous and new distributions?
    int olds = 0, news = 1;

    for(int i = 1; i <= max_m; i++) {
      int max_j = std::min<int>(i, n_unique[pos_n_unique]);
      //return List::create(max_m, max_n, len, len_unique, pos_pair, pos_n_unique, max_j);
      for(int j = 1; j <= max_j; j++) {
        if(j == 1) {
          dists[news][j - 1] = NumericVector(i + 1, 1.0/(i + 1));
        } else {
          int end1 = i * (j - 1);
          int end2 = i * j;
          NumericVector dist(end2 + 1);

          for(int k = 0; k <= end1; k++)
            dist[k] = dists[news][j - 2][k] * j;

          for(int k = j; k <= end2; k++) {
            if(i == j)
              dist[k] += dists[news][j - 2][k - j] * i;
            else
              dist[k] += dists[olds][j - 1][k - j] * i;
          }

          dists[news][j - 1] = dist / (i + j);
        }

        while(pos_pair < len && m[pos_pair] == i && n[pos_pair] == j) {
          out[pos_pair] = dists[news][j - 1];
          pos_pair++;
        }
      }
      olds = 1 - olds;
      news = 1 - news;

      if(pos_n_unique < len_unique && m_unique[pos_n_unique] == i && n_unique[pos_n_unique] == max_j) {
        pos_n_unique++;
      }
    }
  }

  return out;
}
