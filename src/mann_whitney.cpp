#include <Rcpp.h>
using namespace Rcpp;

// Mann-Whitney probabilities for settings without ties
// [[Rcpp::export]]
List mann_whitney_probs_int(
  const IntegerVector m, // must be sorted!!
  const IntegerVector n  // each number corresponds to a value in 'm' and m >= n!!
) {
  // number of distributions to be computed
  uint32_t len = (uint32_t)m.length();
  // largest size of first samples
  uint32_t max_m = (uint32_t)m[len - 1];
  // largest size of second samples
  uint32_t max_n = (uint32_t)max(n);

  // list of results
  List out(len);

  // compute distributions
  if(max_n == 0) {
    // if all sample pairs contain zeros, all distributions are trivial
    for(uint32_t i = 0; i < len; i++) out[i] = NumericVector(1, 1.0);
  } else {
    // determine maximums of second sample sizes for efficiency
    IntegerVector m_unique;
    IntegerVector n_unique;
    uint32_t len_unique = 1;
    if(len > 1) {
      m_unique = IntegerVector(len);
      n_unique = IntegerVector(len);
      m_unique[0] = m[0];
      n_unique[0] = n[0];
      for(uint32_t i = 1; i < len; i++) {
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
    uint32_t pos_pair = 0;
    while(pos_pair < len && (m[pos_pair] == 0 || n[pos_pair] == 0)) {
      out[pos_pair] = NumericVector(1, 1.0);
      pos_pair++;
    }

    // find first unique non-zero second sample size
    uint32_t pos_n_unique = 0;
    while(pos_n_unique < len_unique && n_unique[pos_n_unique] == 0)
      pos_n_unique++;

    // array of NumericVectors for all distributions (1, ..., n) for m
    NumericVector* dist = new NumericVector[max_n];
    for(uint32_t i = 0; i < max_n; i++)
      dist[i] = NumericVector(max_m * max_n + 1);

    // compute non-trivial distributions
    for(uint32_t i = 1; i <= max_m; i++) {
      uint32_t max_j = std::min<uint32_t>(i, (uint32_t)n_unique[pos_n_unique]);
      for(uint32_t j = 1; j <= max_j; j++) {
        if(j == 1) {
          for(uint32_t k = 0; k <= i; k++) dist[0][k] = 1.0/(i + 1);
        } else {
          uint32_t end1 = i * (j - 1);
          uint32_t end2 = i * j;

          if(i == j) {
            for(uint32_t k = end2; k >= j; k--)
              dist[j - 1][k] = dist[j - 2][k - j] * i;
          } else {
            for(uint32_t k = end2; k >= j; k--)
              dist[j - 1][k] = dist[j - 1][k - j] * i;
          }
          for(uint32_t k = 0; k < j; k++)
            dist[j - 1][k] = dist[j - 2][k] * j / (i + j);
          for(uint32_t k = j; k <= end1; k++)
            dist[j - 1][k] += dist[j - 2][k] * j;
          for(uint32_t k = j; k <= end2; k++)
            dist[j - 1][k] /= i + j;
        }

        while(pos_pair < len && m[pos_pair] == i && n[pos_pair] == j) {
          out[pos_pair] = dist[j - 1][Range(0, i * j)];
          pos_pair++;
        }
      }

      if(
        pos_n_unique < len_unique &&
        m_unique[pos_n_unique] == i &&
        n_unique[pos_n_unique] == max_j
      ) {
        pos_n_unique++;
      }
    }

    // garbage collection
    delete[] dist;
  }

  // return results
  return out;
}

#include <Rcpp.h>
using namespace Rcpp;

// Mann-Whitney probabilities
// [[Rcpp::export]]
NumericVector mann_whitney_probs_ties_int(
  const IntegerVector ranks,
  const uint32_t size_x
) {
  IntegerVector scores = clone(ranks);
  std::sort(scores.begin(), scores.end());
  uint32_t len = (uint32_t)scores.size();
  uint32_t size_y = len - size_x;
  uint32_t size_l = std::min<uint32_t>(size_x, size_y);

  std::vector<uint64_t> min_scores(size_l + 1, 0), max_scores(size_l + 1, 0);
  min_scores[1] = scores[0], max_scores[1] = scores[len - 1];
  for(uint32_t i = 2; i <= size_l; i++) {
    min_scores[i] = min_scores[i - 1] + scores[i - 1];
    max_scores[i] = max_scores[i - 1] + scores[len - i];
  }
  uint64_t min_score = min_scores[size_l], max_score = max_scores[size_l];

  std::vector<std::vector<double>> d_list(size_l + 1);
  for(uint32_t i = 0; i <= size_l; i++) {
    d_list[i].assign(max_scores[i] - min_scores[i] + 1, 0.0);
  }
  d_list[0][0] = 1.0;

  std::vector<uint64_t> first(size_l + 1, 0), last(size_l + 1, 0);

  for(uint32_t i = 0; i < len; i++) {
    uint64_t s = scores[i];
    checkUserInterrupt();
    //Rcout << i << ":\n";
    for(uint32_t j = std::min(i + 1, size_l); j > 0; j--) {
      double w_keep = (double)(i + 1 - j) / (double)(i + 1);
      double w_add  = (double)(j) / (double)(i + 1);
      //Rcout << "  " << j - 1 << "\n";
      uint64_t src_lo = first[j - 1] - min_scores[j - 1];
      uint64_t src_hi =  last[j - 1] - min_scores[j - 1];
      uint64_t dst_lo = first[j - 1] - min_scores[j] + s;

      if(first[j]) {
        uint64_t old_lo = first[j] - min_scores[j];
        uint64_t old_hi = last[j]  - min_scores[j];
        for (uint64_t k = old_lo; k <= old_hi; k++)
          d_list[j][k] *= w_keep;
      }

      for(
        uint64_t k = src_lo, l = dst_lo; k <= src_hi; k++, l++
      )
        d_list[j][l] += d_list[j - 1][k] * w_add;

      if(!first[j]) {
        first[j] = first[j - 1] + s;
         last[j] =  last[j - 1] + s;
      } else {
        if(first[j - 1] + s < first[j]) first[j] = first[j - 1] + s;
        if( last[j - 1] + s >  last[j])  last[j] =  last[j - 1] + s;
      }
    }
  }

  uint32_t size_res = max_score - min_score + 1;

  double total = 0;
  for(uint64_t i = 0; i < size_res; i++)
    total += d_list[size_l][i];

  NumericVector res(size_res);
  for(uint64_t i = 0; i < size_res; i++)
    res[i] = d_list[size_l][size_x <= size_y ? i : size_res - i - 1] / total;

  return res;
}
