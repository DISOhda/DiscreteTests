#include <array>
#include <Rcpp.h>
using namespace Rcpp;

// helper function for computing the means of the two samples
static inline std::array<double, 2> means_helper(
  const NumericVector& pool,
  const IntegerVector& chosen,
  const int& size_chosen,
  const int& size_unchosen,
  const double& sum_pool
) {
  //double sum_chosen = 0.0;
  //for(int i = 0; i < size_chosen; i++) sum_chosen += pool[chosen[i]];
  // Kahan summation to prevent rounding problems
  double sum_chosen = 0.0, error = 0.0;
  for (int i = 0; i < size_chosen; i++) {
    double y = pool[chosen[i]] - error;
    double t = sum_chosen + y;
    error = (t - sum_chosen) - y;
    sum_chosen = t;
  }
  return {sum_chosen/size_chosen, (sum_pool - sum_chosen)/size_unchosen};
}

// helper function for computing the "diff_mean" test statistic
static inline double diff_mean_stat(
  const NumericVector& pool,
  const IntegerVector& chosen,
  const int& size_chosen,
  const int& size_unchosen,
  const float& sign,
  const double& sum_pool
) {
  std::array<double, 2> means = means_helper(
    pool, chosen, size_chosen, size_unchosen, sum_pool
  );
  return sign * (means[0] - means[1]);
}

// helper function to obtain indices of values that were NOT selected
static inline IntegerVector complement_index(
  const IntegerVector& chosen,
  const int& size_chosen,
  const int& size_unchosen
) {
  IntegerVector unchosen(size_unchosen);
  // int j = 0, k = 0;
  // for(int i = 0; i < size_chosen + size_unchosen; i++) {
  //   if(j < size_chosen && i == chosen[j]) j++; else unchosen[k++] = i;
  // }
  std::vector<bool> use(size_chosen + size_unchosen, true);

  for(int i = 0; i < size_chosen; i++) use[chosen[i]] = false;

  int k = 0;
  for(int j = 0; j < size_chosen + size_unchosen; j++)
    if(use[j]) unchosen[k++] = j;

  return unchosen;
}

// helper function for computing the "diff_median" test statistic
static inline double diff_median_stat(
  const NumericVector& pool,
  const IntegerVector& chosen,
  const int& size_chosen,
  const int& size_unchosen,
  const float& sign
) {
  IntegerVector unchosen = complement_index(chosen, size_chosen, size_unchosen);
  NumericVector a = pool[chosen];
  NumericVector b = pool[unchosen];
  return sign * (median(a) - median(b));
}

// helper function for computing the "diff_hl" test statistic
static inline double diff_hl_stat(
  const NumericVector& pool,
  const IntegerVector& chosen,
  const int& size_chosen,
  const int& size_unchosen,
  const float& sign,
  const NumericMatrix& diffs
) {
  IntegerVector unchosen = complement_index(chosen, size_chosen, size_unchosen);
  NumericVector diffs_sel(size_chosen * size_unchosen);
  for(int i = 0; i < size_chosen; i++)
    for(int j = 0; j < size_unchosen; j++)
      diffs_sel[i * size_unchosen + j] = sign * diffs(chosen[i], unchosen[j]);
  return median(diffs_sel);
}

// helper function for computing the square sums with known intermediate sums
static inline std::array<double, 4> sqsums_helper(
  const NumericVector& pool,
  const IntegerVector& chosen,
  const int& size_chosen,
  const int& size_unchosen,
  const double& sum_pool,
  const double& sq_sum_pool
) {
  // object for output
  std::array<double, 4> out = {0.0, 0.0, 0.0, 0.0};

  // Kahan summation to prevent rounding problems
  double error_sq = 0.0, error_lin = 0.0;
  for (int i = 0; i < size_chosen; i++) {
    double v = pool[chosen[i]];

    double y = v * v - error_sq;
    double t = out[0] + y;
    error_sq = (t - out[0]) - y;
    out[0] = t;

    y = v - error_lin;
    t = out[2] + y;
    error_lin = (t - out[2]) - y;
    out[2] = t;
  }
  out[1] = sq_sum_pool - out[0];
  out[3] = sum_pool    - out[2];

  // compute variance square sums
  out[0] -= out[2] * out[2] / size_chosen;
  out[1] -= out[3] * out[3] / size_unchosen;

  // compute means
  out[2] /= size_chosen;
  out[3] /= size_unchosen;

  // return results
  return out;
}

// helper function for computing the "diff_t" test statistic
static inline double diff_t_stat(
  const NumericVector& pool,
  const IntegerVector& chosen,
  const int& size_chosen,
  const int& size_unchosen,
  const float& sign,
  const double& sum_pool,
  const double& sq_sum_pool
) {
  // compute means and square sums
  std::array<double, 4> sqs_means = sqsums_helper(
    pool, chosen, size_chosen, size_unchosen, sum_pool, sq_sum_pool
  );

  return sign * (sqs_means[2] - sqs_means[3]) /
    std::sqrt(1.0/size_chosen + 1.0/size_unchosen) /
    std::sqrt(
      (sqs_means[0] + sqs_means[1]) / (size_chosen + size_unchosen - 2)
    );
}

// helper function for computing the "diff_welch" test statistic
static inline double diff_welch_stat(
    const NumericVector& pool,
    const IntegerVector& chosen,
    const int& size_chosen,
    const int& size_unchosen,
    const float& sign,
    const double& sum_pool,
    const double& sq_sum_pool
) {
  // compute means and square sums
  std::array<double, 4> sqs_means = sqsums_helper(
    pool, chosen, size_chosen, size_unchosen, sum_pool, sq_sum_pool
  );

  return sign * (sqs_means[2] - sqs_means[3]) /
    std::sqrt(
      sqs_means[0]/size_chosen/(size_chosen - 1) +
        sqs_means[1]/size_unchosen/(size_unchosen - 1)
    );
}

// helper function for computing the "ratio_var" and "ratio_sd" test statistics
static inline double ratio_stats(
  const NumericVector& pool,
  const IntegerVector& chosen,
  const int& size_chosen,
  const int& size_unchosen,
  const float& sign,
  const double& sum_pool,
  const double& sq_sum_pool,
  const bool sd = false
) {
  // compute means and square sums
  std::array<double, 4> sqs_means = sqsums_helper(
    pool, chosen, size_chosen, size_unchosen, sum_pool, sq_sum_pool
  );

  // compute statistic
  double stat = std::pow(
    sqs_means[0]/(size_chosen - 1) * (size_unchosen - 1)/sqs_means[1],
    sign
  );
  if(sd) stat = std::sqrt(stat);

  // return result
  return stat;
}

static inline double compute_stat(
  const NumericVector& pool,
  const IntegerVector& chosen,
  const int& size_chosen,
  const int& size_unchosen,
  const double& sign,
  const double& sum_pool,
  const double& sq_sum_pool,
  const NumericMatrix& diffs,
  const std::string& method
) {
  double stat = 0;
  if(method == "diff_median")
    stat = diff_median_stat(pool, chosen, size_chosen, size_unchosen, sign);
  else if(method == "diff_hl")
    stat = diff_hl_stat(pool, chosen, size_chosen, size_unchosen, sign, diffs);
  else if(method == "diff_mean")
    stat = diff_mean_stat(
      pool, chosen, size_chosen, size_unchosen, sign, sum_pool
    );
  else if(method == "diff_t")
    stat = diff_t_stat(
      pool, chosen, size_chosen, size_unchosen, sign, sum_pool, sq_sum_pool
    );
  else if(method == "diff_welch")
    stat = diff_welch_stat(
      pool, chosen, size_chosen, size_unchosen, sign, sum_pool, sq_sum_pool
    );
  else
    stat = ratio_stats(
      pool, chosen, size_chosen, size_unchosen, sign,
      sum_pool, sq_sum_pool, method == "ratio_sd"
    );

  return stat;
}

// [[Rcpp::export]]
List perm_test_run(
  const NumericVector x,
  const NumericVector y,
  const double mu,
  const std::string method = "diff_mean",
  const bool exact = true,
  const int sims = 100000
) {
  // input sizes
  int size_x = x.size();
  int size_y = y.size();

  // combined observations
  NumericVector pool = NumericVector(size_x + size_y);
  int size_pool = size_x + size_y;
  if(method == "ratio_sd")
    pool[Range(0, size_x - 1)] = x / mu;
  else if(method == "ratio_var")
    pool[Range(0, size_x - 1)] = x / std::sqrt(mu);
  else pool[Range(0, size_x - 1)] = x - mu;
  pool[Range(size_x, size_pool - 1)] = y;

  // determine minimum and maximum size
  int size_min, size_max;
  float sign;
  if(size_x >= size_y) {
    size_min = size_y;
    size_max = size_x;
    sign = -1.0;
  } else {
    size_min = size_x;
    size_max = size_y;
    sign = 1.0;
  }

  // compute all differences or the sum of observations, depending on method
  double observed, sum_pool, sq_sum_pool;
  NumericMatrix diffs_pool(size_pool, size_pool);
  if(method == "diff_median") {
    observed = median(x) - median(y) - mu;
  } else if(method == "diff_hl") {
    for(int i = 0; i < size_pool; i++)
      for(int j = 0; j < size_pool; j++)
        diffs_pool(i, j) = pool[i] - pool[j];

    NumericMatrix sel = diffs_pool(
      Range(0, size_x - 1), Range(size_x, size_pool - 1)
    );
    observed = median(sel);
  } else if(
    method == "diff_mean"  ||
    method == "diff_t"     ||
    method == "diff_welch" ||
    method == "ratio_var"  ||
    method == "ratio_sd"
  ) {
    sum_pool = sum(pool);
    double sum_x = sum(x);
    observed = sum_x/size_x - (sum_pool + size_x * mu - sum_x)/size_y - mu;
    if(method != "diff_mean") {
      sq_sum_pool = sum(pow(pool, 2.0));
      double vx = var(x), vy = var(y);
      if(method == "diff_t") {
        double sd_pooled = std::sqrt(
          ((size_x - 1) * vx + (size_y - 1) * vy) / (size_pool - 2)
        );
        observed = observed / sd_pooled / std::sqrt(1.0/size_x + 1.0/size_y);
      } else if(method == "diff_welch") {
        observed = observed / std::sqrt(vx/size_x + vy/size_y);
      } else {
        observed = vx / mu / vy;
        if(method == "ratio_sd")
          observed = std::sqrt(observed / mu);
      }
    }
  }

  // statistics map
  std::map<double, int> statistics;

  int size_stats;
  if(exact) {
    // determine number of all possible combinations
    size_stats = 1;
    for(int i = 0; i < size_min; i++) {
      size_stats *= size_pool - i;
      size_stats /= i + 1;
    }

    // determine statistics via computing all combinations
    IntegerVector chosen(size_min, 0);
    std::vector<int> next(size_min, 0);
    int depth = 0;

    while(depth >= 0) {
      checkUserInterrupt();
      if(depth == size_min) {
        // compute statistic
        double stat = compute_stat(
          pool, chosen, size_min, size_max, sign,
          sum_pool, sq_sum_pool, diffs_pool, method
        );

        // count statistic
        ++statistics[stat];

        --depth;
        continue;
      }

      if(next[depth] < size_pool) {
        chosen[depth] = next[depth]++;
        ++depth;
        if(depth < size_min) next[depth] = chosen[depth - 1] + 1;
      } else {
        --depth;
      }
    }
  } else {
    size_stats = sims;

    // count observed statistic
    ++statistics[observed];

    // run simulations
    for(int i = 1; i < sims; i++) {
      // random index set of chosen indices
      IntegerVector chosen = sample(size_pool, size_min, false) - 1;
      //std::sort(chosen.begin(), chosen.end());

      // compute statistic
      double stat = compute_stat(
        pool, chosen, size_min, size_max, sign,
        sum_pool, sq_sum_pool, diffs_pool, method
      );

      // count statistic
      ++statistics[stat];
    }
  }

  // extract sums and frequencies
  NumericVector sums(statistics.size());
  IntegerVector freqs(statistics.size());
  NumericVector percs(statistics.size());
  int i = 0;
  for(const auto& [sum, freq] : statistics) {
    sums[i] = sum;
    //freqs[i] = freq;
    percs[i] = (double)freq/size_stats;
    ++i;
  }

  // return output
  return List::create(
    Named("observed") = observed,
    Named("statistics") = sums,
    //Named("frequencies") = freqs,
    Named("probabilities") = percs
  );
}
