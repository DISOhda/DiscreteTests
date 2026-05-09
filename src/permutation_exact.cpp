#include <Rcpp.h>
using namespace Rcpp;

// static void permutation_backtrack(
//   NumericVector& set,
//   const int size_choices,
//   const int size_others,
//   int depth,
//   NumericVector& current,
//   std::vector<bool>& used,
//   std::map<double, int>& statistics,
//   const double set_sum,
//   const short sign
// ) {
//   checkUserInterrupt();
//
//   if(depth == size_choices) {
//     double current_sum = sum(current), others_sum = set_sum - current_sum;
//     double stat = sign * (current_sum/size_choices - others_sum/size_others);
//     ++statistics[stat];
//     return;
//   }
//
//   for(int i = 0; i < set.size(); i++) {
//     if(!used[i]) {
//       used[i] = true;
//       current[depth] = set[i];
//       permutation_backtrack(set, size_choices, size_others, depth + 1,
//                             current, used, statistics, set_sum, sign);
//       used[i] = false;
//     }
//   }
// }

IntegerVector complement_index(
  const IntegerVector& chosen,
  const int& size_min,
  const int& size_max,
  const int& size_xy
) {
  IntegerVector index(size_max);
  int j = 0, k = 0;
  for(int i = 0; i < size_xy; i++) {
    if(j < size_min && i == chosen[j]) j++; else index[k++] = i;
  }
  return index;
}

// [[Rcpp::export]]
List perm_test_all_combs(
  const NumericVector x,
  const NumericVector y,
  const double mu,
  const std::string method = "diff_mean"
) {
  // input sizes
  int size_x = x.size();
  int size_y = y.size();

  // combined observations
  NumericVector xy = NumericVector(size_x + size_y);
  int size_xy = size_x + size_y;
  if(method == "ratio_sd")
    xy[Range(0, size_x - 1)] = x / mu;
  else if(method == "ratio_var")
    xy[Range(0, size_x - 1)] = x / std::sqrt(mu);
  else xy[Range(0, size_x - 1)] = x - mu;
  xy[Range(size_x, size_xy - 1)] = y;

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

  // statistics map
  int size_stats = 1;
  for(int i = 0; i < size_min; i++) {
    size_stats *= size_xy - i;
    size_stats /= i + 1;
  }
  std::map<double, int> statistics;

  // determine statistics via computing all permutations
  IntegerVector chosen(size_min);
  std::vector<int> next(size_min, 0);

  // observed test statistic
  double observed, sum_xy;
  NumericMatrix diffs_xy(size_xy, size_xy);
  if(method == "diff_median") {
    //observed = median(x) - median(y) - mu;
  } else if(method == "diff_hl") {
    for(int i = 0; i < size_xy; i++)
      for(int j = 0; j < size_xy; j++)
        diffs_xy(i, j) = xy[i] - xy[j];

    //NumericMatrix sel = diffs_xy(
    //  Range(0, size_x - 1), Range(size_x, size_xy - 1)
    //);
    //observed = median(sel);
  } else if(
      method == "diff_mean" ||
      method == "diff_t"    ||
      method == "ratio_var" ||
      method == "ratio_sd"
  ) {
    // observed = mean(x) - mean(y) - mu;
    sum_xy = sum(xy);
    // if(method != "diff_mean") {
    //   double vx = var(x), vy = var(y);
    //   if(method == "diff_t") {
    //     double sd_pooled = std::sqrt(
    //       ((size_x - 1) * vx + (size_y - 1) * vy) / (size_xy - 2)
    //     );
    //     observed = observed / sd_pooled / std::sqrt(1.0/size_x + 1.0/size_y);
    //   } else {
    //     observed = vx / mu / vy;
    //     if(method == "ratio_sd")
    //       observed = std::sqrt(observed / mu);
    //   }
    // }
  } else
    stop("Unknown method '" + method + "'");

  int depth = 0;
  while(depth >= 0) {
    checkUserInterrupt();
    if(depth == size_min) {
      double stat = 0;
      if(method == "diff_median") {
        IntegerVector rest = complement_index(chosen, size_min, size_max, size_xy);
        NumericVector a = xy[chosen];
        NumericVector b = xy[rest];
        stat = sign * (median(a) - median(b));
      } else if(method == "diff_hl") {
        IntegerVector rest = complement_index(chosen, size_min, size_max, size_xy);
        NumericVector sel(size_min * size_max);
        for(int i = 0; i < size_min; i++)
          for(int j = 0; j < size_max; j++)
            sel[i * size_max + j] = sign * diffs_xy(chosen[i], rest[j]);
        stat = median(sel);
      } else { // method == "diff_mean"
        double current_sum = 0;
        for(int i = 0; i < size_min; i++) current_sum += xy[chosen[i]];
        double mx = current_sum/size_min, my = (sum_xy - current_sum)/size_max;
        stat = sign * (mx - my);
        if(method != "diff_mean") {
          IntegerVector rest = complement_index(chosen, size_min, size_max, size_xy);
          double ssx = 0, ssy = 0;
          for(int i = 0; i < size_min; i++)
            ssx += std::pow(xy[chosen[i]] - mx, 2.0);
          for(int i = 0; i < size_max; i++)
            ssy += std::pow(xy[rest[i]] - my, 2.0);
          if(method == "diff_t")
            stat = stat / std::sqrt((ssx + ssy) / (size_xy - 2)) /
              std::sqrt(1.0/size_min + 1.0/size_max);
          else {
            stat = std::pow(ssx * (size_max - 1) / ssy / (size_min - 1), sign);
            if(method == "ratio_sd") stat = std::sqrt(stat);
          }
        }
      }

      ++statistics[stat];

      --depth;
      continue;
    }

    if(next[depth] < size_xy) {
      chosen[depth] = next[depth]++;
      ++depth;
      if(depth < size_min) next[depth] = chosen[depth - 1] + 1;
    } else {
      --depth;
    }
  }

  // extract sums and frequencies
  NumericVector sums(statistics.size());
  IntegerVector freqs(statistics.size());
  NumericVector percs(statistics.size());
  int i = 0;
  for(const auto& [sum, freq] : statistics) {
    sums[i] = sum;
    freqs[i] = freq;
    percs[i] = (double)freq/size_stats;
    ++i;
  }

  // return output
  return List::create(
    //Named("observed") = observed,
    Named("statistics") = sums,
    //Named("frequencies") = freqs,
    Named("probabilities") = percs
  );
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
perm_test_all_combs(c(8, 4, 10), c(3, 5, 2), 0, "diff_mean")
*/
