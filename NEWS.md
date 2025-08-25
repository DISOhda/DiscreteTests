# DiscreteTests 0.3.0

* Added Wilcoxon's one- and two-sample sign rank tests (new functions:
  `wilcox_single_test_pv()` and `signed_rank_test_pv()`.
* Added Mann-Whitney U test (equivalent to Wilcoxon's rank sum test).
* Extensions to `DiscreteTestResults` class to store observations as lists (e.g.
  to store vectors of different sizes) as well as to include the tests'
  statistics, the name of the null distribution and whether the *p*-values were
  computed exactly.


# DiscreteTests 0.2.1

* `get_pvalues()` method of `DiscreteTestResults` class now includes the (row or
  cell) names of the input observations (if present). Output of names can be
  suppressed by a parameter.
* `print()` method of `DiscreteTestResults` class now also prints the (row or
  cell) names of the input observations (if present).


# DiscreteTests 0.2.0

* Migration to snake case for all function and parameter names. Previously
  existing functions are kept for compatibility, but are marked as deprecated
  and will be removed in a future version
* `print()` output of `DiscreteTestResults` objects now wraps output of all 
  input data, parameters and hypotheses correctly


# DiscreteTests 0.1.3

* Minor bug fix for `fisher.test.pv()` and  `mcnemar.test.pv()` functions:
  malformed input data was sometimes still accepted, but caused errors later


# DiscreteTests 0.1.2

* Minor improvements of `print()` method in `DiscreteTestResults` class
* Additions to examples in `DiscreteTestResultsSummary` class and
  `summary.DiscreteTestResults()` function


# DiscreteTests 0.1.1

* Improved `print()` method in `DiscreteTestResults` class
  - long lines are now wrapped and indented with `strwrap()`
  - new parameter `test_idx` for restricting print output to selected hypotheses
  - new parameter `limit` for limiting the number of hypotheses to be printed
    (default: 10); user is informed if limit is reached
* Fixed examples in documentation of `mcnemar.test.pv()`
* Updated references


# DiscreteTests 0.1.0

* Initial release.
* Added a `NEWS.md` file to track changes to the package.
* Added GitHub.
