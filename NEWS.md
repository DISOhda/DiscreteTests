## DiscreteTests 0.1.3

* Minor bug fix for `fisher.test.pv()` and  `mcnemar.test.pv()` functions:
  malformed input data was sometimes still accepted, but caused errors later


## DiscreteTests 0.1.2

* Minor improvements of `print()` method in `DiscreteTestResults` class
* Additions to examples in `DiscreteTestResultsSummary` class and
  `summary.DiscreteTestResults()` function


## DiscreteTests 0.1.1

* Improved `print()` method in `DiscreteTestResults` class
  - long lines are now wrapped and indented with `strwrap()`
  - new parameter `test_idx` for restricting print output to selected hypotheses
  - new parameter `limit` for limiting the number of hypotheses to be printed
    (default: 10); user is informed if limit is reached
* Fixed examples in documentation of `mcnemar.test.pv()`
* Updated references


## DiscreteTests 0.1.0

* Initial release.
* Added a `NEWS.md` file to track changes to the package.
* Added GitHub.
