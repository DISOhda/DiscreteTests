#' @return
#' If `simple.output = TRUE`, a vector of computed p-values is returned.
#' Otherwise, the output is a [`DiscreteTestResults`] R6 class object, which
#' also includes the p-value supports and testing parameters. These have to be
#' accessed by public methods, e.g. `$get_pvalues()`.
