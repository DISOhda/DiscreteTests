% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class-results.R
\name{DiscreteTestResults}
\alias{DiscreteTestResults}
\title{Discrete Test Results Class}
\description{
This is the class used by \code{DiscreteTests} to output more detailed
results, if the \code{simple.output} parameter is set to \code{FALSE}. All
data members are private to prevent causing inconsistencies by deliberate or
inadvertent changes. However, the results can be read by public methods.
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-DiscreteTestResults-new}{\code{DiscreteTestResults$new()}}
\item \href{#method-DiscreteTestResults-get_pvalues}{\code{DiscreteTestResults$get_pvalues()}}
\item \href{#method-DiscreteTestResults-get_inputs}{\code{DiscreteTestResults$get_inputs()}}
\item \href{#method-DiscreteTestResults-get_scenario_supports}{\code{DiscreteTestResults$get_scenario_supports()}}
\item \href{#method-DiscreteTestResults-get_scenario_indices}{\code{DiscreteTestResults$get_scenario_indices()}}
\item \href{#method-DiscreteTestResults-summary}{\code{DiscreteTestResults$summary()}}
\item \href{#method-DiscreteTestResults-print}{\code{DiscreteTestResults$print()}}
\item \href{#method-DiscreteTestResults-clone}{\code{DiscreteTestResults$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-DiscreteTestResults-new"></a>}}
\if{latex}{\out{\hypertarget{method-DiscreteTestResults-new}{}}}
\subsection{Method \code{new()}}{
Creates a new \code{DiscreteTestResults} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DiscreteTestResults$new(
  test_name,
  inputs,
  alternative,
  p_values,
  scenario_supports,
  scenario_indices,
  data_name
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{test_name}}{a single character string with the name of the
test(s).}

\item{\code{inputs}}{a named list of \strong{exactly two} elements
containing the tests parameters. The first
element holds the observed data, e.g. a vector
or a matrix. The second one is a named list of
the \strong{unique} parameter combinations of the
tests, given by vectors.}

\item{\code{alternative}}{a single character string with the testing
alternative, e.g. "two.sided".}

\item{\code{p_values}}{a numeric vector of the p-values calculated by
each discrete testing scenario.}

\item{\code{scenario_supports}}{a list of \strong{unique} numeric vectors
containing all p-values the respective
discrete test setting can produce.}

\item{\code{scenario_indices}}{a list of numeric vectors containing the test
indices that indicates to which individual
testing scenario each unique parameter set and
each unique support belongs.}

\item{\code{data_name}}{a single character string with the name of the
variable that contains the observed data.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-DiscreteTestResults-get_pvalues"></a>}}
\if{latex}{\out{\hypertarget{method-DiscreteTestResults-get_pvalues}{}}}
\subsection{Method \code{get_pvalues()}}{
Returns the computed p-values.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DiscreteTestResults$get_pvalues()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A numeric vector of the p-values of each null hypothesis.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-DiscreteTestResults-get_inputs"></a>}}
\if{latex}{\out{\hypertarget{method-DiscreteTestResults-get_inputs}{}}}
\subsection{Method \code{get_inputs()}}{
Return the list of the test inputs. It can be chosen, if only unique
parameter sets are needed.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DiscreteTestResults$get_inputs(unique = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{unique}}{integer value that indicates whether only unique
parameter sets are to be returned. If \code{unique = FALSE}
(the default), the returned supports may be duplicated.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A list of two elements. The first one contains the observations for each
tested null hypothesis. The second is another list with the parameter
sets.  If \code{unique = TRUE}, only unique parameter sets are returned.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-DiscreteTestResults-get_scenario_supports"></a>}}
\if{latex}{\out{\hypertarget{method-DiscreteTestResults-get_scenario_supports}{}}}
\subsection{Method \code{get_scenario_supports()}}{
Returns the testing scenario supports. It can be chosen, if only unique
supports are needed.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DiscreteTestResults$get_scenario_supports(unique = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{unique}}{integer value that indicates whether only unique supports
are to be returned. If \code{unique = FALSE} (the default),
the returned supports may be duplicated.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A list of numeric vectors. Each one contains all observable p-values of
the respective null hypothesis. If \code{unique = TRUE}, only unique supports
are returned.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-DiscreteTestResults-get_scenario_indices"></a>}}
\if{latex}{\out{\hypertarget{method-DiscreteTestResults-get_scenario_indices}{}}}
\subsection{Method \code{get_scenario_indices()}}{
Returns the indices that indicate to which testing scenario each
unique support belongs.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DiscreteTestResults$get_scenario_indices()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A list of numeric vectors. Each one contains the indices of the null
hypotheses to which the respective support and/or parameter set belongs.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-DiscreteTestResults-summary"></a>}}
\if{latex}{\out{\hypertarget{method-DiscreteTestResults-summary}{}}}
\subsection{Method \code{summary()}}{
Returns the summary table which is printed by \code{print()}. For better
readability, the packages functions pass non-syntactic parameter names.
As a result, the returned summary table may have non syntactic column
names.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DiscreteTestResults$summary()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A data frame that lists all inputs and the resulting p-values. Each row
represents a testing scenario, i.e. the data of a single null hypothesis.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-DiscreteTestResults-print"></a>}}
\if{latex}{\out{\hypertarget{method-DiscreteTestResults-print}{}}}
\subsection{Method \code{print()}}{
Returns the computed p-values.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DiscreteTestResults$print(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{further arguments passed to \code{print.data.frame}.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Prints a summary of the tested null hypotheses. The object itself is
invisibly returned.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-DiscreteTestResults-clone"></a>}}
\if{latex}{\out{\hypertarget{method-DiscreteTestResults-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DiscreteTestResults$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
