% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DiscrimOD_assist.r
\name{getFEDInfo}
\alias{getFEDInfo}
\title{Generation function of Fedorov-Wynn algorithm parameter settings}
\usage{
getFEDInfo(
  FED_MAXIT = 200,
  FED_TRIM = 5,
  FED_TRIM_EPS = 0.001,
  freeRun = 1,
  FED_EPS = 1e-06,
  FED_ALPHA_GRID = 20
)
}
\arguments{
\item{FED_MAXIT}{A integer number of maximal iterations. The default is 200.}

\item{FED_TRIM}{A integer number that indicates the timing of design trimming (pruning).
The default is 5, meaning that, to trim the design every 5 iterations.}

\item{FED_TRIM_EPS}{A small value \eqn{\delta} with default \eqn{\delta=10^{-3}}.
In the trimming process, when the distance between two design points is smaller than \eqn{\delta}, the process will merge these two points.}

\item{freeRun}{A number between \eqn{[0,1]} that controls the percentage of updating iterations which are free from examining the
stopping criterion, \eqn{|f'-f|<\varepsilon}
where \eqn{f'} and \eqn{f} are the objective function values in the previous and current iterations, respectively.
The default is 1.0 implying the algorithm will completely ignore the stopping criterion.
Otherwise, the algorithm checks the stopping criterion after free iterations.}

\item{FED_ALPHA_GRID}{A integer number that specifies the grid size of candidate weights in the update-one-point process of the Fedorov-Wynn algorithm.
The default is 20.}

\item{tol}{A small value for the tolerance, \eqn{\varepsilon}, in the stopping criterion.
For \code{freeRun} smaller than 1.0, the default is \code{1e-6}. Otherwise, this value would not affect the algorithm.}
}
\value{
The list of Fedorov-Wynn parameter settings.
}
\description{
Create a list with Fedorov-Wynn algorithm parameters for optimal discrimination design search.
}
\examples{
# Get default settings for Fedorov-Wynn algorithm.
FED_INFO <- getFEDInfo(FED_MAXIT = 200)

}
