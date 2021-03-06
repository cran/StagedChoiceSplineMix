% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twostglogitregmixEM.r
\name{twostglogitregmixEM}
\alias{twostglogitregmixEM}
\title{Performs an EM algorithm for mixtures of two-stage logistic regressions}
\usage{
twostglogitregmixEM(yb = NULL, idb = NULL, yw = NULL, idw = NULL,
  xb = NULL, xw = NULL, xb.class = NULL, xw.class = NULL,
  sp.cols = NULL, num.knots = NULL, sp.knots = NULL, betab = NULL,
  betaw = NULL, lambda = NULL, k = NULL, epsilon = 1e-06, maxit = 500,
  maxrestarts = 100, maxer = 20, verb = FALSE)
}
\arguments{
\item{yb}{The 1st stage binary variable (browsing). See \code{\link{StagedChoiceSplineMix}} for details.}

\item{idb}{Corresponding id for yb. See \code{\link{StagedChoiceSplineMix}} for details.}

\item{yw}{The 2nd stage binary variable (writing). See \code{\link{StagedChoiceSplineMix}} for details.}

\item{idw}{Corresponding id for yw. See \code{\link{StagedChoiceSplineMix}} for details.}

\item{xb}{Corresponding covariance matrix for yb. See \code{\link{StagedChoiceSplineMix}} for details.}

\item{xw}{Corresponding covariance matrix for yw. See \code{\link{StagedChoiceSplineMix}} for details.}

\item{xb.class}{Corresponding latent classes for yb.}

\item{xw.class}{Corresponding latent classes for yw.}

\item{sp.cols}{See \code{\link{StagedChoiceSplineMix}} for details.}

\item{num.knots}{See \code{\link{StagedChoiceSplineMix}} for details.}

\item{sp.knots}{See \code{\link{StagedChoiceSplineMix}} for details.}

\item{betab}{See \code{\link{StagedChoiceSplineMix}} for details.}

\item{betaw}{See \code{\link{StagedChoiceSplineMix}} for details.}

\item{lambda}{See \code{\link{StagedChoiceSplineMix}} for details.}

\item{k}{See \code{\link{StagedChoiceSplineMix}} for details.}

\item{epsilon}{See \code{\link{StagedChoiceSplineMix}} for details.}

\item{maxit}{See \code{\link{StagedChoiceSplineMix}} for details.}

\item{maxrestarts}{See \code{\link{StagedChoiceSplineMix}} for details.}

\item{maxer}{See \code{\link{StagedChoiceSplineMix}} for details.}

\item{verb}{See \code{\link{StagedChoiceSplineMix}} for details.}
}
\description{
A sub-function of \code{\link{StagedChoiceSplineMix}}. This function performs an EM algorithm for mixtures of two-stage logistic regressions.
}
\seealso{
\code{\link{StagedChoiceSplineMix}} \cr
\code{\link[mixtools]{regmixEM}}
}

