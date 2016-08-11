#' @importFrom stats glm binomial rnorm runif
#' @title Generates initial values for mixtures of logistic regressions
#' @description A sub-function of \code{\link{StagedChoiceSplineMix}}. This function generates initial values for mixtures of logistic regressions. This function is used if starting points for parameters are not specified by the user or when the EM algorithm needs to be initialized due to errors.
#' @param y See \code{\link{StagedChoiceSplineMix}} for details.
#' @param x See \code{\link{StagedChoiceSplineMix}} for details.
#' @param k See \code{\link{StagedChoiceSplineMix}} for details.
#' @param er The total number of errors. See \code{\link{StagedChoiceSplineMix}} for details.
#' @seealso
#' \code{\link{StagedChoiceSplineMix}} \cr
#' "mixtools" package version 1.0.3
#' @export
gen.init<-
function (y, x, k, er)
{
  p <- ncol(x[[1]])
  beta.hyp<-matrix(0, p, k)
  for (j in 1:k) {
    beta.hyp[, j] <- as.matrix(glm(y~ as.matrix(x[[j]])-1, family = binomial())$coef)
  }

  beta.hyp[is.na(beta.hyp)] <- 0
  sd.hyp <- (er+1)*0.2
  beta <- matrix(0, p, k)
  for (j in 1:k) {
    beta[, j] <- rnorm(p, mean = as.vector(beta.hyp[, j]), sd = sd.hyp)
  }

  lambda <- runif(k,0.001,0.999)
  lambda <- lambda/sum(lambda)

  list(lambda = lambda, beta = beta, k = k)
}
