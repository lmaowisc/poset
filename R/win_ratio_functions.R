


# Two-sample win ratio ----------------------------------------------------

#' Two-sample win ratio (net benefit) analysis
#'
#' @description Estimate and make inference on win ratio (net benefit) comparing a treatment to a control group.
#'
#' @param Y1 \eqn{K}-variate response data on \eqn{n_1} subjects in treatment (\eqn{n_1\times K} matrix).
#' @param Y0 \eqn{K}-variate response data on \eqn{n_0} subjects in control (\eqn{n_0\times K} matrix).
#' @param fun User-specified win function for pairwise comparison.
#'      It takes two arguments \eqn{y_1}
#'      and \eqn{y_0} (both \eqn{K}-vectors) and returns 1 if \eqn{y_1} wins,
#'      -1 if \eqn{y_0} wins, and 0 if tied. The default is \code{\link{wprod}}
#'      for the product order of multivariate ordinal data.
#' @return An object of class \code{wrtest} with the following components:
#' \item{theta}{A bivariate vector of win/loss fractions.}
#' \item{lgwr, lgwr_se, lgwr_pval}{Log-win ratio estimate (\code{log(theta[1]/theta[2])}), standard error, and p-value.}
#' \item{nb, nb_se, nb_pval}{Net benefit estimate (\code{theta[1]-theta[2]}), standard error, and p-value.}
#' @importFrom stats pnorm
#' @references
#' Mao, L. (2024). Win ratio for partially ordered data.
#' \emph{Statistica Sinica}, Under revision.
#'
#'  Buyse, M. (2010).  \href{\doi{10.1002/sim.3923}}{Generalized pairwise
#'  comparisons of prioritized outcomes in the two-sample problem.}
#'  \emph{Statistics in Medicine}, 29, 3245-3257.
#' @keywords wrtest
#' @importFrom stats pnorm
#' @importFrom stats complete.cases
#' @export
#' @aliases wrtest
#' @seealso \code{\link{wprod}}, \code{\link{print.wrtest}}.
#' @examples
#' head(liver)
#' ## compare bivariate ratings by fibrosis stage
#' ## lower score is better
#' Y1 <- liver[liver$AF, c("R1NASH", "R2NASH")] # advanced
#' Y0 <- liver[!liver$AF, c("R1NASH", "R2NASH")] # not advanced
#' obj <- wrtest(Y1, Y0)
#' obj
wrtest <- function(Y1, Y0, fun = wprod) {
  # convert to matrix
  dat1 <- as.matrix(Y1)
  dat0 <- as.matrix(Y0)
  # remove missing data
  dat1 <- as.matrix(dat1[complete.cases(dat1), ])
  dat0 <- as.matrix(dat0[complete.cases(dat0), ])

  n1 <- nrow(dat1)
  n0 <- nrow(dat0)

  # win-loss matrix
  # WL_ij=I(Y_1i\succ Y_0j)-I(Y_0j\succ Y_1i)
  WL <- matrix(0, n1, n0)
  for (i in 1:n1) {
    y1 <- dat1[i, ]
    for (j in 1:n0) {
      y0 <- dat0[j, ]
      r <- fun(y1, y0) # r=I(y1\succ y0)-I(y0\succ y1)
      WL[i, j] <- r
    }
  }
  win.mat <- (WL == 1)
  loss.mat <- (WL == -1)
  # win-loss probabilities
  theta1 <- mean(win.mat)
  theta0 <- mean(loss.mat)
  theta <- c(theta1, theta0)

  # influence functions for win-loss probabilities
  # of treatment and control arms
  w1 <- cbind(rowMeans(win.mat) - theta1, rowMeans(loss.mat) - theta0) # n1 x 2
  w0 <- cbind(colMeans(win.mat) - theta1, colMeans(loss.mat) - theta0) # n0 x 2
  # colMeans(w1): should be c(0,0)

  # variance of the win-loss estimators
  Sigma <- n1^{-2} * t(w1) %*% w1 + n0^{-2} * t(w0) %*% w0

  # compute the log-win ratio (NB) and its
  # standard error

  # log-win ratio -------------------
  lgwr <- log(theta[1] / theta[2])
  # derivative of log-WR
  f <- c(1 / theta[1], - 1 / theta[2])
  # standard error of log-win ratio
  # by delta method
  lgwr_se <- sqrt(t(f) %*% Sigma %*% f)
  lgwr_pval <- 2 * (1 - pnorm(abs(lgwr / lgwr_se)))

  # net benefit  -------------------
  nb <- theta[1] - theta[2]
  # derivative of nb
  f <- c(1, - 1)
  # standard error of log-win ratio
  # by delta method
  nb_se <- sqrt(t(f) %*% Sigma %*% f)
  nb_pval <- 2 * (1 - pnorm(abs(nb / nb_se)))

  # combine results
  obj <- list(lgwr = lgwr, lgwr_se = lgwr_se, lgwr_pval = lgwr_pval,
              nb = nb, nb_se = nb_se, nb_pval = nb_pval,
              theta = theta, Sigma = Sigma, n1 = n1, n0 = n0)
  obj$call <- match.call()
  class(obj) <- "wrtest"
  return(obj)
}

#' Print results from \code{wrtest}
#'
#' @description Print the results for two-sample win ratio (net benefit) analysis,
#' including point estimates, 95\% confidence intervals, and p-values.
#'
#' @param x An object returned by \code{\link{wrtest}}.
#' @param ... Further arguments passed to or from other methods
#' @return No return value, called for side effects.
#' @importFrom stats qnorm
#' @export
#' @seealso \code{\link{wrtest}}.
print.wrtest=function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  # title
  n1 <- x$n1
  n0 <- x$n0
  N <- n1 * n0
  theta <- x$theta
  # number of pairs:
  cat("Two-sample (Y1 vs Y0) win ratio/net benefit analysis\n\n")
  cat("Number of pairs: N1 x N0 = ", n1, "x", n0, " = ", N, "\n")
  cat("  Win: ", round(N * theta[1]), " (", round(100 * theta[1], 1), "%)\n",
      "  Loss: ", round(N * theta[2]), " (", round(100 * theta[2], 1), "%)\n",
      "  Tie: ", round(N * ( 1 - theta[1] - theta[2])), " (", round(100 * ( 1 - theta[1] - theta[2]), 1), "%)\n",
      sep ="")
  cat("\n")
  # print out win ratio / net benefit results
  ## get numeric results
  lgwr <- x$lgwr
  lgwr_se <- x$lgwr_se
  lgwr_pval <- x$lgwr_pval
  nb <- x$nb
  nb_se <- x$nb_se
  nb_pval <- x$nb_pval
  za <- qnorm(0.975)
  cat("Win ratio (95% CI): ", round(exp(lgwr), 2), " (",
      round(exp(lgwr - za * lgwr_se), 2), ", ",
      round(exp(lgwr + za * lgwr_se), 2), "), p-value = ",
      lgwr_pval, "\n", sep = "")
  cat("Net benefit (95% CI): ", round(nb, 3), " (",
      round(nb - za * nb_se, 3), ", ",
      round(nb + za * nb_se, 3), "), p-value = ",
      nb_pval, "\n\n", sep = "")

  # print(desc)
}

#' The product-order win function for multivariate ordinal data
#'
#' @description A common rule of comparison for the \code{fun} argument
#'    in \code{\link{wrtest}} and \code{\link{wreg}}.
#'    A winner has all its components
#'    greater than or equal to those of the loser, and strictly
#'    so for at least one component.
#'
#'
#' @param y1 A \eqn{K}-dimensional vector \eqn{y_1}.
#' @param y0 A \eqn{K}-dimensional vector  \eqn{y_0}.
#' @return An integer in \eqn{{1, 0, -1}}:
#'         \item{1}{If \eqn{y_1 \ge y_0} component-wise, with strict inequality for at least
#'                  one component.}
#'         \item{-1}{If \eqn{y_0 \ge y_1} component-wise, with strict inequality for at least
#'                  one component.}
#'         \item{0}{Otherwise.}
#'
#' @keywords wrtest
#' @export
#' @seealso \code{\link{wrtest}}, \code{\link{wreg}}.
wprod <- function(y1, y0) {
  return(all(y1 >= y0) * any(y1 > y0) - all(y1 <= y0) * any(y1 < y0))
}




# Win ratio regression ----------------------------------------------------

#' Win ratio regression analysis
#'
#' @description Fit a multiplicative win-ratio regression model to
#' partially ordered response against covariates.
#' @param Y An \eqn{n\times K} matrix for \eqn{K}-variate response data on \eqn{n} subjects.
#'  The entries must be numeric.
#'  For pseudo-efficient estimation (without specifying \code{sfun}),
#'  the average score across components (row means)
#' should be compatible with the partial order (i.e., preserve the same order for any two
#' comparable and ordered elements).
#' @param Z An \eqn{n\times p} design matrix for covariates.
#' @param fun User-specified win function for pairwise comparison.
#'      It takes two arguments \eqn{y_1}
#'      and \eqn{y_0} (both \eqn{K}-vectors) and returns 1 if \eqn{y_1} wins,
#'      -1 if \eqn{y_0} wins, and 0 if tied. The default is \code{\link{wprod}}
#'      for the product order of multivariate ordinal data.
#' @param sfun The scoring function used in pseudo-efficient estimation.
#'    The default is to take the row means of \code{Y}.
#' @param ep Convergence criterion in Newton-Raphson algorithm. The default is 1e-6.
#' @return An object of class \code{wreg} with the following components:
#' \item{beta}{A vector of estimated regression coefficients.}
#' \item{var}{Estimated covariance matrix for \code{beta}}
#' \item{l}{Number of Newton-Raphson iterations.}
#' \item{beta_nv}{Naive (non-pseudo-efficient) estimates of \code{beta}.}
#' \item{se_nv}{Estimated standard errors for \code{beta_nv}.}
#' \item{n}{Sample size \eqn{n} of input data with non-missing values.}
#' \item{Nwl}{Number of comparable pairs (those with a win and loss)
#' out of the \eqn{n(n-1)/2} possible ones.}

#' @seealso \code{\link{wprod}}, \code{\link{print.wreg}}, \code{\link{summary.wreg}}.
#' @export
#' @importFrom utils combn
#' @importFrom stats complete.cases
#' @aliases wreg
#' @keywords wreg
#' @references Mao, L. (2024). Win ratio for partially ordered data.
#' \emph{Statistica Sinica}, Under revision.
#' @examples
#' head(liver)
#' # regress bivariate ratings against covariates
#' Y <- 5 - liver[, c("R1NASH", "R2NASH")] # lower score is better
#' Z <- cbind("Female" = liver$Sex == "F",
#'            liver[, c("AF", "Steatosis",   "SSF2",  "LSN")]) # covariates
#' obj <- wreg(Y, Z) # fit model
#' obj
#' summary(obj)

wreg <- function(Y, Z, fun = NULL, sfun = NULL, ep = 1e-6) {

  if (is.null(fun)) {
    fun <- wprod
    win = "mvprod"
  }


  # remove missing daya
  Y <- as.matrix(Y)
  Z <- as.matrix(Z)

  dat <- cbind(Y, Z)
  ind <- complete.cases(dat)

  Y <- as.matrix(Y[ind, ])
  Z <- as.matrix(Z[ind, ])

  # dimension of outcome
  K <- ncol(Y)
  # sample size
  n <- nrow(Y)
  N <- n * (n - 1) / 2 # number of pairs

  # p: number of covariates
  p <- ncol(Z)

  dat <- cbind(Y, Z)


  ####################################################
  # Take the ith and jth rows of the data matrix "dat"
  # and compute delta_ij and Z_i-Z_j
  # Argument: x=c(i,j)
  # implicit: data, K (response dimension)
  ###################################################
  contrast.ij <- function(x) {
    xi <- dat[x[1], ]
    xj <- dat[x[2], ]

    yi <- xi[1:K]
    yj <- xj[1:K]

    dij <- fun(yi, yj)

    l <- length(xi)
    Zdij <- xi[(K + 1):l] - xj[(K + 1):l]

    return(c(dij, Zdij))
  }
  value <- t(combn(1:n, 2, FUN = contrast.ij))


  # pairwise win-loss indicator
  delta <- as.matrix(value[, 1])
  # pairwise difference in covariates
  Zd <- as.matrix(value[, 2:(1 + p)])

  # naive weights
  W <- rep(1, n * (n - 1) / 2)
  #### Estimating function ###########################

  # Newton Raphson
  obj <- NR.MWR(beta = rep(0, p), delta, Zd, Z, W, ep=ep*10) # use a less stringent criteron

  beta_nv <- obj$beta
  se_nv <- sqrt(diag(obj$Sigma))



  #######################################################################
  #   Pseudo-efficient weight                   ##
  #######################################################################

  # input beta, r, and Z

  # compute the scores
  # sfun <- sum
  if (win == "mvprod") {
    mlev <- apply(Y, 2, max) # maximum level
    r <- rowMeans(Y / matrix(rep(mlev, each = n), n, K)) # normalized mean value
  } else {
    if (!is.null(sfun)) {
      r <- apply(Y, 1, sfun)
      r <- r / (max(r, na.rm = T) - min(r, na.rm = T))
    } else {
      r <- NULL
    }
  }

  ## computing gamma ###
  # via pseudo-logistic regression
  eta <- exp(Z %*% beta_nv)

  if (!is.null(r)) {
    # gamma.hat
    ghat <- 0
    err <- 1
    iter <- 0
    ep <- 1e-5
    maxiter <- 100
    # NR algorithm

    while (err > ep && iter < maxiter) {
      etag <- exp(ghat) * eta
      score <- mean(r - etag / (1 + etag))
      Info <- mean(etag / (1 + etag)^2)
      ghatp <- ghat
      ghat <- ghat + score / Info
      err <- abs(ghat - ghatp)
      iter <- iter + 1
    }
    # obtained ghat
    etag <- exp(ghat) * eta
  } else {
    etag <- eta
  }
  # compute efficient weights Weff
  # function dependent on ghat and eta

  Weff <- as.matrix(combn(etag, 2, function(x) {
    ((1 + x[1]) * (1 + x[2]))^{-1}
    }))

  # use efficient weights to refit data
  obj_eff <- NR.MWR(beta = beta_nv, delta, Zd, Z, Weff, ep)
  beta <- obj_eff$beta
  Sigma <- obj_eff$Sigma
  se <- sqrt(diag(Sigma))
  l <- obj_eff$l

  Nwl <- sum(delta != 0)

  obj <- list(beta = beta, var = Sigma,
                 l = l, beta_nv = beta_nv, se_nv = se_nv,
                 n = n, N = n * (n - 1) /2, Nwl = Nwl,
                  Y = Y, Z = Z)
  obj$call <- match.call()
  class(obj) <- "wreg"

  return(obj)
}

#' Print concise model results from \code{wreg}
#'
#' @description Print concise results for win ratio regression.
#'
#' @param x An object returned by \code{\link{wreg}}.
#' @param ... Further arguments passed to or from other methods
#' @return No return value, called for side effects.
#' @export
#' @seealso \code{\link{wreg}}.
print.wreg=function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  ## number of subjects
  cat("n = ", x$n, " subjects with complete data\n", sep ="")
  cat("Comparable (win/loss) pairs: ", x$Nwl, "/", x$N, " = ",
      round(100 * x$Nwl / x$N, 1), "%\n\n", sep ="")
  # print out beta
  beta_tab <- as.matrix(t(x$beta))
  colnames(beta_tab) <- colnames(x$Z)
  rownames(beta_tab) <- ""
  print(as.data.frame(beta_tab))
  # print(desc)
}



#' Summarize model results from \code{wreg}
#'
#' Summarize the inferential results for win ratio regression.
#'
#'
#' @param object An object returned by \code{\link{wreg}}.
#' @param ...  Additional arguments affecting the summary produced.
#'
#' @return An object of class \code{summary.wreg} with components:
#'
#' \item{coefficients}{A matrix of coefficients, standard errors, z-values and p-values.}
#' \item{exp_coef}{A matrix of win ratios (exp(coef)) and 95\% confidence intervals.}
#' \item{wald, wald_pval}{Overall wald test statistic on all covariates and p-value.}
#' @seealso \code{\link{wreg}}.
#' @keywords wreg
#' @importFrom stats qnorm pnorm pchisq
#' @examples
#' #See examples for wreg().
#' @export
summary.wreg=function(object, ...){

  x <- object
  beta <- x$beta
  var <- x$var
  se <- sqrt(diag(var))


  # create summary table
  p <- length(beta)
  coefficients <- matrix(NA, p, 5)
  rownames(coefficients) <- colnames(x$Z)
  colnames(coefficients)=c("coef",  "exp(coef)" , "se(coef)", "z", "Pr(>|z|)")
  ## fill in values
  coefficients[, 1] <- beta
  coefficients[, 2] <- exp(beta)
  coefficients[, 3] <- se
  coefficients[, 4] <- beta / se
  coefficients[, 5] <- 2 * (1 - pnorm(abs(coefficients[, 4])))
  # coefficients
  # printCoefmat(coefficients, P.values = TRUE, has.Pvalue = TRUE, digits = 4,
  #              cs.ind = c(1, 3), tst.ind = 4)

  ## exponentiated
  exp_coef <- matrix(NA, p, 4)
  rownames(exp_coef) <- rownames(coefficients)
  colnames(exp_coef) <- c("exp(coef)", "exp(-coef)", "lower .95", "upper .95")
  za <- qnorm(0.975)
  ## fill in values
  exp_coef[, 1] <- exp(beta)
  exp_coef[, 2] <- 1 / exp_coef[, 1]
  exp_coef[, 3] <- exp(beta - za * se)
  exp_coef[, 4] <- exp(beta + za * se)
  # exp_coef
  ## wald test
  wald <- as.numeric(t(beta) %*% solve(var) %*% beta)
  wald_pval <- 1 - pchisq(wald, p)
  # combine results
  result <- list(coefficients = coefficients, exp_coef = exp_coef, wald = wald,
                 wald_pval = wald_pval, p = p,
                 beta = beta, var = var, n = x$n, N = x$N, Nwl = x$Nwl,
                 Y = x$Y, Z = x$Z, l = x$l, call = x$call)
  class(result)<-"summary.wreg"
  return(result)



}


#' Print method for summary.wreg objects
#'
#' Print summary results for win ratio regression.
#'
#' @param x An object returned by \code{\link{summary.wreg}}.
#' @param ... Further arguments passed to or from other methods
#' @return No return value, called for side effects.
#' @export
#' @importFrom stats printCoefmat
print.summary.wreg <- function(x,...){

  cat("Call:\n")
  print(x$call)
  cat("\n")
  ## number of subjects
  cat("n = ", x$n, " subjects with complete data\n", sep ="")
  cat("Comparable (win/loss) pairs: ", x$Nwl, "/", x$N, " = ",
      round(100 * x$Nwl / x$N, 1), "%\n\n", sep = "")

  cat("Newton-Raphson algoritm converged in ", x$l, " iterations\n\n", sep = "")

  printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE, digits = 4,
                             cs.ind = c(1, 3), tst.ind = 4)
  cat("\n")
  printCoefmat(x$exp_coef)
  cat("\n")

  cat("Overall Wald test = ", round(x$wald, digits = 3), " on ",
      x$p, " df,  p = ", x$wald_pval, sep = "")
}


