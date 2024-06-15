


# Two-sample win ratio ----------------------------------------------------

#' Two-sample win ratio (net benefit) analysis
#'
#' @description Estimate and make inference on win ratio (net benefit) comparing a treatment to a control group.
#'
#' @param Y1 \eqn{K}-variate response data on \eqn{n_1} subjects in treatment (\eqn{n_1\times K} matrix).
#' @param Y0 \eqn{K}-variate response data on \eqn{n_0} subjects in control (\eqn{n_0\times K} matrix).
#' @param fun User-specified win function for pairwise comparison.
#'      It takes two arguments \eqn{y_1}
#'      and \eqn{y_0} (both \eqn{p}-dimensional) and returns 1 if \eqn{y_1} wins,
#'      -1 if \eqn{y_0} wins, and 0 if tied. The default is \code{\link{wprod}}
#'      based on the product order of multivariate ordinal data.
#' @return An object of class \code{wrtest} with the following components:
#' \item{theta}{A bivariate vector of win/loss fractions by LWR.}
#' \item{lgwr, lgwr_se, lgwr_pval}{Log-win ratio estimate, standard error, and p-value.}
#' \item{nb, nb_se, nb_pval}{Net benefit estimate (\code{theta[1]-theta[2]}), standard error, and p-value.}
#' @importFrom stats pnorm
#' @references
#' Mao, L. (2024). Win ratio for partially ordered data.
#' \emph{Statistica Sinica}, Under revision.
#'
#'  Buyse, M. (2010).  \href{https://doi.org/10.1002/sim.3923}{Generalized pairwise
#'  comparisons of prioritized outcomes in the two‐sample problem.}
#'  \emph{Statistics in Medicine}, 29, 3245-3257.
#' @keywords wrtest
#' @importFrom stats pnorm
#' @importFrom stats complete.cases
#' @export
#' @aliases wrtest
#' @seealso \code{\link{print.wrtest}}.
#' @examples
#' head(liver)
#' ### data example
#' ## compare bivariate fibrosis ratings between fibrosis
#' ## stages >= 3 and < 3
#'
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
#' include point estimates, 95\% confidence intervals, and p-values.
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
  cat("Win ratio (95% CI): ", round(exp(lgwr), 1), " (",
      round(exp(lgwr - za * lgwr_se), 1), ", ",
      round(exp(lgwr + za * lgwr_se), 1), "), p-value = ",
      lgwr_pval, "\n", sep = "")
  cat("Net benefit (95% CI): ", round(nb, 3), " (",
      round(nb - za * nb_se, 3), ", ",
      round(nb + za * nb_se, 3), "), p-value = ",
      nb_pval, "\n\n", sep = "")

  # print(desc)
}

#' The product-order win function for multivariate ordinal data
#'
#' @description A commonly used rule to compare
#'    multivariate ordinal data in \code{\link{wrtest}}. A winner has all its components
#'    greater than or equal to those of the loser, and strictly
#'    so for at least one component.
#'
#'
#' @param y1 A \eqn{K}-dimensional vector \eqn{y_1}.
#' @param y2 A \eqn{K}-dimensional vector  \eqn{y_2}.
#' @return An integer in \eqn{{1, 0, -1}}:
#'         \item{1}{If \eqn{y_1 \ge y_2} component-wise, with strict inequality for at least
#'                  one component.}
#'         \item{-1}{If \eqn{y_0 \ge y_2} component-wise, with strict inequality for at least
#'                  one component.}
#'         \item{0}{Otherwise.}
#'
#' @keywords wrtest
#' @export
#' @seealso \code{\link{wrtest}}, \code{\link{print.wrtest}}.
wprod <- function(y1, y2) {
  return(all(y1 >= y2) * any(y1 > y2) - all(y1 <= y2) * any(y1 < y2))
}




# Win ratio regression ----------------------------------------------------

#' Win ratio regression analysis
#'
#' @description Fit a multiplicative win-ratio regression model to
#' partially ordered outcome against covariates.
#' @param Y An \eqn{n\times K} matrix for \eqn{K}-variate response data on \eqn{n} subjects.
#'  The entries must be numeric scores.
#'  For pseudo-efficient estimation with specifying score function \code{sfun},
#'  the average score across components (row means)
#' should be compatible with the partial order (i.e., a monotone function on ordered elements).
#' @param Z An \eqn{n\times p} design matrix for covariates.
#' @param fun User-specified win function for pairwise comparison.
#'      It takes two arguments \eqn{y_1}
#'      and \eqn{y_0} (both \eqn{p}-dimensional) and returns 1 if \eqn{y_1} wins,
#'      -1 if \eqn{y_0} wins, and 0 if tied. The default is \code{\link{wprod}}
#'      based on the product order of multivariate ordinal data.
#' @param sfun The score function used in pseudo-efficient estimation.
#'    The default is the mean of component-wise numeric scores in \code{Y}.
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

#' @seealso \code{\link{wprod}}
#' @export
#' @importFrom utils combn
#' @importFrom stats complete.cases
#' @aliases wreg
#' @keywords wreg
#' @references Mao, L. (2024). Win ratio for partially ordered data.
#' \emph{Statistica Sinica}, Under revision.
#' @examples
#' set.seed(12345)
#' n <- 200
#' Z <- cbind(rnorm(n),rnorm(n))
#' Y1 <- rbinom(n,1,exp(Z[,1])/(1+exp(Z[,1])))
#' obj1 <- wreg(Y1,Z)
#' obj1$beta
#' sqrt(diag(obj1$Sigma))
#' obj1$beta/sqrt(diag(obj1$Sigma))
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
#' @description Print concise results for win ratio regression analysis.
#'
#' @param x An object returned by \code{\link{wreg}}.
#' @param ... Further arguments passed to or from other methods
#' @return No return value, called for side effects.
#' @export
#' @seealso \code{\link{wrtest}}.
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
  colnames(beta_tab) <- colnames(Z)
  rownames(beta_tab) <- ""
  print(as.data.frame(beta_tab))
  # print(desc)
}



#' Summarize model results from \code{wreg}
#'
#' Summarize the inferential results for regression
#'
#'
#' @param x An object returned by \code{\link{wreg}}.
#' @param ...  Additional arguments affecting the summary produced.
#'
#' @return An object of class \code{summary.wreg} with components
#'
#' \item{coefficients}{A \eqn{(J\times 4)}-dimensional matrix containing the inference
#' results for the log-loss rate; Columns include
#' \code{Estimate}, \code{Std.Err}, \code{Z value}, and \code{Pr(>|z|)}.}
#' @seealso \code{\link{wreg}}.
#' @keywords wreg
#' @importFrom stats pchisq pnorm
#' @examples
#' #See examples for wreg().
#' @export
summary.wreg=function(x, ...){

  beta <- x$beta
  var <- x$var

  cat("Call:\n")
  print(x$call)
  cat("\n")
  ## number of subjects
  cat("n = ", x$n, " subjects with complete data\n", sep ="")
  cat("Comparable (win/loss) pairs: ", x$Nwl, "/", x$N, " = ",
      round(100 * x$Nwl / x$N, 1), "%\n\n", sep ="")

  # create summary table
  p <- length(beta)

  coefficients <- matrix(NA,J,4)


  Dtab=matrix(NA,J,4)

  # rownames(LRtab)=c(paste0("Ref (Group ",groups[1],")"),
                    # paste0("Group ",groups[2:J]," vs ",ref))
  # colnames(LRtab)=c("Estimate", "Std.Err", "Z value", "Pr(>|z|)")




  class(result)<-"summary.wreg"
  return(result)



}





#' Print method for summary.LRfit objects
#'
#' Produces a printed summary of the results for the while-alive loss rate
#'
#' @param x An object returned by \code{\link{summary.LRfit}}.
#' @param ... Further arguments passed to or from other methods
#' @return No return value, called for side effects.
#' @export
#' @importFrom stats printCoefmat
print.summary.wreg=function(x,...){

  cat("Call:\n")
  print(x$call)
  cat("\n")
  tau=x$tau
  joint.test=x$joint.test
  J<-x$J
  # p-value for the (J-1)-d.f. test on loss rate
  LRchisq=x$LRchisq
  LRpval=x$LRpval

  # table for loss rate
  LRtab<-x$LRtab
  cat("Analysis of log loss rate (LR) by tau = ",tau, ":\n",sep="")
  printCoefmat(LRtab, P.values=TRUE, has.Pvalue=TRUE)
  cat("\n")

  cat("Test of group difference in while-alive LR\n")
  cat("X-squared = ", LRchisq, ", df = ",J-1,", p = ",LRpval, sep="")

  ## exponentiated table for loss rate ratio
  za<-qnorm(0.975)
  beta<-LRtab[2:J,1]
  se<-LRtab[2:J,2]
  LRR<-cbind(exp(beta),exp(beta-za*se),exp(beta+za*se))
  colnames(LRR)<-c("LR ratio","95% lower CL","95% higher CL")
  rownames(LRR)<-rownames(LRtab)[2:J]
  cat("\n\n")
  cat("Point and interval estimates for the LR ratio:\n")
  print(LRR)
  if (joint.test){
    # table for RMST
    Dtab<-x$Dtab
    LRDchisq<-x$LRDchisq
    LRDpval<-x$LRDpval
    cat("\n\n")
    cat("Analysis of log RMST (restricted mean survival time) by tau = ",tau, ":\n",sep="")
    printCoefmat(Dtab, P.values=TRUE, has.Pvalue=TRUE)
    cat("\n\n")
    cat("Test of group difference in while-alive LR and RMST\n")
    cat("X-squared = ", LRDchisq, ", df = ",2*(J-1),", p = ",LRDpval, sep="")

  }




}

