

###########################################################################
#               General win-ratio analysis function                     ###
###########################################################################

#################### Two-sample win ratio######################################
# Input:
# dat1: n1 x p matrix (trt data)
# dat0: n0 x p matrix (control data)
# winfun: win function taking on two p-dimensional
# vectors y1, y0 and returning
#     1 if y1 wins
#    -1 if y0 wins
#     0 if indeterminate
#
# Output:
# log.WR: log-win ratio
# se: estimated standard error of log.WR
# wl: bivariate vector of win-loss probabilities
###########################################################################


#' Two-sample win ratio analysis
#'
#' @description Perform the win ratio (WR) analysis comparing a treatment to a control group.
#'
#' @param Y1 Treatment data on \eqn{n_1} subjects (\eqn{n_1\times K} matrix).
#' @param Y0 Control data on \eqn{n_0} subjects (\eqn{n_0\times K} matrix).
#' @param wf User-specified win function for pairwise comparison.
#'      It takes two arguments \eqn{y_1}
#'      and \eqn{y_0} (both \eqn{p}-dimensional) and returns 1 if \eqn{y_1} wins,
#'      -1 if \eqn{y_0} wins, and 0 if tied. The default is \code{\link{wprod}}
#'      based on the product order of multivariate ordinal data.
#' @return An object of class \code{WRtest} with the following components:
#' \item{theta}{A bivariate vector of win/loss fractions by LWR.}
#' \item{logWR,se}{Log-win ratio estimate and its standard error.}
#' \item{pval}{\eqn{p}-value of testing WR=1.}
#' @importFrom stats pnorm
#' @references Mao, L. and Wang, T. (2022+). Win ratio for partially ordered data.
#' Biometrika, In press.
#' @examples
#' # load the HF-ACTION trial data
#'
#' @keywords WRtest
#' @export
#' @aliases WRtest
#' @seealso \code{\link{print.WRtest}},
#' \code{\link{summary.WRtest}}.
WRtest <- function(Y1, Y0, wf = wprod) {
  dat1 <- as.matrix(Y1)
  dat0 <- as.matrix(Y0)

  n1 <- nrow(dat1)
  n0 <- nrow(dat0)

  # win-loss matrix
  # WL_ij=I(Y_1i\succ Y_0j)-I(Y_0j\succ Y_1i)
  WL <- matrix(0, n1, n0)

  for (i in 1:n1) {
    y1 <- dat1[i, ]
    for (j in 1:n0) {
      y0 <- dat0[j, ]
      r <- wf(y1, y0) # r=I(y1\succ y0)-I(y0\succ y1)
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
  Sigma <- n1^
    {
      -2
    } * t(w1) %*% w1 + n0^
      {
        -2
      } * t(w0) %*% w0

  # compute the log-win ratio and its
  # standard error


  # log-win ratio
  logWR <- log(theta[1] / theta[2])

  # derivative of log-WR
  f <- c(1 / theta[1], -1 / theta[2])

  # standard error of log-win ratio
  # by delta method
  se <- sqrt(t(f) %*% Sigma %*% f)
  pval <- 2 * (1 - pnorm(abs(logWR / se)))

  return(list(logWR = logWR, se = se, pval = pval, theta = theta))
}



#' The product-order win function for multivariate ordinal data
#'
#' @description A commonly used rule to compare
#'    multivariate ordinal data in \code{\link{WRtest}}. A winner has all its components
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
#' @keywords WRtest
#' @export
#' @aliases WRtest
#' @seealso \code{\link{WRtest}}, \code{\link{print.WRtest}},
#' \code{\link{summary.WRtest}}.
wprod <- function(y1, y2) {
  return(all(y1 >= y2) * any(y1 > y2) - all(y1 <= y2) * any(y1 < y2))
}






#################################################################################
##    Multiplicative win ratio regression model for multivariate ordinal data   #
#################################################################################


# dim(value)
# value[1:50,]

# n*(n-1)/2 x (1+p)


# (1, 2)
# (1, 3)
# ...
# (2 ,3)
# (2, 4)
# ...
# (n-1, n)

########################################
# For a matrix mat: n*(n-1)/2 x p
# whose rows are the p-vector valued
# anti-symmetric function f_ij
# in the order
# (1, 2)
# (1, 3)
# ...
# (2 ,3)
# (2, 4)
# ...
# (n-1, n)
#
# The following function extracts
# f_ij from mat
########################################
antifun_ij <- function(i, j, n, mat) {
  if (i < j) {
    k <- (i - 1) * (2 * n - i) / 2 + j - i
    return(mat[k, ])
  } else {
    if (i > j) {
      k <- (j - 1) * (2 * n - j) / 2 + i - j
      return(-mat[k, ])
    } else {
      return(rep(0, ncol(mat)))
    }
  }
}


###################
# Test antifun_ij #
###################
# antifun_ij(6,5,n,Zd)
# data[1:10,]




######################################
# Ef(x_i, X) for *symmetric* function f
# whose values are contained in the
# rows of mat
######################################

mean_i <- function(i, n, mat) {
  mat <- as.matrix(mat)
  k.list <- NULL
  if (i > 1) {
    j <- 1:(i - 1)
    k.list <- c(k.list, (j - 1) * (2 * n - j) / 2 + i - j)
  }
  if (i < n) {
    j <- (i + 1):n
    k.list <- c(k.list, (i - 1) * (2 * n - i) / 2 + j - i)
  }
  return(colSums(as.matrix(mat[k.list, ])) / (n - 1))
}

###############
# Test mean_i #
###############
# mean_i(i=9,n=n,mat=delta*Zd)
#
# tmp=0
# for (i in 1:n){
#   tmp=tmp+winp(data[9,1:2],data[i,1:2])*(data[9,3]-data[i,3])
# }
# tmp/(n-1)

###################################
# Newton-Raphson algorithm
###################################
NR.MWR <- function(beta, delta, Zd, Z, W, ep = 1e-8) {
  N <- nrow(Zd)
  n <- nrow(Z)
  i.list <- rep(NA, N)
  j.list <- rep(NA, N)
  for (i in 1:(n - 1)) {
    i.list[((i - 1) * (2 * n - i) / 2 + 1):(i * (2 * n - i - 1) / 2)] <- i
    j.list[((i - 1) * (2 * n - i) / 2 + 1):(i * (2 * n - i - 1) / 2)] <- (i + 1):n
  }


  p <- length(beta)

  err <- 1
  l <- 0

  maxiter <- 100
  # NR algorithm

  while (err > ep && l < maxiter) {
    eta <- exp(Z %*% beta)

    eta.i <- eta[i.list]
    eta.j <- eta[j.list]

    Wresid <- W * delta * ifelse(delta == 1, eta.j, ifelse(delta == -1, eta.i, 0))

    Uf.mat <- Zd * matrix(rep(Wresid, p), N, p)

    Uf <- colMeans(Uf.mat)

    wz <- abs(delta) * W / (1 / eta.i + 1 / eta.j)


    Zdw <- Zd * matrix(rep(wz, p), N, p)

    Omega <- t(Zdw) %*% Zd / N

    betap <- beta

    beta <- beta + solve(Omega) %*% Uf

    err <- sum(abs(beta - betap))
    l <- l + 1
  }

  ## variance estimation

  # Influence function

  if.fun <- function(i) {
    return(mean_i(i, n, Uf.mat))
  }

  psi.mat <- sapply(1:n, if.fun)
  if (is.vector(psi.mat)){
    psi.mat <- t(as.matrix(psi.mat))
  }
  # cat(dim(psi.mat),"\n")
  # dim(psi.mat)
  invOmega <- solve(Omega)
  Sigma <- 4 * invOmega %*% (psi.mat %*% t(psi.mat)) %*% invOmega / (n * (n - p))


  return(list(beta = beta, Sigma = Sigma, err = err, l = l))
}


#' Multiplicative win ratio (MWR) regression analysis
#'
#' @description Fit a multiplicative win ratio (MWR) regression model to
#' partially ordered outcome against covariates
#'
#' @param Y An \eqn{n\times K} matrix containing \eqn{n} observations
#' of a \eqn{K}-dimensional outcome variable
#' (whose components must be of the same variable type).
#' @param Z An \eqn{n\times p} design matrix for covariates.
#' @param wf User-specified win function for pairwise comparison.
#'      It takes two arguments \eqn{y_1}
#'      and \eqn{y_0} (both \eqn{K}-dimensional) and returns 1 if \eqn{y_1} wins,
#'      -1 if \eqn{y_0} wins, and 0 if tied. The default is \code{\link{wprod}}
#'      based on the product order of multivariate ordinal data.
#' @param sfun The score function used in pseudo-efficient estimation (can be bypassed).
#' @param ep Convergence criterion in Newton-Raphson algorithm. The default is 1e-6.
#' @return An object of class \code{MWRreg} with the following components:
#'
#' \item{beta}{A vector of estimated regression coefficients.}
#' \item{Sigma}{Estimated covariance matrix for \code{beta}}
#' \item{l}{Number of Newton-Raphson iterations.}
#' \item{beta_nv}{Naive (non-pseudo-efficient) estimates of \code{beta}.}
#' \item{se_nv}{Estimated standard errors for \code{beta_nv}.}
#' \item{n}{Sample size \eqn{n} of input data.}
#' \item{Nwl}{Number of comparable pairs (those with a winner and loser)
#' out of the \eqn{n(n+1)/2} possible ones.}

#' @seealso \code{\link{wprod}}
#' @export
#' @importFrom utils combn
#' @importFrom stats complete.cases
#' @aliases MWRreg
#' @keywords MWRreg
#' @references Mao, L. (2024). Win ratio for partially ordered data.
#' Statistica Sinica, Under revision.
#' @examples
#' set.seed(12345)
#' n <- 200
#' Z <- cbind(rnorm(n),rnorm(n))
#' Y1 <- rbinom(n,1,exp(Z[,1])/(1+exp(Z[,1])))
#' obj1 <- MWRreg(Y1,Z)
#' obj1$beta
#' sqrt(diag(obj1$Sigma))
#' obj1$beta/sqrt(diag(obj1$Sigma))

MWRreg <- function(Y, Z, wf = NULL, sfun = NULL, ep = 1e-6) {


  if (is.null(wf)) {
    wf <- wprod
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

    dij <- wf(yi, yj)

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

  Weff <- as.matrix(combn(eta, 2, function(x) {
    ((1 + x[1]) * (1 + x[2]))^{-1}
    }))

  # use efficient weights to refit data
  obj_eff <- NR.MWR(beta = beta_nv, delta, Zd, Z, Weff, ep)
  beta <- obj_eff$beta
  Sigma <- obj_eff$Sigma
  se <- sqrt(diag(Sigma))
  l <- obj_eff$l

  Nwl <- sum(delta != 0)



  result <- list(beta = beta, Sigma = Sigma, l = l, beta_nv = beta_nv, se_nv = se_nv, n = n, Nwl = Nwl)

  class(result) <- "MWRreg"

  return(result)
}
