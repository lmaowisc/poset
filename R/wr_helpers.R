
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


  return(list(beta = as.vector(beta), Sigma = Sigma, err = err, l = l))
}

