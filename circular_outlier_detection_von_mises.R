# circular_outlier_detection_von_mises.R
# ------------------------------------------------------------
# Project: Circular Outlier Detection Using Mixture of Von Mises Distributions
# Author: Khadijeh Abdi
#
# This script is a GitHub-ready cleaned transcription of the R code
# from the scanned thesis appendix/program pages.
#
# IMPORTANT:
# The original source was a scanned PDF, so some symbols and variables
# should be manually checked against the original thesis before formal use.
# ------------------------------------------------------------

if (!requireNamespace("CircStats", quietly = TRUE)) {
  install.packages("CircStats")
}
library(CircStats)

# ------------------------------------------------------------
# 1. COVRATIO calculation
# ------------------------------------------------------------

calculate_covratio <- function(u, v, order = 1) {
  n <- length(u)
  data <- cbind(u, v)

  order.matrix <- t(matrix(rep(c(1:order), n), ncol = n))
  cos.u <- cos(u * order.matrix)
  sin.u <- sin(u * order.matrix)

  ones <- matrix(1, n, 1)
  U <- cbind(ones, cos.u, sin.u)

  V1 <- cos(v)
  V2 <- sin(v)

  M <- t(U) %*% solve(t(U) %*% U) %*% t(U)
  I <- diag(n)

  RO <- matrix(0, nrow = 2, ncol = 2)
  RO[1, 1] <- t(V1) %*% (I - M) %*% V1
  RO[2, 2] <- t(V2) %*% (I - M) %*% V2
  RO[1, 2] <- t(V1) %*% (I - M) %*% V2
  RO[2, 1] <- t(V2) %*% (I - M) %*% V1

  covrao <- (1 / (n - 2 * (order + 1))) * RO
  det_full <- det(covrao)

  Deterr <- c()

  for (j in 1:n) {
    uu <- u[-j]
    vv <- v[-j]

    order.matrix <- t(matrix(rep(c(1:order), n - 1), ncol = n - 1))
    cos.u <- cos(uu * order.matrix)
    sin.u <- sin(uu * order.matrix)

    ones <- matrix(1, n - 1, 1)
    U <- cbind(ones, cos.u, sin.u)

    V1 <- cos(vv)
    V2 <- sin(vv)

    M <- t(U) %*% solve(t(U) %*% U) %*% t(U)
    I <- diag(n - 1)

    RO <- matrix(0, nrow = 2, ncol = 2)
    RO[1, 1] <- t(V1) %*% (I - M) %*% V1
    RO[2, 2] <- t(V2) %*% (I - M) %*% V2
    RO[1, 2] <- t(V1) %*% (I - M) %*% V2
    RO[2, 1] <- t(V2) %*% (I - M) %*% V1

    covrao_j <- (1 / ((n - 1) - 2 * (order + 1))) * RO
    Deterr[j] <- det(covrao_j)
  }

  COVRATIO <- matrix(0, nrow = n)
  for (j in 1:n) {
    COVRATIO[j] <- Deterr[j] / det_full
  }

  p <- abs(COVRATIO - 1)
  maxP <- max(p)

  return(list(
    COVRATIO = COVRATIO,
    absolute_deviation = p,
    max_deviation = maxP
  ))
}

# ------------------------------------------------------------
# 2. EM algorithm for a two-component Von Mises mixture model
# ------------------------------------------------------------

run_em_von_mises <- function(iteration = 500,
                             n = 10000,
                             h = 10000,
                             k = 2,
                             meanobs1 = 0.5 * pi,
                             meanobs2 = pi,
                             kappaobs1 = 10.3,
                             kappaobs2 = 12.8,
                             pobs = 0.01) {

  Result <- matrix(1, h, 5)
  colnames(Result) <- c("mu1", "mu2", "kappa1", "kappa2", "p1")

  Result[1, ] <- c(0.035 * pi, 1.01 * pi, 10.3, 10, 0.1)

  for (e in 2:h) {
    p <- c(Result[(e - 1), 5], 1 - Result[(e - 1), 5])
    mu <- Result[(e - 1), 1:2]
    kappa <- Result[(e - 1), 3:4]

    z <- rbinom(n, 1, pobs)
    randata1 <- rvm(n, meanobs1, kappaobs1)
    randata2 <- rvm(n, meanobs2, kappaobs2)
    theta <- (z * randata1 + (1 - z) * randata2) %% (2 * pi)

    f <- matrix(0, n, k)
    for (j in 1:k) {
      f[, j] <- dvm(theta, mu[j], kappa[j])
    }

    tt <- matrix(0, n, k)
    for (i in 1:n) {
      for (j in 1:k) {
        tt[i, j] <- (p[j] * f[i, j]) / (p[1] * f[i, 1] + p[2] * f[i, 2])
      }
    }

    for (j in 1:k) {
      sorat <- sum(tt[, j] * sin(theta))
      makhrj <- sum(tt[, j] * cos(theta))
      mu[j] <- atan(sorat / makhrj)

      if (makhrj < 0) {
        mu[j] <- pi + atan(sorat / makhrj)
      } else {
        mu[j] <- ifelse(sorat > 0,
                        atan(sorat / makhrj),
                        2 * pi + atan(sorat / makhrj))
      }
    }

    Result[e, 1:2] <- c(mu[1], mu[2])

    d <- c()
    for (j in 1:k) {
      d[j] <- sum(tt[, j] * cos(theta - mu[j])) / sum(tt[, j])
      kappa[j] <- -1 * (20 * d[j] - 53 * d[j]^3) / (d[j]^2 - 1)
    }

    Result[e, 3:4] <- c(kappa[1], kappa[2])
    Result[e, 5] <- sum(tt[, 1]) / n

    if ((Result[e, ] - Result[(e - 1), ]) %*% (Result[e, ] - Result[(e - 1), ]) < 0.00001) {
      break
    }

    print(e)
  }

  M <- Result[e, ]
  print(M)

  zz <- max.col(tt)
  index <- 1:n
  index <- index[zz == 1]

  m <- c(meanobs1, meanobs2, kappaobs1, kappaobs2, pobs)

  MSE <- matrix(1, nrow = 1, ncol = 5)
  colnames(MSE) <- c("MSE.mu1", "MSE.mu2", "MSE.kappa1", "MSE.kappa2", "MSE.p1")

  for (j in 1:iteration) {
    for (i in 1:5) {
      MSE[1, i] <- (1 / j) * (sum((M[i] - m[i])^2))
    }
  }

  return(list(
    Result = Result[1:e, ],
    Final_Estimates = M,
    MSE = MSE,
    detected_index = index
  ))
}

# ------------------------------------------------------------
# Example usage
# ------------------------------------------------------------
# em_output <- run_em_von_mises(iteration = 500)
# print(em_output$Final_Estimates)
# print(em_output$MSE)

# u <- runif(100, 0, 2*pi)
# v <- runif(100, 0, 2*pi)
# covratio_output <- calculate_covratio(u, v)
# print(covratio_output$max_deviation)
