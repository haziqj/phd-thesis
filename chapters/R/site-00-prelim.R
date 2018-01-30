## ---- prelim ----
library(iprobit)
library(iprior)
library(ipriorBVS)
library(ggplot2)
library(directlabels)
gg_colour_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
# library(gganimate)
library(animation)
library(reshape2)
library(directlabels)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stan2coda <- function(fit) {
  mcmc.list(lapply(1:ncol(fit), function(x) mcmc(as.array(fit)[,x,])))
}
library(ggmcmc)
library(coda)

# Function to specify decimal places
decPlac <- function(x, k = 2) format(round(x, k), nsmall = k)

# Function to determine even numbers
isEven <- function(x) x %% 2 == 0

# SE kernel
fnH4 <- function(x, y = NULL, l = 1) {
  x <- scale(x, scale = FALSE)
  if (is.vector(x))
    x <- matrix(x, ncol = 1)
  n <- nrow(x)
  A <- matrix(0, n, n)
  index.mat <- upper.tri(A)
  index <- which(index.mat, arr.ind = TRUE)
  xcrossprod <- tcrossprod(x)
  if (is.null(y)) {
    tmp1 <- diag(xcrossprod)[index[, 1]]
    tmp2 <- diag(xcrossprod)[index[, 2]]
    tmp3 <- xcrossprod[index]
    A[index.mat] <- tmp1 + tmp2 - 2 * tmp3
    A <- A + t(A)
    tmp <- exp(-A / (2 * l ^ 2))
  } else {
    if (is.vector(y))
      y <- matrix(y, ncol = 1)
    else y <- as.matrix(y)
    y <- sweep(y, 2, attr(x, "scaled:center"), "-")
    m <- nrow(y)
    B <- matrix(0, m, n)
    indexy <- expand.grid(1:m, 1:n)
    ynorm <- apply(y, 1, function(z) sum(z ^ 2))
    xycrossprod <- tcrossprod(y, x)
    tmp1 <- ynorm[indexy[, 1]]
    tmp2 <- diag(xcrossprod)[indexy[, 2]]
    tmp3 <- as.numeric(xycrossprod)
    B[, ] <- tmp1 + tmp2 - 2 * tmp3
    tmp <- exp(-B / (2 * l ^ 2))
  }
  tmp
}

# Polynomial kernel
poly_kernel_deg <- 2
fnH5 <- function(x, y = NULL, c = 1, d = poly_kernel_deg, lambda = 1) {
  rownames(x) <- colnames(x) <- rownames(y) <- colnames(y) <- NULL
  x <- scale(x, scale = FALSE)  # centre the variables
  if (is.null(y)) {
    tmp <- (lambda * tcrossprod(x) + c) ^ d
  } else {
    if (is.vector(y)) y <- matrix(y, ncol = ncol(x))
    else y <- as.matrix(y)
    y <- sweep(y, 2, attr(x ,"scaled:center"), "-")
    tmp <- (lambda * tcrossprod(y, x) + c) ^ d
  }
  # class(tmp) <- "Canonical"
  tmp
}

dev_SEkern_iprior <- function(theta, y = y) {
  alpha <- mean(y)
  lambda <- exp(theta[1])
  psi <- exp(theta[2])
  n <- length(y)
  H <- fnH4(x, l = exp(theta[3]))
  tmp <- eigen(lambda * H)
  u <- psi * tmp$val ^ 2 + 1 / psi
  V <- tmp$vec
  res <- -(n / 2) * log(2 * pi) - (1 / 2) * sum(log(u)) -
    (1 / 2) * ((y - alpha) %*% V) %*% ((t(V) / u) %*% (y - alpha))
  as.numeric(-2 * res)
}

dev_FBMkern_iprior <- function(theta, y = y) {
  alpha <- mean(y)
  lambda <- exp(theta[1])
  psi <- exp(theta[2])
  # theta[3] <- 0
  n <- length(y)
  H <- fnH3(x, gamma = pnorm(theta[3]))
  tmp <- eigen(lambda * H)
  u <- psi * tmp$val ^ 2 + 1 / psi
  V <- tmp$vec
  res <- -(n / 2) * log(2 * pi) - (1 / 2) * sum(log(u)) -
    (1 / 2) * ((y - alpha) %*% V) %*% ((t(V) / u) %*% (y - alpha))
  as.numeric(-2 * res)
}

dev_Cankern_iprior <- function(theta, y = y) {
  alpha <- mean(y)
  lambda <- exp(theta[1])
  psi <- exp(theta[2])
  n <- length(y)
  H <- fnH2(x)
  tmp <- eigen(lambda * H)
  u <- psi * tmp$val ^ 2 + 1 / psi
  V <- tmp$vec
  res <- -(n / 2) * log(2 * pi) - (1 / 2) * sum(log(u)) -
    (1 / 2) * ((y - alpha) %*% V) %*% ((t(V) / u) %*% (y - alpha))
  as.numeric(-2 * res)
}

dev_Polykern_iprior <- function(theta, y = y) {
  alpha <- mean(y)
  lambda <- exp(theta[1])
  psi <- exp(theta[2])
  c <- exp(theta[3])
  n <- length(y)
  H <- fnH5(x, c = c, lambda = lambda)
  tmp <- eigen(H)
  u <- psi * tmp$val ^ 2 + 1 / psi
  V <- tmp$vec
  res <- -(n / 2) * log(2 * pi) - (1 / 2) * sum(log(u)) -
    (1 / 2) * ((y - alpha) %*% V) %*% ((t(V) / u) %*% (y - alpha))
  as.numeric(-2 * res)
}

move_fig_to_site <- function() {
  files <- list.files("image/")
  file.copy(file.path("image", files),
            file.path("/Users/haziqjamil/Desktop/Research/phd-poster/docs/assets/images", files),
            overwrite = TRUE)
}
