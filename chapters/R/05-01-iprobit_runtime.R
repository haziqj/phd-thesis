# Chapter 5
# Estimated run-time of iprobit models
source("00-prelim.R")

n <- c(30, seq(50, 500, by = 50))
m <- c(2, 3, 4, 5)
tab <- matrix(NA, nrow = length(n), ncol = length(m))
rownames(tab) <- paste("n =", n)
colnames(tab) <- paste("m =", m)

for (i in seq_along(n)) {
  for (j in seq_along(m)) {
    dat <- gen_mixture(n = n[i], m = m[j])
    mod <- iprobit(y ~ X1, dat, control = list(maxit = 5), silent = TRUE)
    tab[i, j] <- mod$time$time
    cat(i, j, " ")
    print(mod$time)
    cat("\n")
  }
}
