# Chapter 8
# Bayesian Variable Selection using I-priors simulation study
source("00-prelim.R")
library(ipriorBVS)
library(doSNOW)
library(foreach)
library(pushoverr)
library(glmnet)
library(BAS)
library(cowplot)
userID <- "uyq2g37vnityt1b3yvpyicv6o9h456"
appToken <- "avxnrig1qppsgsw9woghwwmxsobo4a"

## ---- sim.function ----
bvs_sim <- function(snr = 0.90, gamma = snr, n.sim = 25,
                    save.intermediate = TRUE) {
  res1 <- res2 <- matrix(NA, nrow = n.sim, ncol = 4)
  colnames(res1) <- colnames(res2) <- c("false.inc", "false.exc", "false", "brier")

  for (sim in seq_len(n.sim)) {
    cat("\n")
    cat("------ SIMULATION NUMBER", sim, paste0("(SNR: ", snr, ")"), "------")
    cat("\n")
    dat <- gen_bvs(n = 150, p = 100, snr = snr)

    # BVS I-prior (first stage) ------------------------------------------------
    runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE, modules = "lecuyer")
    mod <- ipriorBVS(y ~ ., dat, n.samp = 10000, n.burnin = 500,
                     n.adapt = 125, gamma.prob = gamma,
                     priors = "lambda ~ dunif(0, 100)", two.stage = FALSE)
    tmp <- unlist(predict(mod, dat$truth))[1:4]
    res1[sim, ] <- tmp
    cat("FIRST STAGE\n")
    print(res1[sim, ])
    cat("\n")

    # BVS I-prior (second stage) -----------------------------------------------
    runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE, modules = "lecuyer")
    mod <- ipriorBVS(y ~ ., dat, n.samp = 10000, n.burnin = 500,
                     n.adapt = 125, gamma.prob = gamma * get_mpm(mod),
                     priors = "lambda ~ dunif(0, 100)", two.stage = FALSE)
    tmp <- unlist(predict(mod, dat$truth))[1:4]
    res2[sim, ] <- tmp
    cat("SECOND STAGE\n")
    print(res2[sim, ])
    cat("\n")

    if (isTRUE(save.intermediate)) {
      save(res1, file = paste0("data/bvs-stage1-", snr, "-gam", gamma))
      save(res2, file = paste0("data/bvs-stage2-", snr, "-gam", gamma))
    }
  }

  push.message <- paste0(
    "BVS simulations (SNR: ", snr, ", gamma:", gamma, ") COMPLETED."
  )
  pushoverr::pushover(message = push.message, user = userID, app = appToken)

  list(stage1 = res1, stage2 = res2, snr = snr, gamma = gamma)
}

## ---- bvs.sims.gamma90 ----
res.90 <- bvs_sim(snr = 0.90, gamma = 0.9)
res.75 <- bvs_sim(snr = 0.75, gamma = 0.9)
res.50 <- bvs_sim(snr = 0.50, gamma = 0.9)
res.25 <- bvs_sim(snr = 0.25, gamma = 0.9)
res.10 <- bvs_sim(snr = 0.10, gamma = 0.9)
res.gam90 <- list("0.90" = res.90[1:2],
                  "0.75" = res.75[1:2],
                  "0.50" = res.50[1:2],
                  "0.25" = res.25[1:2],
                  "0.10" = res.10[1:2])
save(res.gam90, file = "data/bvs-res-gam90")

## ---- bvs.sims.gamma75 ----
res.90 <- bvs_sim(snr = 0.90, gamma = 0.75)
res.75 <- bvs_sim(snr = 0.75, gamma = 0.75)
res.50 <- bvs_sim(snr = 0.50, gamma = 0.75)
res.25 <- bvs_sim(snr = 0.25, gamma = 0.75)
res.10 <- bvs_sim(snr = 0.10, gamma = 0.75)
res.gam75 <- list("0.90" = res.90[1:2],
                  "0.75" = res.75[1:2],
                  "0.50" = res.50[1:2],
                  "0.25" = res.25[1:2],
                  "0.10" = res.10[1:2])
save(res.gam75, file = "data/bvs-res-gam75")

## ---- bvs.sims.gamma50 ----
res.90 <- bvs_sim(snr = 0.90, gamma = 0.5)
res.75 <- bvs_sim(snr = 0.75, gamma = 0.5)
res.50 <- bvs_sim(snr = 0.50, gamma = 0.5)
res.25 <- bvs_sim(snr = 0.25, gamma = 0.5)
res.10 <- bvs_sim(snr = 0.10, gamma = 0.5)
res.gam50 <- list("0.90" = res.90[1:2],
                  "0.75" = res.75[1:2],
                  "0.50" = res.50[1:2],
                  "0.25" = res.25[1:2],
                  "0.10" = res.10[1:2])
save(res.gam50, file = "data/bvs-res-gam50")

## ---- bvs.sims.gamma25 ----
res.90 <- bvs_sim(snr = 0.90, gamma = 0.25)
res.75 <- bvs_sim(snr = 0.75, gamma = 0.25)
res.50 <- bvs_sim(snr = 0.50, gamma = 0.25)
res.25 <- bvs_sim(snr = 0.25, gamma = 0.25)
res.10 <- bvs_sim(snr = 0.10, gamma = 0.25)
res.gam25 <- list("0.90" = res.90[1:2],
                  "0.75" = res.75[1:2],
                  "0.50" = res.50[1:2],
                  "0.25" = res.25[1:2],
                  "0.10" = res.10[1:2])
save(res.gam25, file = "data/bvs-res-gam25")

## ---- bvs.sims.gamma10 ----
res.90 <- bvs_sim(snr = 0.90, gamma = 0.1)
res.75 <- bvs_sim(snr = 0.75, gamma = 0.1)
res.50 <- bvs_sim(snr = 0.50, gamma = 0.1)
res.25 <- bvs_sim(snr = 0.25, gamma = 0.1)
res.10 <- bvs_sim(snr = 0.10, gamma = 0.1)
res.gam10 <- list("0.90" = res.90[1:2],
                  "0.75" = res.75[1:2],
                  "0.50" = res.50[1:2],
                  "0.25" = res.25[1:2],
                  "0.10" = res.10[1:2])
save(res.gam10, file = "data/bvs-res-gam10")

## ---- combine ----
# Only if need to combine two or more saved results
for (i in seq_along(res.gam10)) {
  res.gam25[[i]] <- mapply(rbind, res.gam25[[i]], res.gam25.2[[i]],
                           SIMPLIFY = FALSE)
}

## ---- bvs.sims.res ----
count_false_choices <- function(x) {
  x.0to2 <- mean(x <= 2)
  x.6 <- mean(x >= 6)
  x.3to5 <- 1 - x.0to2 - x.6
  ress <- c(x.0to2, x.3to5, x.6)
  names(ress) <- c("0-2", "3-5", ">5")
  ress
}

se_false_choices <- function(x) {
  n <- length(x)
  x.0to2 <- (sd(x <= 2) * (n - 1) / n) / sqrt(n)
  x.6 <- (sd(x >= 6) * (n - 1) / n) / sqrt(n)
  x.3to5 <- (x > 2) & (x < 6)
  x.3to5 <- (sd(x.3to5) * (n - 1) / n) / sqrt(n)
  ress <- c(x.0to2, x.3to5, x.6)
  names(ress) <- c("0-2", "3-5", ">5")
  ress
}

bvs_res <- function(x = res, type = mean, models = "all") {
  if (models == "all") {
    res <- lapply(x, function(y) t(sapply(y, function(z) apply(z, 2, type))))

    if (ncol(res[[1]]) > 3) {
      res <- lapply(res, function(y) {
        y <- y[, 7:9];
        colnames(y) <- c("0-2", "3-5", ">5");
        y
      })
    }
  } else {
    if (models == "iprior") tmp.res <- lapply(x, function(y) y[[1]])
    if (models == "flat_prior") tmp.res <- lapply(x, function(y) y[[2]])
    if (models == "gprior") tmp.res <- lapply(x, function(y) y[[3]])
    if (models == "lasso") tmp.res <- lapply(x, function(y) y[[4]])


    res <- t(sapply(tmp.res, function(x) apply(x, 2, type)))
    if (ncol(res) > 3) {
      res <- res[, 7:9]
      colnames(res) <- c("0-2", "3-5", ">5")
    }
  }
  res
}
t(bvs_res(models = "lasso", type = count_false_choices))
t(bvs_res(models = "lasso", type = se_false_choices))




## ---- bvs.results ----
plot.iprior.df <- rbind(
  data.frame(res[[1]][[1]], SNR = "90%"),
  data.frame(res[[2]][[1]], SNR = "75%"),
  data.frame(res[[3]][[1]], SNR = "50%"),
  data.frame(res[[4]][[1]], SNR = "25%"),
  data.frame(res[[5]][[1]], SNR = "10%")
)

plot.flat.df <- rbind(
  data.frame(res[[1]][[2]], SNR = "90%"),
  data.frame(res[[2]][[2]], SNR = "75%"),
  data.frame(res[[3]][[2]], SNR = "50%"),
  data.frame(res[[4]][[2]], SNR = "25%"),
  data.frame(res[[5]][[2]], SNR = "10%")
)

plot.gprior.df <- rbind(
  data.frame(res[[1]][[3]], SNR = "90%"),
  data.frame(res[[2]][[3]], SNR = "75%"),
  data.frame(res[[3]][[3]], SNR = "50%"),
  data.frame(res[[4]][[3]], SNR = "25%"),
  data.frame(res[[5]][[3]], SNR = "10%")
)

plot.lasso.df <- rbind(
  data.frame(res[[1]][[4]], SNR = "90%"),
  data.frame(res[[2]][[4]], SNR = "75%"),
  data.frame(res[[3]][[4]], SNR = "50%"),
  data.frame(res[[4]][[4]], SNR = "25%"),
  data.frame(res[[5]][[4]], SNR = "10%")
)

theme_set(theme_bw())

ggplot(plot.iprior.df, aes(false, ..density.., col = SNR)) +
  geom_freqpoly(binwidth = 1, size = 0.8) +
  scale_y_continuous(limits = c(0, 0.85)) +
  scale_x_continuous(limits = c(0, 90)) +
  ggtitle("I-prior") +
  labs(x = "Number of false choices", y = "Proportion") -> p1
# ggplot(plot.iprior.df, aes(brier, ..scaled.., fill = SNR, col = SNR)) +
#   geom_density(alpha = 0.3)

ggplot(plot.flat.df, aes(false, ..density.., col = SNR)) +
  geom_freqpoly(binwidth = 1, size = 0.8) +
  scale_y_continuous(limits = c(0, 0.85)) +
  scale_x_continuous(limits = c(0, 90)) +
  ggtitle("Independent prior") +
  labs(x = "Number of false choices", y = "Proportion") -> p2
ggplot(plot.flat.df, aes(brier, ..scaled.., fill = SNR, col = SNR)) +
  geom_density(alpha = 0.3)

ggplot(plot.gprior.df, aes(false, ..density.., col = SNR)) +
  geom_freqpoly(binwidth = 1, size = 0.8) +
  scale_y_continuous(limits = c(0, 0.85)) +
  scale_x_continuous(limits = c(0, 90)) +
  ggtitle("g-prior") +
  labs(x = "Number of false choices", y = "Proportion") -> p3
ggplot(plot.gprior.df, aes(brier, ..scaled.., fill = SNR, col = SNR)) +
  geom_density(alpha = 0.3)

ggplot(plot.lasso.df, aes(false, ..density.., col = SNR)) +
  geom_freqpoly(binwidth = 1, size = 0.8) +
  scale_y_continuous(limits = c(0, 0.85)) +
  scale_x_continuous(limits = c(0, 90)) +
  ggtitle("Lasso") +
  labs(x = "Number of false choices", y = "Proportion") -> p4

cowplot::plot_grid(p1, p2, p3, p4, ncol = 1)
