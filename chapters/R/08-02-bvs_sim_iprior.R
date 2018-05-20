# Chapter 8
# Bayesian Variable Selection using I-priors simulation study (no analysis)
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
# load("data/bvs-res-gam10")
# res.gam10.2 <- res.gam10
# load("data/gam/bvs-res-gam10")
# for (i in seq_along(res.gam10)) {
#   res.gam10[[i]] <- mapply(rbind, res.gam10[[i]], res.gam10.2[[i]],
#                            SIMPLIFY = FALSE)
# }
