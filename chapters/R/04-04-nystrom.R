# Chapter 4
# Nystrom example
source("00-prelim.R")

## ---- nystrom.data ----
dat <- gen_smooth(n = 2000, xlim = c(-1, 5.5), seed = 1)
head(dat)
ggplot2::qplot(X, y, data = dat)
# NOTE: Fitting full model with n = 2000 takes roughly 15 minutes.

## ---- nystrom.mod.full ----
(mod.full <- iprior(y ~ X, dat, kernel = "fbm",
                    control = list(silent = TRUE)))

## ---- nystrom.mod ----
(mod.nys <- iprior(y ~ X, dat, kernel = "fbm", nystrom = 50,
                   control = list(silent = TRUE)))

## ---- nystrom.size ----
get_time(mod.full); get_size(mod.full, units = "MB")
get_time(mod.nys); get_size(mod.nys)

## ---- nystrom.plot ----
plot(mod.full)
plot(mod.nys)
