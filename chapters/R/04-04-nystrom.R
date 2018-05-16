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
get_time(mod.full); get_size(mod.full, "MB"); get_prederror(mod.full)
get_time(mod.nys); get_size(mod.nys); get_prederror(mod.nys)

## ---- nystrom.size.2 ----
# tab <- data.frame(
#   "Full" = c(
#     paste(dec_plac(get_time(mod.full)$time), get_time(mod.full)$unit),
#     capture.output(get_size(mod.full, units = "MB")),
#     dec_plac(get_prederror(mod.full), 3)
#     ),
#   "Nystrom" = c(
#     paste(dec_plac(get_time(mod.nys)$time), get_time(mod.nys)$unit),
#     capture.output(get_size(mod.nys)),
#     dec_plac(get_prederror(mod.nys), 3)
#   )
# )
# colnames(tab) <- c("Full model", "Nyström model")
# rownames(tab) <- c("Time taken", "Object size", "Training RMSE")
# kable(tab, align = "r", format = "latex", booktabs = TRUE) %>%
#   kable_styling(position = "center")

## ---- nystrom.plot ----
plot(mod.full) + ggtitle("Full model")
plot(mod.nys) + ggtitle("Nyström model")
