# Chapter 5, examples
# Prediction of cardiac arrhythmia
library(iprobit)
library(tidyverse)
library(reshape2)

## ---- data.cardiac ----
# 451 observations with 194 continuous covariates
load("data/Arrh194.RData")
experiment.name <- "Cardiac data"
X.orig <- ArrhDataNew$x
y <- ArrhDataNew$y
y <- y - 1  # convert to 0 and 1
y <- as.factor(y)
levels(y) <- c("Normal", "Arrhythmia")
N <- length(y)
n <- c(50, 100, 200)  # subsamples
X <- scale(X.orig)  # standardise data
X[, 36] <- rep(0, N); X[, 181] <- rep(0, N)  # all zeroes
summary(y)

## ---- data_sex ----
sex <- read.csv("data/arrhythmia.data.txt")[, 2]
sex <- factor(sex)
levels(sex) <- c("Male", "Female")

## ---- plot.cardiac ----
# set.seed(123)
# n.plot.size <- 3
# n.plot <- sort(sample(1:nrow(X), size = n.plot.size, FALSE))
# X.arr <- X[y == "Arrhythmia", ]; X.arr <- apply(X.arr, 1, mean)
# X.nor <- X[y == "Normal", ]; X.nor <- apply(X.nor, 1, mean)
n.plot <- c(1, 8)
plot.df <- data.frame(x = 1:ncol(X), t(X[n.plot, ]))
colnames(plot.df) <- c("x", "Arrhythmia", "Normal")
plot.df <- reshape2::melt(plot.df, id.vars = "x")
ggplot(plot.df, aes(x = x, y = value, col = variable)) +
  geom_line(alpha = 0.8) +
  directlabels::geom_dl(aes(label = variable),
                        method = list("last.bumpup", cex = 0.8,
                                      directlabels::dl.trans(x = x + 0.1))) +
  scale_x_continuous(limits = c(0, 210), breaks = NULL) +
  labs(y = "Standardised attribute values", x = NULL) +
  theme_bw() +
  theme(legend.position = "none")

## ---- mod.full.cardiac ----
set.seed(123)
(mod <- iprobit(y, X, kernel = "fbm", control = list(maxit = 100)))
# (mod.sex <- iprobit(y, X, sex, kernel = "fbm"))
# for (i in 1:20) {
#   samp <- sample(1:451, size = 200)
#   yy <- y[samp]
#   XX <- X[samp, ]
#   sexx <- sex[samp]
#   mod.sex <- iprobit(yy, XX, sexx, kernel = "fbm",
#                      control = list(maxit = 5 , silent = TRUE))
#   mod.pred <- predict(mod.sex, list(X[-samp, ], sex[-samp]), y[-samp])
#   print(mod.pred$error.rate)
# }

## ---- cardiac.mod.full.plot ----
# tmp <- cowplot::plot_grid(iplot_lb(mod, lab.pos = "down"), NULL,
#                           rel_widths = c(1, 0.06))
# cowplot::plot_grid(tmp, iplot_error(mod), nrow = 2)
iplot_lb(mod, lab.pos = "down")
iplot_error(mod)

## ---- simulations.cardiac ----
source("07-03-classification_simulation.R")
set.seed(456)
res.canonical <- my_iprobit_sim(nsim = 100, kernel = "canonical")
res.fbm       <- my_iprobit_sim(nsim = 100, kernel = "fbm")
res.se        <- my_iprobit_sim(nsim = 100, kernel = "se")

(tab <- tab_res("I-probit (linear)"  = res.canonical,
               "I-probit (fBm-0.5)" = res.fbm,
               "I-probit (SE-1.0)"  = res.se))

# Other results
knn.mean        <- c(40.64, 38.94, 35.76)
knn.se          <- c(0.33, 0.33, 0.36)
knn             <- mean_and_se(knn.mean, knn.se)
svm.mean        <- c(36.16, 35.64, 35.20)
svm.se          <- c(0.47, 0.39, 0.35)
svm             <- mean_and_se(svm.mean, svm.se)
svm.radial.mean <- c(48.39, 47.24, 46.85)
svm.radial.se   <- c(0.49, 0.46, 0.43)
svm.radial      <- mean_and_se(svm.radial.mean, svm.radial.se)
gpc.radial.mean <- c(37.28, 33.80, 29.31)
gpc.radial.se   <- c(0.42, 0.40, 0.35)
gpc.radial      <- mean_and_se(gpc.radial.mean, gpc.radial.se)
rf.mean         <- c(31.65, 26.72, 22.40)
rf.se           <- c(0.39, 0.29, 0.31)
rf              <- mean_and_se(rf.mean, rf.se)
nsc.mean        <- c(34.98, 33.00, 31.08) # nearest shrunken centroids
nsc.se          <- c(0.46, 0.440, 0.41)
nsc             <- mean_and_se(nsc.mean, nsc.se)
penlog.mean     <- c(34.92, 30.48, 26.12) # l1-penalised logistic regression
penlog.se       <- c(0.42, 0.34, 0.27)
penlog          <- mean_and_se(penlog.mean, penlog.se)

other.tab <- rbind(
  "k-nn"              = knn,
  "SVM"            = svm,
  # "Radial Support Vector Machine"   = svm.radial,
  "GP (radial)"        = gpc.radial,
  "Random forests"                    = rf,
  "NSC"        = nsc,
  "L-1 logistic"      = penlog
)
colnames(other.tab) <- colnames(tab$tab)

# Calculate ranks
tab.mean <- rbind(tab$tab.mean,
                  "k-nn"              = knn.mean,
                  "SVM"            = svm.mean,
                  # "Radial Support Vector Machine"   = svm.radial.mean,
                  "GP (radial)"        = gpc.radial.mean,
                  "Random forests"                    = rf.mean,
                  "NSC"        = nsc.mean,
                  "L-1 logistic"      = penlog.mean
)
tab.se <- rbind(tab$tab.se,
                "k-nn"              = knn.se,
                "SVM"            = svm.se,
                # "Radial Support Vector Machine"   = svm.radial.se,
                "GP (radial)"        = gpc.radial.se,
                "Random forests"                    = rf.se,
                "NSC"        = nsc.se,
                "L-1 logistic"      = penlog.se
)
tab.ranks <- tab_rank(tab.mean, tab.se)

# Tabulate results
tab.all <- cbind(rbind(tab$tab, other.tab), Rank = tab.ranks)

# Save results
save(tab.mean, file = "data/cardiac_mean_results")
save(tab.se, file = "data/cardiac_se_results")
save(tab.ranks, file = "data/cardiac_rank_results")
save(tab.all, file = "data/cardiac_results")

## ---- sim.res.cardiac ----
load("data/cardiac_mean_results")
load("data/cardiac_se_results")
load("data/cardiac_rank_results")
load("data/cardiac_results")
# knitr::kable(tab.all, align = "r")

## ---- plot_cardiac ----
plot.df <- cbind(tab.mean, id = rownames(tab.mean))
plot.df <- reshape2::melt(plot.df)
id2 <- plot.df$id
tmp <- levels(id2)
tmp[grepl("I-probit", tmp)] <- "I-probit"
levels(id2) <- tmp
plot.se <- reshape2::melt(tab.se)
plot.df <- cbind(plot.df, se = plot.se[, 2], id2 = id2)
plot.df$id <- factor(plot.df$id,
                     levels = names(sort(tab.ranks, decreasing = TRUE)))

ggplot(plot.df, aes(x = value / 100, y = id, col = id2, label = dec_plac(value))) +
  geom_point(size = 2) +
  geom_text(nudge_y = 0.25, size = 3) +
  geom_errorbarh(aes(xmin = (value - 1.96 * se) / 100,
                     xmax = (value + 1.96 * se) / 100,
                     height = 0)) +
  facet_grid(. ~ variable) +
  labs(x = "Misclassification rate", y = "Method") + guides(col = FALSE) +
  scale_x_continuous(labels = scales::percent) +
  coord_cartesian(xlim = c(0.21, 0.50)) +
  theme_bw()
