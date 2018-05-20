# Chapter 8
# Bayesian Variable Selection using I-priors simulation study
source("00-prelim.R")
chapter <- "06"

# ---- bvs.ozone ----

data(Ozone, package = "mlbench")
colnames(Ozone) <- c("Month", "DayMonth", "DayWeek", "Ozone", "PresVand",
                     "WindLAX", "HumLAX", "TempSand", "TempElMon", "ibhLAX",
                     "PresGrad", "ibtLAX", "VisLAX")
Ozone <- Ozone[complete.cases(Ozone), ]
Ozone <- as.data.frame(lapply(Ozone, as.numeric))
y <- Ozone$Ozone; y <- scale(y)
X <- Ozone[, -4]; X <- scale(X)
X.sq <- X ^ 2
colnames(X.sq) <- paste0(colnames(X), "2")
X.int.names <- X.int <- NULL
for (i in seq_len(ncol(X))) {
  for (j in seq_len(ncol(X))) if (j > i) {
    X.int <- cbind(X.int, X[, j] * X[, i])
    X.int.names <- c(X.int.names, paste0(colnames(X)[j], "_", colnames(X)[i]))
  }
}
colnames(X.int) <- X.int.names
X2 <- cbind(X, X.sq, X.int)
n <- nrow(Ozone)

# Simulations
n.sim <- 50
res1 <- res2 <- data.frame(matrix(NA, nrow = n.sim, ncol = 4))
runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
for (i in seq_len(n.sim)) {
  test.samp <- sample(1:n, size = 25)
  train.samp <- (1:n)[-test.samp]

  y.test <- y[test.samp]
  y.train <- y[train.samp]
  X.test <- X[test.samp, ]
  X.train <- X[train.samp, ]
  X2.test <- X2[test.samp, ]
  X2.train <- X2[train.samp, ]

  mod1 <- ipriorBVS(y.train, X.train, two.stage = TRUE)
  mod2 <- ipriorBVS(y.train, X2.train, two.stage = TRUE)

  hpm1 <- get_pmps(mod1)[1]
  hpm2 <- get_pmps(mod2)[1]
  pred1 <- get_predict(mod1, X.test, y.test)
  pred2 <- get_predict(mod2, X2.test, y.test)

  res1[i, ] <- c(names(hpm1), as.numeric(hpm1), get_R2(mod1), sqrt(pred1$mse))
  res2[i, ] <- c(names(hpm2), as.numeric(hpm2), get_R2(mod2), sqrt(pred2$mse))

  cat("SIMULATION ", i, "\n")
}
colnames(res1) <- colnames(res2) <- c("Model", "PostProb", "R2", "RMSE")
save(res1, file = "data/res1")
save(res2, file = "data/res2")

res1[, -1] <- as.data.frame(lapply(res1[, -1], as.numeric))
res2[, -1] <- as.data.frame(lapply(res2[, -1], as.numeric))

library(tibble)
res1 %>%
  group_by(Model) %>%
  summarise(prob = mean(`Post. Prob`), R2 = mean(R2), RMSE = mean(RMSE),
            prop = n() / 50) %>%
  arrange(desc(prob)) -> res1
res2 %>%
  group_by(Model) %>%
  summarise(prob = mean(`Post. Prob`), R2 = mean(R2), RMSE = mean(RMSE),
            prop = n() / 50) %>%
  arrange(desc(prob)) -> res2

# HPM for mod1
colnames(X)[as.logical(as.numeric(strsplit(res1[[1]][[1]], "")[[1]]))]

# HPM for mod2
colnames(X2)[as.logical(as.numeric(strsplit(res2[[1]][[1]], "")[[1]]))]


# ---- bvs.pollution ----
data("pollution", package = "iprior")
pollution.stand <- as.data.frame(scale(pollution))
mod <- ipriorBVS(Mortality ~ ., pollution.stand, stand.x = FALSE, stand.y = FALSE,
                 two.stage = TRUE)
coef(mod)


# ---- bvs.boston ----
data("BostonHousing2", package = "mlbench")
boston <- BostonHousing2[complete.cases(BostonHousing2), -c(1, 2, 5)]
boston$chas <- as.numeric(boston$chas)
mod <- ipriorBVS(cmedv ~ ., boston)

# ---- bvs.aerobic ----
data(aerobic, package = "ipriorBVS")
colnames(aerobic)[-1] <- paste0("X", 1:6)
(mod1 <- ipriorBVS(Oxygen ~ ., aerobic))
plot_coef2(mod1)
ggsave("figure/06-aerobic_coef.pdf", plot_coef2(mod1), "pdf",
       width = 7, height = 7)
move_fig_to_chapter()
(mod2 <- ipriorBVS(Oxygen ~ ., aerobic, two.stage = TRUE))
