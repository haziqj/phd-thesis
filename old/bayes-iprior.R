## ---- prelim ----
library(mvtnorm)
library(iprior)
library(coda)
library(R2jags)
library(runjags)
library(superdiag)
library(mcmcplots)
library(ggmcmc)
library(lattice)
library(bayesplot)
library(reshape2)
library(rstan)
library(dplyr)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
stan2coda <- function(fit) {
  mcmc.list(lapply(1:ncol(fit), function(x) mcmc(as.array(fit)[,x,])))
}
decimal <- function(x, k = 2) format(round(x, k), nsmall = k)

## ---- sim.data ----
data(cats, package = "MASS")
str(cats)

## ---- iprior.mod ----
(mod.iprior <- iprior(Hwt ~ Bwt, data = cats,
                      control = list(silent = TRUE)))
logLik(mod.iprior)

## ---- iprior.jags ----
mod <- function() {
  for (i in 1:n) {
    w[i] ~ dnorm(0, psi)
    mu[i] <- theta * inprod(H[i, ], w[1:n])
    Y[i] ~ dnorm(mu[i], psi)
  }

  # Priors
  psi ~ dgamma(0.001, 0.001)
  theta ~ dunif(0, 10)
  lambda <- theta / psi
}
mod.data <- list(Y = cats$Hwt - mean(cats$Hwt), H = fnH2(cats$Bwt),
                 n = length(cats$Hwt))
mod.params <- c("lambda", "psi")

## ---- iprior.jags.fit ----
system.time(
  mod.fit <- jags.parallel(
    data = mod.data, inits = NULL, parameters.to.save = mod.params,
    model.file = mod, n.chains = 8, n.iter = 4000, n.thin = 20,
    DIC = FALSE
  )
)

## ---- iprior.jags.res ----
mod.fit
tmp <- as.mcmc(mod.fit, end = 150, thin = 20)
mod.fit.ggs <- ggs(tmp)
ggs_autocorrelation(mod.fit.ggs) + scale_x_continuous(breaks = c(0, 25, 50))
ggs_traceplot(mod.fit.ggs) + scale_x_continuous(labels = seq(2000, 4000, length = 5))

## ---- iprior.stan1 ----
# |lambda| ~ Cauchy(0, 5)
stan.iprior.mod1 <- "
  data {
    int<lower=0> n; // number of data
    vector[n] y; // responses
    matrix[n, n] H; // centred canonical kernel
  }
  parameters {
    real<lower=0> lambda;
    real<lower=0> psi;
  }
  transformed parameters {
    cov_matrix[n] Vy;
    Vy = psi * lambda * H * lambda * H + diag_matrix(rep_vector(1 / psi, n));
  }
  model {
    target += multi_normal_lpdf(y | rep_vector(0, n), Vy);
    target += cauchy_lpdf(lambda | 0, 5);
  }
"
m1 <- rstan::stan_model(model_code = stan.iprior.mod1)
m1@model_name <- "iprior.cauchy-flat"
stan.iprior.dat <- list(y = cats$Hwt - mean(cats$Hwt),
                        H = fnH2(cats$Bwt), n = length(cats$Hwt))

## ---- iprior.stan1.fit ----
system.time(
  fit.stan1 <- sampling(m1, data = stan.iprior.dat, pars = c("lambda", "psi"),
                        iter = 2000, chains = 8, thin = 2)
)

## ---- iprior.stan1.res ----
print(fit.stan1, pars = c("lambda", "psi"))
tmp <- stan2coda(fit.stan1)
mod.fit.ggs <- ggs(tmp[, c("lambda", "psi")])
ggs_autocorrelation(mod.fit.ggs) + scale_x_continuous(breaks = c(0, 25, 50))
ggs_traceplot(mod.fit.ggs) + scale_x_continuous(labels = seq(1000, 2000, length = 6))

## ---- iprior.stan2 ----
# |lambda| ~ Cauchy(1, 0.5)
stan.iprior.mod2 <- "
  data {
    int<lower=0> n; // number of data
    vector[n] y; // responses
    matrix[n, n] H; // centred canonical kernel
  }
  parameters {
    real<lower=0> lambda;
    real<lower=0> psi;
  }
  transformed parameters {
    cov_matrix[n] Vy;
    Vy = psi * lambda * H * lambda * H + diag_matrix(rep_vector(1 / psi, n));
  }
  model {
    target += multi_normal_lpdf(y | rep_vector(0, n), Vy);
    target += cauchy_lpdf(lambda | 1, 0.5);
  }
"
m2 <- rstan::stan_model(model_code = stan.iprior.mod2)
m2@model_name <- "iprior.cauchy-tight"
stan.iprior.dat <- list(y = cats$Hwt - mean(cats$Hwt),
                        H = fnH2(cats$Bwt), n = length(cats$Hwt))
system.time(
  fit.stan2 <- sampling(m2, data = stan.iprior.dat, pars = c("lambda", "psi"),
                        iter = 2000, chains = 8, thin = 2)
)
print(fit.stan2, pars = c("lambda", "psi"))
tmp <- stan2coda(fit.stan2)
mod.fit.ggs <- ggs(tmp[, c("lambda", "psi")])
ggs_autocorrelation(mod.fit.ggs) + scale_x_continuous(breaks = c(0, 25, 50))
ggs_traceplot(mod.fit.ggs) + scale_x_continuous(labels = seq(1000, 2000, length = 6))

## ---- iprior.stan3 ----
# |lambda| ~ N(0, 5^2)
stan.iprior.mod3 <- "
  data {
    int<lower=0> n; // number of data
    vector[n] y; // responses
    matrix[n, n] H; // centred canonical kernel
  }
  parameters {
    real<lower=0> lambda;
    real<lower=0> psi;
  }
  transformed parameters {
    cov_matrix[n] Vy;
    Vy = psi * lambda * H * lambda * H + diag_matrix(rep_vector(1 / psi, n));
  }
  model {
    target += multi_normal_lpdf(y | rep_vector(0, n), Vy);
    target += normal_lpdf(lambda | 0, 5);
  }
"
m3 <- rstan::stan_model(model_code = stan.iprior.mod3)
m3@model_name <- "iprior.normal-flat"
stan.iprior.dat <- list(y = cats$Hwt - mean(cats$Hwt),
                        H = fnH2(cats$Bwt), n = length(cats$Hwt))
system.time(
  fit.stan3 <- sampling(m3, data = stan.iprior.dat, pars = c("lambda", "psi"),
                        iter = 2000, chains = 8, thin = 2)
)
print(fit.stan3, pars = c("lambda", "psi"))
tmp <- stan2coda(fit.stan3)
mod.fit.ggs <- ggs(tmp[, c("lambda", "psi")])
ggs_autocorrelation(mod.fit.ggs) + scale_x_continuous(breaks = c(0, 25, 50))
ggs_traceplot(mod.fit.ggs) + scale_x_continuous(labels = seq(1000, 2000, length = 6))



## ---- iprior.stan4 ----
# |lambda| ~ N(1, 0.5^2)
stan.iprior.mod4 <- "
  data {
    int<lower=0> n; // number of data
    vector[n] y; // responses
    matrix[n, n] H; // centred canonical kernel
  }
  parameters {
    real<lower=0> lambda;
    real<lower=0> psi;
  }
  transformed parameters {
    cov_matrix[n] Vy;
    Vy = psi * lambda * H * lambda * H + diag_matrix(rep_vector(1 / psi, n));
  }
  model {
    target += multi_normal_lpdf(y | rep_vector(0, n), Vy);
    target += normal_lpdf(lambda | 1, 0.5);
  }
"
m4 <- rstan::stan_model(model_code = stan.iprior.mod4)
m4@model_name <- "iprior.normal-tight"
stan.iprior.dat <- list(y = cats$Hwt - mean(cats$Hwt),
                        H = fnH2(cats$Bwt), n = length(cats$Hwt))
system.time(
  fit.stan4 <- sampling(m4, data = stan.iprior.dat, pars = c("lambda", "psi"),
                        iter = 2000, chains = 8, thin = 2)
)
print(fit.stan4, pars = c("lambda", "psi"))
tmp <- stan2coda(fit.stan4)
mod.fit.ggs <- ggs(tmp[, c("lambda", "psi")])
ggs_autocorrelation(mod.fit.ggs) + scale_x_continuous(breaks = c(0, 25, 50))
ggs_traceplot(mod.fit.ggs) + scale_x_continuous(labels = seq(1000, 2000, length = 6))



## ---- iprior.stan5 ----
# lambda ~ Unif(0,10)
stan.iprior.mod5 <- "
  data {
    int<lower=0> n; // number of data
    vector[n] y; // responses
    matrix[n, n] H; // centred canonical kernel
  }
  parameters {
    real<lower=0,upper=10> lambda;
    real<lower=0> psi;
  }
  transformed parameters {
    cov_matrix[n] Vy;
    Vy = psi * lambda * H * lambda * H + diag_matrix(rep_vector(1 / psi, n));
  }
  model {
    target += multi_normal_lpdf(y | rep_vector(0, n), Vy);
  }
"
m5 <- rstan::stan_model(model_code = stan.iprior.mod5)
m5@model_name <- "iprior.unif-flat"
stan.iprior.dat <- list(y = cats$Hwt - mean(cats$Hwt),
                        H = fnH2(cats$Bwt), n = length(cats$Hwt))
system.time(
  fit.stan5 <- sampling(m5, data = stan.iprior.dat, pars = c("lambda", "psi"),
                        iter = 2000, chains = 8, thin = 2)
)
print(fit.stan5, pars = c("lambda", "psi"))
tmp <- stan2coda(fit.stan5)
mod.fit.ggs <- ggs(tmp[, c("lambda", "psi")])
ggs_autocorrelation(mod.fit.ggs) + scale_x_continuous(breaks = c(0, 25, 50))
ggs_traceplot(mod.fit.ggs) + scale_x_continuous(labels = seq(1000, 2000, length = 6))



## ---- iprior.stan6 ----
# lambda ~ Unif(0.5,1.5)
stan.iprior.mod6 <- "
  data {
    int<lower=0> n; // number of data
    vector[n] y; // responses
    matrix[n, n] H; // centred canonical kernel
  }
  parameters {
    real<lower=0.5,upper=1.5> lambda;
    real<lower=0> psi;
  }
  transformed parameters {
    cov_matrix[n] Vy;
    Vy = psi * lambda * H * lambda * H + diag_matrix(rep_vector(1 / psi, n));
  }
  model {
    target += multi_normal_lpdf(y | rep_vector(0, n), Vy);
  }
"
m6 <- rstan::stan_model(model_code = stan.iprior.mod6)
m6@model_name <- "iprior.unif-flat"
stan.iprior.dat <- list(y = cats$Hwt - mean(cats$Hwt),
                        H = fnH2(cats$Bwt), n = length(cats$Hwt))
system.time(
  fit.stan6 <- sampling(m6, data = stan.iprior.dat, pars = c("lambda", "psi"),
                        iter = 2000, chains = 8, thin = 2)
)
print(fit.stan6, pars = c("lambda", "psi"))
tmp <- stan2coda(fit.stan6)
mod.fit.ggs <- ggs(tmp[, c("lambda", "psi")])
ggs_autocorrelation(mod.fit.ggs) + scale_x_continuous(breaks = c(0, 25, 50))
ggs_traceplot(mod.fit.ggs) + scale_x_continuous(labels = seq(1000, 2000, length = 6))



## ---- iprior.stan7 ----
# LKJ prior
stan.iprior.mod7 <- "
  data {
    int<lower=0> n; // number of data
    vector[n] y; // responses
    matrix[n,n] H; // centred canonical kernel
  }
  parameters {
    real<lower=0> lambda;
    real<lower=0> psi;
  }
  transformed parameters {
    corr_matrix[n] R;
    cholesky_factor_corr[n] L;
    cov_matrix[n] Vy;
    vector<lower=0>[n] tau;
    Vy = psi * lambda * H * lambda * H + diag_matrix(rep_vector(1 / psi, n));
    tau = sqrt(diagonal(Vy));
    R = diag_matrix(1 ./ tau) * Vy * diag_matrix(1 ./ tau);
    L = cholesky_decompose(R);
  }
  model {
    target += multi_normal_lpdf(y | rep_vector(0, n), Vy);
    target += cauchy_lpdf(tau | 0, 5);
    target += lkj_corr_cholesky_lpdf(L | 1);
  }
"
m7 <- rstan::stan_model(model_code = stan.iprior.mod7)
m7@model_name <- "iprior.LKJ"
stan.iprior.dat <- list(y = cats$Hwt - mean(cats$Hwt),
                        H = fnH2(cats$Bwt), n = length(cats$Hwt))

## ---- iprior.stan7.fit ----
system.time(
  fit.stan7 <- sampling(m7, data = stan.iprior.dat, pars = c("lambda", "psi"),
                        iter = 2000, chains = 8, thin = 2)
)

## ---- iprior.stan7.res ----
print(fit.stan7, pars = c("lambda", "psi"))
tmp <- stan2coda(fit.stan7)
mod.fit.ggs <- ggs(tmp[, c("lambda", "psi")])
ggs_autocorrelation(mod.fit.ggs) + scale_x_continuous(breaks = c(0, 25, 50))
ggs_traceplot(mod.fit.ggs) + scale_x_continuous(labels = seq(1000, 2000, length = 6))

## ---- iprior.stan8 ----
# No priors, MAP = MLE
stan.iprior.mod8 <- "
  data {
    int<lower=0> n; // number of data
    vector[n] y; // responses
    matrix[n, n] H; // centred canonical kernel
  }
  parameters {
    real lambda;
    real<lower=0> psi;
  }
  transformed parameters {
    cov_matrix[n] Vy;
    Vy = psi * lambda * H * lambda * H + diag_matrix(rep_vector(1 / psi, n));
  }
  model {
    target += multi_normal_lpdf(y | rep_vector(0, n), Vy);
  }
"
m8 <- rstan::stan_model(model_code = stan.iprior.mod8)
m8@model_name <- "iprior.map"
stan.iprior.dat <- list(y = cats$Hwt - mean(cats$Hwt),
                        H = fnH2(cats$Bwt), n = length(cats$Hwt))

## ---- iprior.stan8.fit ----
mod.stan8 <- optimizing(m8, data = stan.iprior.dat, hessian = TRUE)
mod.stan8$par[names(mod.stan8$par) == c("lambda", "psi")]
mod.stan8$value

## ---- stan.compare ----
tab <- matrix(NA, ncol = 6, nrow = 9)

# EM algorithm
tmp <- summary(mod.iprior)
lambda <- mod.iprior$lambda
lambda.se <- tmp$coefficients[2, 2]
lambda.int <- lambda + lambda.se * c(-1.96, 1.96)
psi <- mod.iprior$psi
psi.se <- tmp$psi.and.se[2]
psi.int <- psi + psi.se * c(-1.96, 1.96)
tab[1, ] <- c(lambda, lambda.int, psi, psi.int)

# Models 1-7
for (i in 1:7) {
  eval(call("<-", "currentmod", as.name(paste0("fit.stan", i))))
  tmp <- summary(currentmod)$summary
  lambda <- tmp[1, 1]
  lambda.int <- tmp[1, c(4, 8)]
  psi <- tmp[2, 1]
  psi.int <- tmp[2, c(4, 8)]
  tab[i + 1, ] <- c(lambda, lambda.int, psi, psi.int)
}

# MAP (Model 8)
se <- diag(solve(-mod.stan8$hessian))
lambda <- mod.stan8$par[names(mod.stan8$par) == c("lambda")]
lambda.int <- lambda + se[1] * c(-1.96, 1.96)
psi <- mod.stan8$par[names(mod.stan8$par) == c("psi")]
psi.int <- psi + se[2] * c(-1.96, 1.96)
tab[9, ] <- c(lambda, lambda.int, psi, psi.int)

## ---- tecator ----
data(tecator, package = "caret")
endpoints <- as.data.frame(endpoints)
colnames(endpoints) <- c("water", "fat", "protein")
# n = 215, use first 160 for training
absorpTrain <- -t(diff(t(absorp)))	# this takes first differences using diff()
absorpTest <- absorpTrain[161:215,]
absorpTrain <- absorpTrain[1:160,]
# Responses
fatTrain <- endpoints$fat[1:160]
fatTest <- endpoints$fat[161:215]
# Other variables
waterTrain <- endpoints$water[1:160]
waterTest <- endpoints$water[161:215]

## ---- tecator.iprior1 ----
# Canonical kernel
mod1 <- kernL(y = fatTrain, absorpTrain)
mod1.fit <- ipriorOptim(mod1, control = list(silent = TRUE))
RMSE.Train1 <- mod1.fit$sigma
fatTestPredicted1 <- predict(mod1.fit, list(absorpTest))
RMSE.Test1 <- sqrt(mean((fatTestPredicted1 - fatTest) ^ 2))

## ---- tecator.iprior2 ----
# FBM kernel
mod2 <- kernL(y = fatTrain, absorpTrain, model = list(kernel = "FBM"))
mod2.fit <- ipriorOptim(mod2, control = list(silent = TRUE))
RMSE.Train2 <- mod2.fit$sigma
fatTestPredicted2 <- predict(mod2.fit, list(absorpTest))
RMSE.Test2 <- sqrt(mean((fatTestPredicted2 - fatTest) ^ 2))

## ---- tecator.iprior3 ----
# FBM kernel + extra covariate
mod3 <- kernL(y = fatTrain, absorpTrain, waterTrain,
              model = list(kernel = c("FBM", "Canonical")))
mod3.fit <- ipriorOptim(mod3, control = list(silent = TRUE))
RMSE.Train3 <- mod3.fit$sigma
fatTestPredicted3 <- predict(mod3.fit, list(absorpTest, waterTest))
RMSE.Test3 <- sqrt(mean((fatTestPredicted3 - fatTest) ^ 2))

## ---- tecator.stan1 ----
mod1.stan <- "
  data {
    int<lower=0> n; // number of data
    int<lower=0> m; // number of prediction points
    vector[n] y; // data y
    matrix[n,n] H; // the kernel matrix
    matrix[m,n] Hpred; // the kernel matrix for prediction
  }
  parameters {
    real alpha;
    real<lower=0> psi;
    real<lower=0> lambda;
  }
  transformed parameters {
    cov_matrix[n] Vy;
    Vy = psi * (lambda * H) * (lambda * H) + diag_matrix(rep_vector(1 / psi, n));
  }
  model {
    target += multi_normal_lpdf(y | rep_vector(alpha, n), Vy);
    target += cauchy_lpdf(lambda | 0, 10);
  }
  generated quantities {
    vector[n] w;
    matrix[n,n] Vw;
    vector[m] yhat;
    vector[m] mu;
    Vw = inverse(Vy);
    w = multi_normal_rng(psi * lambda * H * (Vw * (y - alpha)), Vw);
    mu = alpha + lambda * Hpred * w;
    for (i in 1:m)
      yhat[i] = normal_rng(mu[i], 1 / sqrt(psi));
  }
"
m <- stan_model(model_code = mod1.stan)
dat <- list(y = mod1$Y, H = mod1$Hl[[1]], Hpred = fnH2(mod1$x[[1]], absorpTest),
            n = mod1$n, m = length(fatTest))
## ---- tecator.stan1.fit ----
mod1.stan.fit <- sampling(m, data = dat, iter = 2000, chains = 4, thin = 2,
                          pars = c("alpha", "lambda", "psi", "yhat"))

## ---- tecator.stan1.res ----
# Diagnostics
# fit.ggs <- ggs(stan2coda(mod1.stan.fit)[, c("lambda", "psi")])
# ggs_running(fit.ggs)
# ggs_autocorrelation(fit.ggs)
# ggs_compare_partial(fit.ggs)
theta.hat <- summary(mod1.stan.fit)$summary[, 1]
wherey <- grep("yhat\\[", names(theta.hat))
fatTestPredicted1B <- theta.hat[wherey]
RMSE.Train1B <- 1 / sqrt(theta.hat["psi"])
(RMSE.Test1B <- sqrt(mean((fatTestPredicted1B - fatTest) ^ 2)))
# PPC
ppc_dens_overlay(fatTest, as.matrix(mod1.stan.fit)[, wherey]) +
  ggplot2::theme_gray() +
  labs(x = "value", y = "density") +
  theme(legend.text.align = 0)

## ---- tecator.compare ----
tab <- c(coef(mod1.fit)[-1], RMSE.Test1)
tab <- rbind(tab, c(theta.hat[c(2, 3)], RMSE.Test1B))
colnames(tab) <- c("lambda", "psi", "RMSE Test")
rownames(tab) <- c("MLE", "Bayes")
tab

## ---- tecator.stan2 ----
# FBM RKHS
mod.stan <- "
  data {
    int<lower=0> n; // number of data
    int<lower=0> m; // number of prediction points
    vector[n] y; // data y
    matrix[n,n] H; // the kernel matrix
    matrix[m,n] Hpred; // the kernel matrix for prediction
  }
  parameters {
    real alpha;
    real<lower=0> psi;
    real<lower=0> lambda;
  }
  transformed parameters {
    cov_matrix[n] Vyinv;
    vector[n] u;
    matrix[n,n] V;
    u = eigenvalues_sym(lambda * H);
    V = eigenvectors_sym(lambda * H);
    Vyinv = diag_pre_multiply(1 ./ (psi * square(u) + 1 / psi), V) * V';
  }
  model {
    target += cauchy_lpdf(lambda | 0, 5);
    // target += cauchy_lpdf(psi | 0, 5);
    target += multi_normal_prec_lpdf(y | rep_vector(alpha, n), Vyinv);
  }
  generated quantities {
  //  vector[n] w;
  //  cov_matrix[n] Vw;
  //  vector[m] yhat;
  //  vector[m] mu;
  //  Vw = inverse(psi * (lambda * H) * (lambda * H) + diag_matrix(rep_vector(1 / psi, n)));
  //  w = multi_normal_rng(psi * lambda * H * (Vw * (y - alpha)), Vw);
  //  mu = alpha + lambda * Hpred * w;
  //  for (i in 1:m)
  //  yhat[i] = normal_rng(mu[i], 1 / sqrt(psi));
  }
"
m <- stan_model(model_code = mod.stan)
dat <- list(y = mod2$Y, H = mod2$Hl[[1]], Hpred = fnH3(mod2$x[[1]], absorpTest),
            n = mod2$n, m = length(fatTest))
mod2.stan.fit <- sampling(m, data = dat, iter = 200, chains = 1, thin = 1,
                          pars = c("alpha", "lambda", "psi"))
# Diagnostics
fit.ggs <- ggs(stan2coda(mod2.stan.fit)[, c("lambda", "psi")])
ggs_running(fit.ggs)
ggs_autocorrelation(fit.ggs)
ggs_compare_partial(fit.ggs)
theta.hat <- summary(mod2.stan.fit)$summary[, 1]
wherey <- grep("yhat\\[", names(theta.hat))
fatTestPredicted2B <- theta.hat[wherey]
RMSE.Train2B <- 1 / sqrt(theta.hat["psi"])
(RMSE.Test2B <- sqrt(mean((fatTestPredicted2B - fatTest) ^ 2)))
# PPC
ppc_dens_overlay(fatTest, as.matrix(mod2.stan.fit)[, wherey]) +
  ggplot2::theme_gray() +
  labs(x = "value", y = "density") +
  theme(legend.text.align = 0)


## ---- tecator.stan3 ----
# FBM RKHS
mod.stan <- "
  data {
    int<lower=0> n; // number of data
    vector[n] y; // data y
    matrix[n,n] H1; // the kernel matrix
    matrix[n,n] H2; // the kernel matrix for prediction
  }
  parameters {
    real alpha;
    real<lower=0> psi;
    vector<lower=0>[2] lambda;
  }
  transformed parameters {
    cov_matrix[n] Vy;
    cholesky_factor_cov[n] L;
    matrix[n,n] Hlam;
    Hlam = lambda[1] * H1 + lambda[2] * H2;
    Vy = psi * Hlam * Hlam + diag_matrix(rep_vector(1 / psi, n));
    L = cholesky_decompose(Vy);
  }
  model {
    target += cauchy_lpdf(lambda | 0, 5);
    target += cauchy_lpdf(psi | 0, 5);
    target += multi_normal_cholesky_lpdf(y | rep_vector(alpha, n), L);
  }
"
m <- stan_model(model_code = mod.stan)
dat <- list(y = mod2$Y, H1 = mod3$Hl[[1]], H2 = mod3$Hl[[2]], n = mod2$n)
mod3.stan.fit <- sampling(m, data = dat, iter = 200, chains = 4, thin = 1,
                          pars = c("alpha", "lambda", "psi"))
# Diagnostics
fit.ggs <- ggs(stan2coda(mod3.stan.fit)[, c("lambda[1]", "lambda[2]", "psi")])
ggs_running(fit.ggs)
ggs_autocorrelation(fit.ggs)
ggs_compare_partial(fit.ggs)
theta.hat <- summary(mod2.stan.fit)$summary[, 1]
wherey <- grep("yhat\\[", names(theta.hat))
fatTestPredicted2B <- theta.hat[wherey]
RMSE.Train2B <- 1 / sqrt(theta.hat["psi"])
(RMSE.Test2B <- sqrt(mean((fatTestPredicted2B - fatTest) ^ 2)))
# PPC
ppc_dens_overlay(fatTest, as.matrix(mod2.stan.fit)[, wherey]) +
  ggplot2::theme_gray() +
  labs(x = "value", y = "density") +
  theme(legend.text.align = 0)


## ---- stabilise.iprior ----
model_code <- "
functions {
  real iprior_lpdf(vector y, int n, real alpha, real lambda, real psi, matrix H) {
    vector[n] u;
    matrix[n,n] V;
    vector[n] us;
    real logdet;
    vector[n] a;
    u = eigenvalues_sym(lambda * H);
    V = eigenvectors_sym(lambda * H);
    for (i in 1:n)
      us[i] = psi * u[i] * u[i] + 1 / psi;
    logdet = sum(log(us));
    a = V * diag_matrix(1 ./ us) * V' * (y - alpha);
    return -0.5 * n * log(2 * pi()) - 0.5 * logdet - 0.5 * ((y - alpha)' * a);
  }
  }
model {}
"
expose_stan_functions(stanc(model_code = model_code))

## ---- explore.condition ----
ipriorCon <- function(lambda, psi) {
  H <- mod2$Hl[[1]]
  u <- eigen(lambda * H, only.values = TRUE)[[1]]
  us <- psi * u ^ 2 + 1 / psi
  list(max = max(us), min = min(us), con = max(us) / min(us))
}



## ---- concrete ----
# Model: SLUMP.cm ~ X
concrete <- read.table("/Users/haziqjamil/Dropbox/LSE/Year 3 (PhD)/bayes-iprior/slump_test.data.txt", sep = ",", header = TRUE)
concrete <- concrete[, -1]
str(concrete)
# Training set observations 1 to 78
Xtrain <- scale(concrete[1:78, 1:7])
ytrain <- concrete[1:78, 8]
# Test set observations 79 to 103
Xtest <- scale(concrete[79:103, 1:7])
ytest <- concrete[79:103, 8]

## ---- concrete.iprior1 ----
# Canonical kernel, single lambda
mod1 <- kernL(y = ytrain, X = Xtrain)
mod1.fit <- iprior(mod1)
RMSE.train1 <- mod1.fit$sigma
slump.test1 <- predict(mod1.fit, list(Xtest))
RMSE.test1 <- sqrt(mean((slump.test1 - ytest) ^ 2))

## ---- concrete.iprior2 ----
# FBM kernel, single lambda
mod2 <- kernL(y = ytrain, X = Xtrain, model = list(kernel = "FBM"))
mod2.fit <- iprior(mod2)
RMSE.train2 <- mod2.fit$sigma
slump.test2 <- predict(mod2.fit, list(Xtest))
RMSE.test2 <- sqrt(mean((slump.test2 - ytest) ^ 2))

## ---- concrete.iprior3 ----
# FBM kernel, single lambda
mod3.fit <- fbmOptim(mod2)
RMSE.train3 <- mod3.fit$sigma
slump.test3 <- predict(mod3.fit, list(Xtest))
RMSE.test3 <- sqrt(mean((slump.test3 - ytest) ^ 2))

## ---- concrete.stan1 ----
mod1.stan <- "
  data {
    int<lower=0> n; // number of data
    int<lower=0> m; // number of prediction points
    vector[n] y; // data y
    matrix[n,n] H; // the kernel matrix
    matrix[m,n] Hpred; // the kernel matrix for prediction
  }
  parameters {
    real alpha;
    real<lower=0> psi;
    real<lower=0> lambda;
  }
  transformed parameters {
    cov_matrix[n] Vy;
    cholesky_factor_cov[n] L;
    Vy = psi * (lambda * H) * (lambda * H) + diag_matrix(rep_vector(1 / psi, n));
    L = cholesky_decompose(Vy);
  }
  model {
    target += multi_normal_cholesky_lpdf(y | rep_vector(alpha, n), L);
    target += cauchy_lpdf(lambda | 0, 10);
  }
  generated quantities {
    vector[n] w;
    matrix[n,n] Vw;
    vector[m] yhat;
    vector[m] mu;
    Vw = inverse(Vy);
    w = multi_normal_rng(psi * lambda * H * (Vw * (y - alpha)), Vw);
    mu = alpha + lambda * Hpred * w;
    for (i in 1:m)
    yhat[i] = normal_rng(mu[i], 1 / sqrt(psi));
  }
"
m <- stan_model(model_code = mod1.stan)
dat <- list(y = ytrain, H = fnH2(Xtrain), Hpred = fnH2(Xtrain, Xtest),
            n = length(ytrain), m = length(ytest))
mod1.stan.fit <- sampling(m, data = dat, iter = 2000, chains = 4, thin = 2,
                          pars = c("alpha", "lambda", "psi", "yhat"))
# Diagnostics
# fit.ggs <- ggs(stan2coda(mod1.stan.fit)[, c("lambda", "psi")])
# ggs_running(fit.ggs)
# ggs_autocorrelation(fit.ggs)
# ggs_compare_partial(fit.ggs)
theta.hat <- summary(mod1.stan.fit)$summary[, 1]
wherey <- grep("yhat\\[", names(theta.hat))
slump.test1B <- theta.hat[wherey]
RMSE.Train1B <- 1 / sqrt(theta.hat["psi"])
(RMSE.Test1B <- sqrt(mean((slump.test1B - ytest) ^ 2)))
# PPC
ppc_dens_overlay(ytest, as.matrix(mod1.stan.fit)[, wherey]) +
  ggplot2::theme_gray() +
  labs(x = "value", y = "density") +
  theme(legend.text.align = 0)


## ---- concrete.stan2 ----
mod2.stan <- "
  data {
    int<lower=0> n; // number of data
    int<lower=0> m; // number of prediction points
    vector[n] y; // data y
    matrix[n,n] H; // the kernel matrix
    matrix[m,n] Hpred; // the kernel matrix for prediction
  }
  parameters {
    real alpha;
    real<lower=0> psi;
    real<lower=0> lambda;
  }
  transformed parameters {
    cov_matrix[n] Vy;
    cholesky_factor_cov[n] L;
    Vy = psi * (lambda * H) * (lambda * H) + diag_matrix(rep_vector(1 / psi, n));
    L = cholesky_decompose(Vy);
  }
  model {
    target += multi_normal_cholesky_lpdf(y | rep_vector(alpha, n), L);
    target += cauchy_lpdf(lambda | 0, 10);
  }
  generated quantities {
    vector[n] w;
    matrix[n,n] Vw;
    vector[m] yhat;
    vector[m] mu;
    Vw = inverse(Vy);
    w = multi_normal_rng(psi * lambda * H * (Vw * (y - alpha)), Vw);
    mu = alpha + lambda * Hpred * w;
    for (i in 1:m)
    yhat[i] = normal_rng(mu[i], 1 / sqrt(psi));
  }
"
m <- stan_model(model_code = mod2.stan)
dat <- list(y = ytrain, H = fnH3(Xtrain), Hpred = fnH3(Xtrain, Xtest),
            n = length(ytrain), m = length(ytest))
mod2.stan.fit <- sampling(m, data = dat, iter = 2000, chains = 4, thin = 2,
                          pars = c("alpha", "lambda", "psi", "yhat"))
# Diagnostics
fit.ggs <- ggs(stan2coda(mod2.stan.fit)[, c("lambda", "psi")])
ggs_running(fit.ggs)
ggs_autocorrelation(fit.ggs)
ggs_compare_partial(fit.ggs)
theta.hat <- summary(mod2.stan.fit)$summary[, 1]
wherey <- grep("yhat\\[", names(theta.hat))
slump.test2B <- theta.hat[wherey]
RMSE.Train2B <- 1 / sqrt(theta.hat["psi"])
(RMSE.Test2B <- sqrt(mean((slump.test2B - ytest) ^ 2)))
# PPC
ppc_dens_overlay(ytest, as.matrix(mod2.stan.fit)[, wherey]) +
  ggplot2::theme_gray() +
  labs(x = "value", y = "density") +
  theme(legend.text.align = 0)

## ---- concrete.stan3 ----
mod3.stan <- "
  functions {
    matrix fnH3stan(vector x, real gamma, int n) {
      matrix[n,n] A;
      vector[n] r;
      real s;
      matrix[n,n] H;
      for (i in 1:n)
      for (j in 1:n)
      A[i,j] = pow(fabs(x[i] - x[j]), 2 * gamma);
      for(i in 1:n)
      r[i] = sum(A[,i]);
      s = sum(r);
      H = A - r * rep_vector(1, n)' / n - rep_vector(1, n) * r' / n;
      return -0.5 * (H + (s / square(n)));
    }
  }
  data {
    int<lower=0> n; // number of data
    int<lower=0> p; // number of covariates
    vector[n] y; // data y
    matrix[n,p] X; // the kernel matrix
  }
  parameters {
    real alpha;
    real<lower=0> psi;
    real<lower=0> lambda;
    real<lower=0,upper=1> gamma;
  }
  transformed parameters {
    matrix[n,n] H;
    cov_matrix[n] Vy;
    cholesky_factor_cov[n] L;
    H = fnH3stan(col(X, 1), gamma, n);
    for (j in 2:p)
      H = H + fnH3stan(col(X, j), gamma, n);
    Vy = psi * (lambda * H) * (lambda * H) + diag_matrix(rep_vector(1 / psi, n));
    L = cholesky_decompose(Vy);
  }
  model {
    target += multi_normal_cholesky_lpdf(y | rep_vector(alpha, n), L);
    target += cauchy_lpdf(lambda | 0, 10);
  }
  generated quantities {
    vector[n] u;
    vector[n] w;
    vector[n] mw;
    mw = psi * lambda * H * mdivide_left_spd(Vy, y - alpha);
    for (i in 1:n)
      u[i] = normal_rng(0, 1);
    w = mw +  inverse(L)' * u;
  }
"
m <- stan_model(model_code = mod3.stan)
dat <- list(y = ytrain, X = Xtrain, p = ncol(Xtrain), n = length(ytrain))
mod3.stan.fit <- sampling(m, data = dat, iter = 2000, chains = 4, thin = 2,
                          pars = c("alpha", "lambda", "psi", "gamma", "w"))
# Diagnostics
# fit.ggs <- ggs(stan2coda(mod3.stan.fit)[, c("lambda", "psi", "gamma")])
# ggs_running(fit.ggs)
# ggs_autocorrelation(fit.ggs)
# ggs_compare_partial(fit.ggs)
theta.hat <- summary(mod3.stan.fit)$summary[, 1]
draws <- extract(mod3.stan.fit)
ytestB <- matrix(NA, nrow = 2000, ncol = length(ytest))
for (i in 1:2000) {
  H <- fnH3(Xtrain, Xtest, gamma = draws$gamma[i])
  mu <- as.numeric(draws$alpha[i] + draws$lambda[i] * H %*% draws$w[i, ])
  ytestB[i, ] <- mu + rnorm(length(ytest), 0, 1 / sqrt(draws$psi[i]))
}
RMSE.Train3B <- 1 / sqrt(theta.hat["psi"])
(RMSE.Test3B <- sqrt(mean((slump.test3B - ytest) ^ 2)))
# PPC
ppc_dens_overlay(ytest, ytestB) +
  ggplot2::theme_gray() +
  labs(x = "value", y = "density") +
  theme(legend.text.align = 0)

## ---- fbmsim ----
data(datfbm)
str(datfbm)
fx <- function(x) 65 * dnorm(x, mean = 2) + 35 * dnorm(x, mean = 7, sd = 1.5)
y.true <- as.numeric(fx(datfbm$x))
datfbm <- cbind(datfbm, y.true)
ggplot(datfbm, aes(x = x, y = y)) +
  geom_point() +
  geom_line(aes(y = y.true), col = rainbow(1), linetype = 2)

## ---- datfbm.iprior1 ----
mod1.fit <- iprior(y ~ x, datfbm, model = list(kernel = "FBM"))


## ---- datfbm.iprior2 ----
mod2.fit <- fbmOptim(kernL(y ~ x, datfbm, model = list(kernel = "FBM")))

## ---- datfbm.stan1 ----
# FBM RKHS
mod.stan <- "
  data {
    int<lower=0> n; // number of data
    vector[n] y; // data y
    matrix[n,n] H; // the kernel matrix
  }
  parameters {
    real alpha;
    real<lower=0> psi;
    real<lower=0> lambda;
  }
  transformed parameters {
    cov_matrix[n] Vy;
    Vy = psi * lambda * H * lambda * H + diag_matrix(rep_vector(1 / psi, n));
  }
  model {
    target += cauchy_lpdf(lambda | 0, 5);
    target += cauchy_lpdf(psi | 0, 5);
    target += multi_normal_lpdf(y | rep_vector(alpha, n), Vy);
  }
  generated quantities {
    vector[n] w;
    vector[n] mw;
    cov_matrix[n] Vw;
    vector[n] yhat;
    vector[n] mu;
    Vw = inverse_spd(Vy);
    mw = psi * lambda * H * mdivide_left_spd(Vy, (y - alpha));
    w = multi_normal_rng(mw, Vw);
    mu = alpha + lambda * H * w;
    for (i in 1:n)
      yhat[i] = normal_rng(mu[i], 1 / sqrt(psi));
  }
"
m <- stan_model(model_code = mod.stan)
m@model_name <- "iprior.fbm"
dat <- list(y = datfbm$y, H = fnH3(datfbm$x), n = length(datfbm$y))
mod.stan.fit <- sampling(m, data = dat, iter = 2000, chains = 4, thin = 2,
                         pars = c("alpha", "lambda", "psi", "yhat"))

## ---- datfbm.stan1.res ----
# print(mod.stan.fit, pars = c("lambda", "psi"))
# Diagnostics
# fit.ggs <- ggs(stan2coda(mod.stan.fit)[, c("lambda", "psi")])
# ggs_running(fit.ggs)
# ggs_autocorrelation(fit.ggs)
# ggs_compare_partial(fit.ggs)
# ggs_density(fit.ggs)
theta.hat <- summary(mod.stan.fit)$summary[, 1]
wherey <- grep("yhat\\[", names(theta.hat))
# PPC
ppc_dens_overlay(datfbm$y, as.matrix(mod.stan.fit)[, wherey]) +
  ggplot2::theme_gray() +
  labs(x = "value", y = "density") +
  theme(legend.text.align = 0)
# PPC2
# predictions <- as.data.frame(mod.stan.fit) %>%
#   select(contains("yhat"))
# predictions <- t(predictions) %>%
#   melt()
# predictions <- cbind(predictions, x = datfbm$x)
datfbm.pred <- cbind(datfbm, summary(mod.stan.fit)$summary[wherey, c(1, 4, 8)])
names(datfbm.pred) <- c("y", "x", "y.true", "y.hat", "lwr", "upr")
ggplot(datfbm.pred, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_point() +
  geom_line(aes(y = y.hat, col = "1")) +
  geom_line(aes(y = y.true, col = "2"),  linetype = 2) +
  scale_colour_manual(name = "", values = c("1" = "black", "2" = rainbow(1)),
                      label = c("Fitted regression curve", "True regression curve"))

## ---- datfbm.stan2 ----
# FBM RKHS estimate gamma
mod.stan <- "
  functions {
    matrix fnH3(vector x, real gamma, int n) {
      matrix[n,n] A;
      vector[n] r;
      real s;
      matrix[n,n] H;
      for (i in 1:n)
        for (j in 1:n)
          A[i,j] = pow(fabs(x[i] - x[j]), 2 * gamma);
      for(i in 1:n)
        r[i] = sum(A[,i]);
      s = sum(r);
      H = A - r * rep_vector(1, n)' / n - rep_vector(1, n) * r' / n;
      return -0.5 * (H + (s / square(n)));
    }
  }
  data {
    int<lower=0> n; // number of data
    vector[n] y; // data y
    vector[n] X; // data x
  }
  parameters {
    real alpha;
    real<lower=0> psi;
    real<lower=0> lambda;
    real<lower=0,upper=1> gamma;
  }
  transformed parameters {
    cov_matrix[n] Vy;
    matrix[n,n] H;
    H = fnH3(X, gamma, n);
    Vy = psi * lambda * H * lambda * H + diag_matrix(rep_vector(1 / psi, n));
  }
  model {
    target += cauchy_lpdf(lambda | 0, 5);
    target += cauchy_lpdf(psi | 0, 5);
    target += multi_normal_lpdf(y | rep_vector(alpha, n), Vy);
  }
  generated quantities {
    vector[n] w;
    vector[n] mw;
    cov_matrix[n] Vw;
    vector[n] yhat;
    vector[n] mu;
    Vw = inverse_spd(Vy);
    mw = psi * lambda * H * mdivide_left_spd(Vy, (y - alpha));
    w = multi_normal_rng(mw, Vw);
    mu = alpha + lambda * H * w;
    for (i in 1:n)
      yhat[i] = normal_rng(mu[i], 1 / sqrt(psi));
  }
"
m <- stan_model(model_code = mod.stan)
dat <- list(y = datfbm$y, X = datfbm$x, n = length(datfbm$y))
mod.stan.fit <- sampling(m, data = dat, iter = 2000, chains = 4, thin = 2,
                         pars = c("alpha", "lambda", "psi", "gamma", "yhat"))
m@model_name <- "iprior.fbm.gamma"
print(mod.stan.fit, pars = c("lambda", "psi", "gamma"))
# Diagnostics
fit.ggs <- ggs(stan2coda(mod.stan.fit)[, c("lambda", "psi", "gamma")])
ggs_running(fit.ggs)
ggs_autocorrelation(fit.ggs)
ggs_compare_partial(fit.ggs)
ggs_density(fit.ggs)
theta.hat2 <- summary(mod.stan.fit)$summary[, 1]
wherey2 <- grep("yhat\\[", names(theta.hat2))
# PPC
ppc_dens_overlay(datfbm$y, as.matrix(mod.stan.fit)[, wherey2]) +
  ggplot2::theme_gray() +
  labs(x = "value", y = "density") +
  theme(legend.text.align = 0)
# PPC2
# predictions <- as.data.frame(mod.stan.fit) %>%
#   select(contains("yhat"))
# predictions <- t(predictions) %>%
#   melt()
# predictions <- cbind(predictions, x = datfbm$x)
datfbm.pred <- cbind(datfbm, summary(mod.stan.fit)$summary[wherey2, c(1, 4, 8)])
names(datfbm.pred) <- c("y", "x", "y.true", "y.hat", "lwr", "upr")
ggplot(datfbm.pred, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.3) +
  geom_point() +
  geom_line(aes(y = y.hat)) +
  geom_line(aes(y = y.true), col = rainbow(1), linetype = 2)

## ---- compare.fbm ----
mse.ml <- mean((fitted(mod1.fit) - y.true) ^ 2)
mse.mlg <- mean((fitted(mod2.fit) - y.true) ^ 2)
mse.bayes <- mean((theta.hat[wherey] - y.true) ^ 2)
mse.bayesg <- mean((theta.hat2[wherey2] - y.true) ^ 2)
tab <- c(coef(mod1.fit)[-1], 0.5, mse.ml)
tab <- rbind(tab, c(coef(mod2.fit)[-1], mod2.fit$ipriorKernel$model$Hurst, mse.mlg))
tab <- rbind(tab, c(theta.hat[-c(1, wherey, 104)], 0.5, mse.bayes))
tab <- rbind(tab, c(theta.hat2[-c(1, wherey2, 105)], mse.bayesg))
rownames(tab) <- c("MLE", "", "Bayes", "")
colnames(tab) <- c("lambda", "psi", "Hurst", "MSE")
tab














# prob <- c(0.05, 0.25, 0.55, 0.075, 0.075)
# grade <- c("A", "B", "C", "D", "E")
# value <- c(7.5, 15, 70, 95)
# x <- sample(grade, 10000, prob = prob, replace = TRUE)
# scoreFn <- function(x) {
#   if (x == "A") y <- rnorm(1, (95 + 100) / 2, 5)
#   else if (x == "B") y <- rnorm(1, (70 + 95) / 2, 10)
#   else if (x == "C") y <- rnorm(1, (15 + 70) / 2, 15)
#   else if (x == "D") y <- rnorm(1, (7.5 + 15) / 2, 10)
#   else y <- rnorm(1, 7.5 / 2, 10)
#   y
# }
# score <- rep(NA, 10000)
# for (i in 1:10000) score[i] <- scoreFn(x[i])
# score <- 100 / (max(score) - min(score)) * (score - max(score)) + 100
# dat <- data.frame(grade = x, score = score)
# ggplot(dat, aes(x = score)) + geom_density() #+ geom_segment(aes(x=95,y=0,xend=100,yend=0.1))
# ggplot() + geom_bar(aes(x = grade, y = prob), stat = "identity")
# qplot(grade, prob, data = dat, geom = "bar", stat = "identity")
