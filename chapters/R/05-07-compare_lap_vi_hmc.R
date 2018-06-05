# Chapter 7
# Comparing VI, Laplace and HMC
library(iprobit)
library(parallel)
library(doSNOW)
library(tidyverse)
library(plotly)
library(ggmcmc)
library(coda)
library(mvtnorm)
library(MCMCglmm)  # for posterior.mode()
library(runjags)  # for combine.mcmc()
library(rstan)
stan2coda <- function(fit) {
  mcmc.list(lapply(1:ncol(fit), function(x) mcmc(as.array(fit)[,x,])))
}
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
chapter.no <- "05"

## ---- spiral.data ----
dat <- iprobit::gen_spiral(n = 300, seed = 123)  # generate binary toy example data set

## ---- iprobit.vi ----
mod.vi <- iprobit(y ~ X1 + X2, dat, one.lam = TRUE, kernel = "fbm")

## ---- iprobit.lap ----
# mod.lap <- iprobit(y ~ X1 + X2, dat, one.lam = TRUE, kernel = "fbm",
#                    method = "laplace")
# save(mod.lap, file = "compare_lap")
load("data/compare_lap")

## ---- iprobit.hmc ----
stan.iprobit.dat <- list(y = as.numeric(dat$y) - 1, H = iprior::kern_fbm(dat$X),
                         n = length(dat$y))
stan.iprobit.mod <- "
data {
  int<lower=0> n; // number of data
  int<lower=0,upper=1> y[n]; // data y
  matrix[n, n] H; // the kernel matrix
}
parameters {
  real alpha;
  real<lower=0> lambda;
  vector[n] w;
}
transformed parameters {
  vector[n] mu;
  vector<lower=0,upper=1>[n] pi;
  mu = alpha + lambda * H * w;
  pi = Phi(mu);
}
model {
  y ~ bernoulli(pi);
  w ~ normal(0, 1);
  lambda ~ normal(0, 10);
}
generated quantities {
  // generate from posterior of y
  vector[n] ypred;
  vector[n] yvec;
  real brierscore;
  real errorrate;
  real logLik;
  for (i in 1:n)
    yvec[i] = y[i];
  brierscore = mean(square(pi - yvec));
  for (i in 1:n)
    if (mu[i] > 0) {
      ypred[i] = 1;
    } else {
      ypred[i] = 0;
    }
  errorrate = mean(square(ypred - yvec)) * 100;
  logLik = bernoulli_lpmf(y|pi);
}
"
# Compile the Stan programme
m <- stan_model(model_code = stan.iprobit.mod)
m@model_name <- "iprobit.fbm"

# Fit stan model
fit.stan <- sampling(m, data = stan.iprobit.dat,
                     pars = c("alpha", "lambda", "brierscore", "errorrate", "w",
                              "logLik"),
                     iter = 200, chains = 8, thin = 1)
print(fit.stan, pars = c("alpha", "lambda", "brierscore", "errorrate", "logLik"))

# Check fit
a <- stan2coda(fit.stan)
fit.ggs <- ggs(a[, c("alpha", "lambda")])
ggs_density(fit.ggs)
ggs_traceplot(fit.ggs)
ggs_running(fit.ggs)
ggs_autocorrelation(fit.ggs)
ggs_compare_partial(fit.ggs)
ggs_geweke(fit.ggs)
ggs_caterpillar(fit.ggs)

postmean <- summary(stan2coda(fit.stan))$stat[, 1]
postsd <- summary(stan2coda(fit.stan))$stat[, 2]
b.alpha <- postmean[grep("alpha", names(postmean))]
b.lambda <- postmean[grep("lambda", names(postmean))]
b.alpha.se <- postsd[grep("alpha", names(postsd))]
b.lambda.se <- postsd[grep("lambda", names(postsd))]
b.w <- postmean[grep("w", names(postmean))]
b.w.se <- postmean[grep("w", names(postsd))]

# Slight hack
mod.hmc <- iprobit(y ~ X1 + X2, dat, one.lam = TRUE, kernel = "fbm",
                   control = list(theta0 = log(b.lambda), int.only = TRUE))
mod.hmc$w <- b.w
mod.hmc$param.full[1, ] <- b.alpha
iplot_predict(mod.hmc)

## ---- compare.lik ----
no.points <- 50
x <- log(get_lambda(mod.vi))
x <- seq(x - 5, x + 5, length = no.points)
y <- get_alpha(mod.vi)
y <- seq(y - 5, y + 5, length = no.points)
tab <- expand.grid(x = x, y = y)

# # Generate z points (takes a long time)
stan.iprobit.mod <- "
data {
int<lower=0> n; // number of data
int<lower=0,upper=1> y[n]; // data y
matrix[n, n] H; // the kernel matrix
real alpha;
real<lower=0> lambda;
}
parameters {
vector[n] w;
}
transformed parameters {
vector<lower=0,upper=1>[n] pi;
pi = Phi(alpha + lambda * H * w);
}
model {
w ~ normal(0, 1);
y ~ bernoulli(pi);
}
generated quantities {
  // generate from posterior of y
real logLik;
logLik = bernoulli_lpmf(y|pi);
}
"
m <- stan_model(model_code = stan.iprobit.mod)
pre.dat <- list(y = as.numeric(dat$y) - 1, H = iprior::kern_fbm(dat$X),
                n = length(dat$y))
stan.iprobit.dat <- c(pre.dat, b.alpha, b.lambda)
fit.stan <- sampling(
  m, data = stan.iprobit.dat, iter = 60, chains = 1, pars = c("logLik"),
  thin = 1, init = list(list(w = rep(0, 300)))
)
# postmean <- summary(stan2coda(fit.stan))$stat[, 1]
# tmp <- data.frame(1 - postmean[-301], postmean[-301])
# colnames(tmp) <- levels(dat$y)
# loglik <- with(tmp, {sum(log(tmp[cbind(seq_along(dat$y), dat$y)]))})

do_it <- function() {
  pb <- txtProgressBar(min = 0, max = nrow(tab), style = 3)
  progress <- function(i) setTxtProgressBar(pb, i)

  cl <- makeCluster(parallel::detectCores())
  registerDoSNOW(cl)
  res <- foreach(i = 1:1250, .combine = rbind,  # seq_len(nrow(tab))
                 .packages = c("iprobit", "rstan", "coda", "purrr"),
                 .export = c("mod.vi", "tab", "m", "pre.dat", "stan2coda"),
                 .options.snow = list(progress = progress)) %dopar% {
                   aa <- tab[i, 2]
                   tt <- tab[i, 1]
                   stan.iprobit.dat <- c(pre.dat, alpha = aa, lambda = exp(tt))
                   fit.stan <- sampling(
                     m, data = stan.iprobit.dat, iter = 60, chains = 1,
                     pars = "logLik", thin = 1,
                     init = list(list(w = rep(0, 300)))
                    )
                   postmean <- summary(stan2coda(fit.stan))$stat[, 1]
                   res.hmc <- postmean[grep("logLik", names(postmean))]
                   # res.var <- logLik(mod.vi, theta = tt, alpha = aa)
                   # res.lap <- lap_bin(mu = c(aa, tt), object = mod.vi$ipriorKernel)
                   # c(res.var, res.lap, res.hmc)
                   res.hmc
                 }
  close(pb)
  stopCluster(cl)
  res
}
z3 <- do_it()
# save(z, file = "data/iprobit_lik")
load("data/iprobit_lik")

# 2D plots
# ggplot(cbind(tab, z = z), aes(x = x, y = y, z = z)) +
#   geom_contour(aes(col = ..level..), binwidth = 1) +
#   theme_bw()

# 3D plots
tab.var <- matrix(z[, 1], nrow = no.points, ncol = no.points)
tab.lap <- matrix(z[, 2], nrow = no.points, ncol = no.points)
tab.hmc <- matrix(z[, 3], nrow = no.points, ncol = no.points)

p1 <- plot_ly(x = x, y = y, z = tab.var) %>%
  add_surface() %>%
  layout(
    title = "I-probit variational likelihood",
    scene = list(
      xaxis = list(title = "log(lambda)"),
      yaxis = list(title = "alpha"),
      zaxis = list(title = "Log-likelihood")
    ))
p2 <- plot_ly(x = x, y = y, z = tab.lap, zauto = FALSE, zmin = min(z[, 1])) %>%
  add_surface() %>%
  layout(
    title = "I-probit Laplace likelihood",
    scene = list(
      xaxis = list(title = "log(lambda)"),
      yaxis = list(title = "alpha"),
      zaxis = list(title = "Log-likelihood", range = c(min(z[, 1]), max(z[, 2])))
    ))
p3 <- plot_ly(x = x, y = y, z = tab.hmc, zauto = FALSE, zmin = min(z[, 1])) %>%
  add_surface() %>%
  layout(
    title = "I-probit MCMC likelihood",
    scene = list(
      xaxis = list(title = "log(lambda)"),
      yaxis = list(title = "alpha"),
      zaxis = list(title = "Log-likelihood")
    ))

## ---- save.plots ----
p1f <- iplot_predict(mod.vi)
p2f <- iplot_predict(mod.lap)
p3f <- iplot_predict(mod.hmc)

ggsave("figure/05-example_data.pdf", plot(dat), "pdf", width = 6, height = 4.5)
ggsave("figure/05-fit_lap.pdf",
       p2f + ggtitle("(a) Laplace approximation"), "pdf",
       width = 6.5 * 0.9, height = 4.5)
ggsave("figure/05-fit_vi.pdf",
       p1f + ggtitle("(b) Variational EM"), "pdf",
       width = 6.5 * 0.9, height = 4.5)
ggsave("figure/05-fit_hmc.pdf",
       p3f + ggtitle("(c) Hamiltonian MC"), "pdf",
       width = 6.5 * 0.9, height = 4.5)
move_fig_to_chapter()
