# Chapter 4
# Figures for post-estimation section
# source("00-prelim.R")
library(iprior)
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
chapter.no <- "04"

# Set ggplot theme
theme_reg <- theme_bw() +
  theme(axis.ticks = element_blank(), axis.text = element_blank(),
        panel.grid = element_blank())

# Generate data set. Note that this data set can also be generated using
# iprior::gen_smooth().
set.seed(123)
N <- 150
f <- function(x, truth = FALSE) {
  35 * dnorm(x, mean = 1, sd = 0.8) +
    65 * dnorm(x, mean = 4, sd = 1.5) +
    (x > 4.5) * (exp((1.25 * (x - 4.5))) - 1) +
    3 * dnorm(x, mean = 2.5, sd = 0.3)
}
x <- c(seq(0.2, 1.9, length = N * 5 / 8), seq(3.7, 4.6, length = N * 3 / 8))
x <- sample(x, size = N)
x <- x + rnorm(N, sd = 0.65)  # adding random fluctuation to the x
x <- sort(x)
y.err <- rt(N, df = 1)
y <- f(x) + sign(y.err) * pmin(abs(y.err), rnorm(N, mean = 4.1))  # adding random terms to the y

# True values
x.true <- seq(min(x), max(x), length = 500)
y.true <- f(x.true, TRUE)

# Data for plot
dat <- data.frame(x, y, points = "Observed")
dat.truth <- data.frame(x = x.true, y = y.true)
dat.fit <- rbind(
  dat[, 1:2],
  data.frame(x = x, y = f(x, TRUE))
)

linlab <- paste0("alpha + ", "italic(f[true](x))")
p <- ggplot() +
  geom_line(data = dat.truth, aes(x = x, y = y), col = "red3", size = 1,
            linetype = "solid") +
  geom_point(data = dat, aes(x = x, y = y)) +
  scale_x_continuous(name = expression(italic(x))) +
  scale_y_continuous(name = expression(italic(y))) +
  annotate(geom = "text", x = 5.8, y = y.true[length(y.true)] - 1.5,
           label = linlab, col = "red3", parse = TRUE) +
  theme_bw(); p

## ---- direct ----
mod1 <- iprior(y ~ x, dat.fit, kernel = "fbm", train.samp = 1:N)

## ---- EM ----
mod2 <- iprior(y ~ x, dat.fit, kernel = "fbm", method = "em", train.samp = 1:N,
               control = list(maxit = 1000))

## ---- MCMC ----
stan.iprior.dat <- list(y = y, H = iprior::kern_fbm(x), n = length(y),
                        ytrue = f(x))
stan.iprior.mod <- "
data {
  int<lower=0> n; // number of data
  vector[n] y; // data y
  matrix[n, n] H; // the kernel matrix
  vector[n] ytrue; // true values
}
parameters {
  real alpha;
  real<lower=0> psi;
  real<lower=0> lambda;
}
transformed parameters {
  cov_matrix[n] Vy;
  Vy = psi * (lambda * H) * (lambda * H);
  for (i in 1:n)
    Vy[i, i] = Vy[i, i] + 1 / psi;
}
model {
  matrix[n, n] L_cov;
  L_cov = cholesky_decompose(Vy);
  y ~ multi_normal_cholesky(rep_vector(alpha, n), L_cov);
  lambda ~ normal(0, 10);
  psi ~ normal(0, 10);
}
generated quantities {
  // generate from posterior of y
  vector[n] w;
  matrix[n,n] Vw;
  vector[n] ynew;
  vector[n] muynew;
  matrix[n,n] Vynew;
  real rmse;
  real logLik;
  Vw = inverse(Vy);
  w = psi * lambda * H * inverse(Vy) * (y - alpha);
  muynew = alpha + lambda * H * w;
  Vynew = square(lambda) * H * Vw * H;
  for (i in 1:n)
    Vynew[i, i] = Vynew[i, i] + 1 / psi;
  ynew = multi_normal_rng(muynew, Vynew);
  rmse = sqrt(mean(square(ynew - ytrue)));
  logLik = multi_normal_lpdf(y|rep_vector(alpha, n), Vy);
}
"
# Compile the Stan programme
m <- stan_model(model_code = stan.iprior.mod)
m@model_name <- "iprior.fbm"

# Fit stan model
fit.stan <- sampling(m, data = stan.iprior.dat,
                     pars = c("alpha", "lambda", "psi", "ynew", "rmse", "logLik"),
                     iter = 2000, chains = 8, thin = 1)
print(fit.stan, pars = c("alpha", "lambda", "psi", "rmse", "logLik"))

# Check fit
a <- stan2coda(fit.stan)
fit.ggs <- ggs(a[, c("alpha", "lambda", "psi")])
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
b.psi <- postmean[grep("psi", names(postmean))]
b.alpha.se <- postsd[grep("alpha", names(postsd))]
b.lambda.se <- postsd[grep("lambda", names(postsd))]
b.psi.se <- postsd[grep("psi", names(postsd))]
ynew <- postmean[grep("ynew", names(postmean))]
b.rmse <- sqrt(mean((ynew - f(x)) ^ 2))
# b.dens <- postmean[grep("lp__", names(postmean))] - (N / 2) * log(2 * pi)
b.dens <- postmean[grep("logLik", names(postmean))]

## ---- compare ----
met1 <- c(get_intercept(mod1), sd(dat$y) / sqrt(nrow(dat)),  get_lambda(mod1), get_se(mod1)[1], get_psi(mod1),
          get_se(mod1)[3], logLik(mod1), mod1$time$time, get_prederror(mod1)[2])
met2 <- c(get_intercept(mod2), sd(dat$y) / sqrt(nrow(dat)), get_lambda(mod2), get_se(mod2)[1], get_psi(mod2),
          get_se(mod2)[3], logLik(mod2), mod2$time$time, get_prederror(mod2)[2])
met3 <- c(b.alpha, b.alpha.se, b.lambda, b.lambda.se, b.psi, b.psi.se, b.dens,
          231.762, b.rmse)
tab <- cbind(met1, met2, met3)
colnames(tab) <- c("Direct optim.", "EM alg.", "HMC")
rownames(tab) <- c("Intercept", "se int", "lambda", "se lam", "psi", "se psi",
                   "Log density", "Time (sec)", "RMSE")
tab

## ---- save.plots ----
ggsave("figure/04-example_data.pdf", p, "pdf", width = 6.5, height = 4)
move_fig_to_chapter()
