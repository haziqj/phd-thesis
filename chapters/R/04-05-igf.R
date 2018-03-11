# Chapter 4
# Random effects model (IGF data set)
source("00-prelim.R")

## ---- IGF.data ----
data(IGF, package = "nlme")
head(IGF)

## ---- IGF.mod.iprior ----
mod.iprior <- iprior(conc ~ age * Lot, IGF, method = "em")
summary(mod.iprior)

## ---- IGF.mod.iprior.plot ----
plot_fitted_multilevel(mod.iprior, facet = 1, cred.bands = FALSE,
                       extrapolate = TRUE, show.legend = FALSE)

## ---- IGF.mod.iprior.const ----
mod.iprior.const <- iprior(conc ~ age * Lot, IGF, method = "em",
                           est.lambda = FALSE, lambda = c(0, 0))
logLik(mod.iprior.const)
(D <- -2 * (logLik(mod.iprior.const) - logLik(mod.iprior)))
pchisq(D, df = 2)

## ---- IGF.mod.lmer ----
(mod.lmer <- lmer(conc ~ age + (age | Lot), IGF))
eigen(VarCorr(mod.lmer)$Lot)

## ---- IGf.prep.plot ----
grp <- unique(as.numeric(IGF$Lot))
beta.iprior <- matrix(NA, nrow = length(grp), ncol = 2)
tmp.df <- data.frame(x = IGF$age, y = fitted(mod.iprior)$y,
                     grp = as.numeric(IGF$Lot))
for (i in seq_along(grp)) {
  beta.iprior[i, ] <- coef(lm(y ~ x, tmp.df[tmp.df$grp == grp[i], ]))
}
beta.lmer <- coef(mod.lmer)$Lot

beta.fixed.iprior <- coef(lm(y ~ x, tmp.df))
beta.fixed.lmer <- mod.lmer@beta
beta.fixed.df <- data.frame(
  beta = c(beta.fixed.iprior, beta.fixed.lmer),
  type = rep(c("Intercept", "Slope")),
  model = rep(c("iprior", "lmer"), each = 2)
)

Sigma.iprior <- cov(beta.iprior)
sigma0.iprior <- sqrt(Sigma.iprior[1, 1])
sigma1.iprior <- sqrt(Sigma.iprior[2, 2])
corr.iprior <- cor(beta.iprior)[1, 2]  # Sigma.iprior[1, 2] / (sigma0.iprior * sigma1.iprior)
Sigma.lmer <- VarCorr(mod.lmer)[[1]]
sigma0.lmer <- attr(Sigma.lmer, "stddev")[1]
sigma1.lmer <- attr(Sigma.lmer, "stddev")[2]
corr.lmer <- attr(Sigma.lmer, "correlation")[1, 2]

plot.df.beta <- data.frame(
  beta = as.numeric(beta.iprior),
  Lot = unique(IGF$Lot),
  type = rep(c("Intercept", "Slope"), each = 10),
  model = "iprior"
)
plot.df.beta <- rbind(plot.df.beta, data.frame(
  beta = unlist(beta.lmer),
  Lot = unique(IGF$Lot),
  type = rep(c("Intercept", "Slope"), each = 10),
  model = "lmer"
))

plot.df.param <- data.frame(
  param = c(sigma0.iprior, sigma1.iprior, corr.iprior,
            sigma0.lmer, sigma1.lmer, corr.lmer),
  name = rep(c("sigma[0]", "sigma[1]", "sigma[01]"), 2),
  model = rep(c("iprior", "lmer"), each = 3)
)
plot.df.param$name <- factor(plot.df.param$name,
                             levels = c("sigma[01]", "sigma[1]", "sigma[0]"))

## ---- IGF.plot.beta ----
plot.df.fake <- plot.df.beta[c(19, 39, 9, 29), ]
plot.df.fake[1, 1] <- 0.02  # slopes
plot.df.fake[2, 1] <- -0.02
plot.df.fake[3, 1] <- 5.45  # intercepts
plot.df.fake[4, 1] <- 5.25
ggplot(plot.df.beta) +
  geom_point(aes(beta, Lot, col = model)) +
  geom_point(data = plot.df.fake, aes(beta, Lot, col = model), alpha = 0) +
  geom_vline(data = beta.fixed.df, aes(xintercept = beta, col = model),
             linetype = "dashed") +
  facet_grid(. ~ type, scales = "free") +
  theme_bw() +
  theme(legend.position = "top")

## ---- IGF.plot.param ----
ggplot(plot.df.param) +
  geom_point(aes(name, param, col = model),
             position = position_dodge(width = 0.2)) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  theme_bw()
