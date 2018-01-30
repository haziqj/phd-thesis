# Chapter 4, Section 2
# Multilevel modelling using I-priors

## ---- 04_01_data ----
data(tutorial, package = "R2MLwiN")
str(tutorial[, c("normexam", "school", "standlrt")])

## ---- 04_01_analysis ----
# as_tibble(tutorial) %>%
#   select(normexam, school, standlrt) -> tutorial.tb
# tutorial.tb %>%
#   group_by(school) %>%
#   summarise(mean_normexam = mean(normexam), mean_standlrt = mean(standlrt))
school_levels <- levels(tutorial$school)
school_unique <- match(school_levels, tutorial$school)

# Varying intercept model
# MLE value obtained by running the EM a total of 8 hours
# 120 iterations to complete, max log-lik = -5503.8503
mod1.fit <- iprior(normexam ~ school, data = tutorial,
                   control = list(lambda = 0.0006998815, psi = 1.1799067261,
                                  maxit = 1, silent = TRUE))

# Standard random effects model
mod1a <- lmer(normexam ~ 1 + (1 | school), tutorial)

# Obtain fitted values
mod_vi_em_unique_fitted <- fitted(mod1.fit)[school_unique]
mod_vi_re_unique_fitted <- coef(mod1a)$school[, 1]

# The two peculiar values school 48 and 54
int1 <- decPlac(mod_vi_em_unique_fitted[48])
int2 <- decPlac(mod_vi_em_unique_fitted[54])
int1a <- decPlac(mod_vi_re_unique_fitted[48])
int2a <- decPlac(mod_vi_re_unique_fitted[54])

# tibble(
#   school = factor(school_levels),
#   iprior = mod_vi_em_unique_fitted,
#   re     = mod_vi_re_unique_fitted
# ) -> mod_vi_comparison
# mod_vi_comparison %>%
#   mutate(diff = abs(iprior - re)) %>%
#   arrange(desc(diff)) %>%
#   select(school) -> top_two

# Varying slopes model
# MLE value obtained by running the EM a total of 24 hours!
# 635 iterations to complete, max log-lik = -4670.3808
mod2.fit <- iprior(normexam ~ school * standlrt, data = tutorial,
                   control = list(lambda = c(0.0004234420, 0.37315515627),
                                  psi = 1.8028197426, maxit = 1, silent = TRUE))

# Standard random effects model
mod2a <- lmer(normexam ~ 1 + standlrt + (1 + standlrt | school), tutorial)

# Obtain slopes
slopes2a <- coef(mod2a)$school[, 2]
slopes2 <- slope(mod2.fit)
indslop <- sort(as.numeric(names(slopes2)), index.return = TRUE)$ix
slopes2 <- slopes2[indslop]

# Capture the output
mod1.output <- capture.output(print(mod1.fit))
mod2.output <- capture.output(print(mod2.fit))

# Delete models otherwise take up too much space in cache
rm(mod1.fit)
rm(mod2.fit)
rm(mod1a)
rm(mod2a)

## ---- 04_01_plot_em ----
tibble(
  school = factor(school_levels),
  iprior = mod_vi_em_unique_fitted,
  re     = mod_vi_re_unique_fitted
) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1, col = "grey30", linetype = "dashed") +
  geom_text(aes(x = re, y = iprior, label = school, col = school)) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "Standard random effects model intercepts",
       y = "I-prior intercepts",
       title = "Varying intercept model")

tibble(
  school = factor(school_levels),
  iprior = slopes2,
  re     = slopes2a
) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1, col = "grey30", linetype = "dashed") +
  geom_text(aes(x = re, y = iprior, label = school, col = school)) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "Standard random effects model slopes",
       y = "I-prior slopes",
       title = "Varying slopes model")

## ---- 04_01_nystrom ----
mod_vi_nys <- ipriorNystrom(normexam ~ school, data = tutorial, size = 250)
mod_vs_nys <- ipriorNystrom(normexam ~ school * standlrt, data = tutorial,
                            size = 400)

## ---- 04_01_nystorm_vi ----
# Varying intercept model comparison
tibble(
  school = factor(school_levels),
  iprior = fitted(mod_vi_nys)[school_unique],
  re = coef(mod1a)[[1]][, 1]
) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1, col = "grey30", linetype = "dashed") +
  geom_text(aes(x = re, y = iprior, label = school, col = school)) +
  theme_bw() +
  theme(legend.position = "none")

## ---- 04_01_nystorm_vs ----
# Varying slopes model comparison
tibble(
  school = factor(tutorial$school),
  x = tutorial$standlrt,
  y = fitted(mod_vs_nys)
) %>%
  group_by(school) %>%
  arrange(x, .by_group = TRUE) %>%
  summarise(iprior = (sum(x * y) - sum(x) * sum(y) / length(y)) / sum((x - mean(x)) ^ 2)) %>%
  mutate(re = coef(mod2a)$school[, 2]) %>%
  ggplot() +
  geom_abline(intercept = 0, slope = 1, col = "grey30", linetype = "dashed") +
  geom_text(aes(x = re, y = iprior, label = school, col = school)) +
  theme_bw() +
  theme(legend.position = "none")
