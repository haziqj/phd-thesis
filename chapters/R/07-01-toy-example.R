# Chapter 7, Section 1
# Illustrating iprobit models with toy examples
source("00-prelim.R")

## ---- spiral_data ----
spiral <- gen_spiral(n = 300, sd = 0.07)
plot(spiral)

## ---- spiral_canonical ----
# Bad results, linear functions not able to predict spirals well
(mod1 <- iprobit(y ~ X1 + X2, spiral, kernel = "Canonical"))
iplot_predict(mod1)

## ---- spiral_fbm ----
# Getting there, but still not nice
(mod2 <- iprobit(y ~ X1 + X2, spiral, kernel = "FBM"))
iplot_predict(mod2)

## ---- spiral_fbm_onelam ----
# Turns out the scale parameters matter here
(mod3 <- iprobit(y ~ X1 + X2, spiral, kernel = "FBM", one.lam = TRUE))
iplot_predict(mod3)

## ---- circle_data ----
circle <- gen_circle(n = 500, m = 4)
plot(circle)

## ---- circle_mod ----
(mod <- iprobit(y ~ X1 + X2, circle, kernel = "FBM", one.lam = TRUE))
iplot_predict(mod)

## ---- mixture_mod ----
mixture <- gen_mixture(n = 500, m = 4, sd = 1.5)
(mod <- iprobit(y ~ X1 + X2, mixture))

## ---- mixture_mod_plot ----
plot(mixture)
iplot_predict(mod)
