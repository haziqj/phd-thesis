# Chapter 4
# The I-prior likelihood ridge
library(iprior)
library(plotly)


mod <- iprior(y ~ X, gen_smooth(), kernel = "fbm")

no.points <- 50
x <- log(get_lambda(mod))
x <- seq(x - 20, x + 15, length = no.points)  # lambda
x <- sort(c(x, log(get_lambda(mod))))
y <- log(get_psi(mod))
y <- seq(y - 10, y + 15, length = no.points)  # psi
y <- sort(c(y, log(get_psi(mod))))
tab <- expand.grid(x = x, y = y)

z <- rep(NA, nrow(tab))
for (i in seq_along(z)) {
  z[i] <- logLik(mod, theta = as.numeric(tab[i, ]))
}

tab.loglik <- matrix(z, nrow = no.points + 1, ncol = no.points + 1)
zmin <- -5000
plot_ly(x = x, y = y, z = tab.loglik, zauto = FALSE, zmin = zmin) %>%
  add_surface() %>%
  # add_trace(x = log(get_lambda(mod)), y = log(get_psi(mod)), z = -353.2314,
  #           type = "scatter3d", mode = "markers") %>%
  layout(
    title = "I-prior log-likelihood",
    scene = list(
      xaxis = list(title = "log(lambda)"),
      yaxis = list(title = "log(psi)"),
      zaxis = list(title = "Log-likelihood",  range = c(zmin, max(z)))
    ))


