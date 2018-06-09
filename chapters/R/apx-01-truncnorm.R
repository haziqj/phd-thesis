library(iprobit)
library(plotly)
library(tidyverse)

no.points <- 100
x <- y <- seq(-3, 3, length = no.points)
tab <- expand.grid(x = x, y = y, z = NA)
tab[, 3] <- mvtnorm::dmvnorm(tab[, 1:2], mean = c(0, 0)) * (tab[, 1] > tab[, 2])
z <- matrix(tab[, 3], nrow = no.points, ncol = no.points)

plot_ly(x = x, y = y, z = z) %>%
  add_surface() %>%
  layout(
    title = "Truncated multivariate normal density surface",
    scene = list(
      xaxis = list(title = "x2"),
      yaxis = list(title = "x2"),
      zaxis = list(title = "p(x1, x2)")
    ))

ggplot(tab, aes(x, y, z = z)) + geom_contour()
