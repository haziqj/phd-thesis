source("site-00-prelim.R")

## ---- points ----
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
x.true <- seq(-2.1, 7, length = 1000)
y.true <- f(x.true, TRUE)

# # New points
# set.seed(123)
# x.up <- runif(10, max(x), max(x.true))
# x.lo <- runif(10, min(x.true), min(x))
# x.new <- c(x.up, x.lo)
# y.err.new <- rt(10, df = 1)
# y.new <- f(x.new) + sign(y.err.new) * pmin(abs(y.err.new), rnorm(10, mean = 4.1))

# Data for plot
dat <- data.frame(x, y, points = "Observed")
# dat <- rbind(dat, data.frame(x = x.new, y = y.new, points = "Unobserved"))
dat.truth <- data.frame(x.true, y.true)

p1 <- ggplot() +
  geom_point(data = dat, aes(x = x, y = y)) +
  scale_x_continuous(
    limits = c(min(x.true), max(x.true)),
    breaks = NULL, name = expression(italic(x))
  ) +
  scale_y_continuous(
    breaks = NULL, name = expression(italic(y))
  ) +
  coord_cartesian(ylim = c(min(y, y) - 5, max(y, y) + 5)) +
  theme_bw()

## ---- main.plot ----
no.of.draws <- 190

# Based of an fBm-0.5 RKHS
mod <- optim(c(1, 1, 0), dev_FBMkern_iprior, method = "L-BFGS", y = y)
n <- length(y)
H <- fnH3(x, gamma = expit(mod$par[3])); class(H) <- NULL
alpha <- mean(y)
lambda <- exp(mod$par[1])
psi <- exp(mod$par[2])
Vy <- psi * (lambda * H) %*% (lambda * H) + diag(1 / psi, n)
w.hat <- psi * lambda * H %*% solve(Vy, y - alpha)
H.star <- fnH3(x = x, y = x.true, gamma = expit(mod$par[3]))
class(H.star) <- NULL
y.fitted <- as.numeric(mean(y) + lambda * H.star %*% w.hat)
y.fitted2 <- as.numeric(mean(y) + lambda * H %*% w.hat)

# Prior variance for f
Vf.pri <- psi * lambda ^ 2 * tcrossprod(H.star)
class(Vf.pri) <- NULL

# Posterior variance for f
Vf.pos <- lambda ^ 2 * H.star %*% solve(Vy, t(H.star))
class(Vf.pos) <- NULL

# Prepare random draws from prior and posterior --------------------------------
draw.pri <- t(mvtnorm::rmvnorm(no.of.draws, mean = rep(alpha, 1000),
                               sigma = Vf.pri))
draw.pos <- t(mvtnorm::rmvnorm(no.of.draws, mean = y.fitted, sigma = Vf.pos))
melted.pos <- melt(data.frame(f = draw.pos, x = x.true), id.vars = "x")
melted.pri <- melt(data.frame(f = draw.pri, x = x.true), id.vars = "x")
melted <- rbind(cbind(melted.pri, type = "Prior"),
                cbind(melted.pos, type = "Posterior"))

p.main <- ggplot() +
  geom_point(data = dat, aes(x = x, y = y), col = "grey50") +
  scale_x_continuous(
    limits = c(min(x.true), max(x.true)),
    breaks = NULL, name = expression(italic(x))
  ) +
  scale_y_continuous(
    breaks = NULL, name = expression(italic(y))
  ) +
  geom_line(data = subset(melted, type == "Prior"),
            aes(x = x, y = value, group = variable),
            col = "steelblue3", size = 0.33, alpha = 0.5) +
  coord_cartesian(ylim = c(min(y) + 5, max(y) - 2),
                  xlim = c(min(dat$x), max(dat$x))) +
  theme_void() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
# p.main

ggsave("image/main_image.png", p.main, "png", width = 12, height = 4, dpi = 150)
move_fig_to_site()
