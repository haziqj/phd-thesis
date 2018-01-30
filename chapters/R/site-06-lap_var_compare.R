source("site-00-prelim.R")

# The pdf to approximate
rx <- function(x) {
  exp(-(x ^ 2) / 2) / (1 + exp(-(20 * x + 4)))
}
const <- 1 / integrate(rx, -Inf, Inf)$value
px <- function(x) {
  rx(x) * const
}
dev <- function(x) -2 * log(px(x))
EX <- integrate(function(x) x * px(x), -Inf, Inf)$value
EX2 <- integrate(function(x) (x ^ 2) * px(x), -Inf, Inf)$value
VarX <- EX2 - EX ^ 2
x.hat <- optim(0, dev, method = "BFGS")$par

# Laplace approximation
tmp <- optim(0, function(x) -log(px(x)), method = "BFGS", hessian = TRUE)
mode.p <- tmp$par
var.p <- as.numeric(1 / tmp$hessian)
lap.approx <- function(x) dnorm(x, mean = mode.p, sd = sqrt(var.p))
lap.approx.dev <- function(x) -2 * log(lap.approx(x))
x.lap <- optim(0, lap.approx.dev, method = "BFGS")$par

# Variational approximation
# Minimise KL(q||p) = E_q[log q(x) - log p(x)]
KL.qp <- function(theta, lim = 23.5947161) {  #23.5947161
  mu <- theta[1]
  sigma <- theta[2]
  log.qx <- function(x, log.qx = TRUE) dnorm(x, mean = mu, sd = sigma, log = log.qx)
  log.px <- function(x) log(px(x))
  integrand <- function(x) {
    (log.qx(x) - log.px(x)) * dnorm(x, mean = mu, sd = sigma)
  }
  integrate(integrand, -lim, lim)$value
}
res <- optim(c(0, 1), KL.qp, method = "L-BFGS-B", lower = c(-Inf, 1e-12))
mu.var <- res$par[1]
sd.var <- res$par[2]
var.approx <- function(x) dnorm(x, mean = mu.var, sd = sd.var)
var.approx.dev <- function(x) -2 * log(var.approx(x))
x.var <- optim(0, var.approx.dev, method = "BFGS")$par

# Plot the densities
x <- seq(-2, 3, length = 1000)
den.df <- data.frame(x = x, Truth = px(x), Laplace = lap.approx(x),
                     Variational = var.approx(x))
den.df <- melt(den.df, id.vars = "x")
max.df <- data.frame(x        = c(x.hat, x.lap, x.var),
                     variable = c("Truth", "Laplace", "Variational"),
                     value    = c(px(x.hat), lap.approx(x.lap), var.approx(x.var)))
dl.list <- list("top.bumptwice",
                dl.move("Truth", -0.21, 0.67),
                dl.move("Laplace", -0.52, 0.55),
                dl.move("Variational", 1.4, 0.6))

ggplot(data = den.df, aes(x = x, y = value, group = variable)) +
  geom_line(aes(col = variable)) +
  geom_vline(xintercept = x.hat, linetype = 2, col = thecol[1]) +
  annotate("text", label = "mode", col = iprior::ggColPal(1), x = x.hat - 0.2, y = -0.03) +
  geom_vline(xintercept = EX, linetype = 2, col = iprior::ggColPal(1)) +
  annotate("text", label = "mean", col = thecol[1], x = EX + 0.2, y = -0.03) +
  geom_area(aes(fill = variable), alpha = 0.3, position = "identity") +
  geom_dl(aes(label = variable, col = variable), method = dl.list) +
  scale_x_continuous(
    breaks = NULL, name = expression(italic(z))
  ) +
  scale_y_continuous(
    limits = c(-0.03, max(px(x.hat), lap.approx(x.lap), var.approx(x.var)) + 0.05),
    breaks = NULL, name = "Density"
  ) +
  labs(title = "Comparison of the Laplace and variational approximation") +
  theme_bw() +
  theme(legend.position = "none") -> p.lap.var

dl.list <- list("top.bumptwice",
                dl.move("Truth", -0.33, 0.67),
                dl.move("Laplace", -0.68, 0.55),
                dl.move("Variational", 1.63, 0.6))

ggplot(data = den.df, aes(x = x, y = value, group = variable)) +
  geom_line(aes(col = variable)) +
  geom_area(aes(fill = variable), alpha = 0.3, position = "identity") +
  geom_dl(aes(label = variable, col = variable), method = dl.list) +
  theme_void() +
  theme(legend.position = "none") -> p.lap.var.th

ggsave("image/compare_lap_var_th.png", p.lap.var.th, "png", width = 4, height = 8 / 3, dpi = 150, scale = 1.1)
ggsave("image/compare_lap_var.png", p.lap.var, "png", width = 9, height = 6)
move_fig_to_site()
