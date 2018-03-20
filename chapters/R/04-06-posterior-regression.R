# Chapter 4
# Figures for post-estimation section
source("00-prelim.R")
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
x.true <- seq(-2.1, 7, length = 1000)
y.true <- f(x.true, TRUE)

# Data for plot
dat <- data.frame(x, y, points = "Observed")
dat.truth <- data.frame(x = x.true, y = y.true)

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

no.of.draws <- 90

# Based on an fBm RKHS
mod <- iprior(y ~ x, dat, kernel = "fbm", est.hurst = FALSE)
n <- length(y)
H.eta <- iprior::get_kern_matrix(mod)
alpha <- iprior::get_intercept(mod)
psi <- iprior::get_psi(mod)
Vy <- psi * H.eta %*% H.eta + diag(1 / psi, n)
w.hat <- mod$w
# H.star <- iprior::get_kern_matrix(mod, newdata = dat.truth)
H.star <- iprior::kern_fbm(x, x.true, gamma = get_hurst(mod)) * iprior::get_lambda(mod) / 2
# reason for this is to not make the prior paths so wild
y.fitted <- predict(mod, newdata = dat.truth)$y  # predicted values
y.fitted2 <- fitted(mod)$y  # fitted values

# Prior variance for f
Vf.pri <- psi * tcrossprod(H.star)

# Posterior variance for f
Vf.pos <- H.star %*% solve(Vy, t(H.star))

# Prepare random draws from prior and posterior --------------------------------
draw.pri <- t(mvtnorm::rmvnorm(no.of.draws, mean = rep(alpha, 1000),
                               sigma = Vf.pri))
draw.pos <- t(mvtnorm::rmvnorm(40, mean = y.fitted, sigma = Vf.pos))
melted.pos <- melt(data.frame(f = draw.pos, x = x.true), id.vars = "x")
melted.pri <- melt(data.frame(f = draw.pri, x = x.true), id.vars = "x")
melted <- rbind(cbind(melted.pri, type = "Prior"),
                cbind(melted.pos, type = "Posterior"))

# Posterior predictive covariance matrix -------------------------------------
varyprior <- abs(diag(Vf.pri)) + 1 / psi
varystar <- abs(diag(Vf.pos)) + 1 / psi
dat.fit <- data.frame(x.true, y.fitted, sdev = sqrt(varystar),
                      type = "95% credible interval")
dat.f <- rbind(data.frame(x = x.true, y = mean(y), sdev = NA, type = "Prior"),
               data.frame(x = x.true, y = y.fitted, sdev = sqrt(varystar), type = "Posterior"))

# Prepare random draws for posterior predictive checks -----------------------
VarY.hat <- H.eta %*% solve(Vy, H.eta) + diag(1 / psi, nrow(Vy))
ppc <- t(mvtnorm::rmvnorm(no.of.draws, mean = y.fitted2, sigma = VarY.hat))
melted.ppc <- melt(data.frame(x = x, ppc = ppc), id.vars = "x")
melted.ppc <- cbind(melted.ppc, type = "Posterior predictive check")

linlab <- paste0("alpha + ", "italic(f(x))")

## ---- prior.and.posterior ----
p <- ggplot() +
  geom_point(data = dat, aes(x = x, y = y), col = "grey60") +
  scale_x_continuous(breaks = NULL, name = expression(italic(x))) +
  scale_y_continuous(breaks = NULL, name = expression(italic(y))) +
  coord_cartesian(
    xlim = c(min(x.true) + 0.45, max(x.true) - 0.45),
    ylim = c(min(y) - 5, max(y) + 5)
  )

p.prior.post <- p +
  geom_line(data = melted,
            aes(x = x, y = value, group = variable, linetype = "1",
                col = "1", size = "1"), alpha = 0.5) +
  geom_line(data = dat.f,
            aes(x = x, y = y, linetype = "2", col = "2", size = "2")) +
  scale_linetype_manual(name = "", values = c("solid", "dashed"),
                        labels = c("Sample paths", "Mean path")) +
  scale_colour_manual(name = "", values = c("steelblue3", "grey20"),
                      labels = c("Sample paths", "Mean path")) +
  scale_size_manual(name = "", values = c(0.19, 0.8),
                    labels = c("Sample paths", "Mean path")) +
  facet_grid(type ~ .) +
  theme_bw() +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.key.width = unit(3, "line"),
    legend.justification = c(1, 0),
    legend.position = c(1 - 0.001, 0 + 0.001),
    legend.text.align = 0,
    legend.background = element_rect(fill = scales::alpha('white', 0))
  )

## ---- credibility.interval ----
ggplot() +
  scale_x_continuous(breaks = NULL, name = expression(italic(x))) +
  scale_y_continuous(breaks = NULL, name = expression(italic(y))) +
  coord_cartesian(
    xlim = c(min(x.true) + 0.45, max(x.true) - 0.45),
    ylim = c(min(y) - 5, max(y) + 5)
  ) +
  geom_ribbon(data = dat.fit, fill = "grey", alpha = 0.5,
              aes(x = x.true, ymin = y.fitted + qnorm(0.05 / 2) * sdev,
                  ymax = y.fitted + qnorm(1 - 0.05 / 2) * sdev)) +
  geom_point(data = dat, aes(x = x, y = y), col = "grey10", shape = 21) +
  geom_line(data = dat.fit, aes(x = x.true, y = y.fitted, col = "1", linetype = "1"), size = 0.8) +
  geom_line(data = dat.truth, aes(x = x, y = y, col = "2", linetype = "2"), size = 0.8) +
  # labs(title = "Fitted regression line with 95% credibility interval for predicted values") +
  # annotate(geom = "text", x = 7 - 0.5, y = dat.fit$y.fitted[1000] - 1.2,
           # label = linlab, col = "red3", parse = TRUE) +
  scale_colour_manual(name = "", values = c("black", "red3"),
                      labels = c("Estimated", "Truth")) +
  scale_linetype_manual(name = "", values = c("solid", "dashed"),
                        labels = c("Estimated", "Truth")) +
  theme_bw() +
  theme(
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    legend.key.width = unit(3, "line"),
    legend.justification = c(1, 0),
    legend.position = c(1 - 0.001, 0 + 0.001),
    legend.text.align = 0,
    legend.background = element_rect(fill = scales::alpha('white', 0))
  ) -> p.cred; p.cred

## ---- ppc ----
p.ppc <- ggplot() +
  scale_x_continuous(breaks = NULL, name = expression(italic(y))) +
  scale_y_continuous(breaks = NULL) +
  geom_line(data = melted.ppc,
            aes(x = value, group = variable, col = "yrep", size = "yrep"),
            stat = "density", alpha = 0.5) +
  geom_line(data = dat, aes(x = y, col = "y", size = "y"), stat = "density") +
  theme(legend.position = "bottom") +
  scale_colour_manual(
    name = NULL, labels = c("Observed", "Replications"),
    values = c("grey10", "steelblue3")
  ) +
  scale_size_manual(
    name = NULL, labels = c("Observed", "Replications"),
    values = c(1.1, 0.19)
  ) +
  # labs(y = "Density", title = "Posterior predictive density check") +
  theme_bw() +
  theme(legend.position = c(0.9, 0.5))

## ---- hist.x ----
p.hist <- ggplot(dat, aes(x = x)) +
  geom_density(fill = "grey80", col = "grey80") +
  scale_x_continuous(limits = c(min(x.true) + 0.8, max(x.true))) +
  theme_void()
gg.marg.size <- 8
# grid.arrange(p.hist, p.prior.post, ncol = 1, heights = c(1, gg.marg.size))
# p.prior.dens <- ggExtra::ggMarginal(
#   p.prior, margins = "x", size = gg.marg.size, fill = "grey80", col = "grey80"
# )

## ---- save.plots ----
ggsave("figure/04-post_reg_prior_post.pdf",
       grid.arrange(p.hist, p.prior.post, ncol = 1, heights = c(1, gg.marg.size)),
       "pdf", width = 6.5, height = 8)
# ggsave("figure/04-post_reg_posterior.pdf", p.posterior, "pdf", width = 6.5, height = 4)
ggsave("figure/04-post_reg_cred.pdf", p.cred, "pdf", width = 6.5, height = 4)
ggsave("figure/04-post_reg_ppc.pdf", p.ppc, "pdf", width = 6.5, height = 4)
move_fig_to_chapter()
