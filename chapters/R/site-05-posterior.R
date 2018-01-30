source("site-00-prelim.R")

theme_reg <- theme_bw() +
  theme(axis.ticks = element_blank(), axis.text = element_blank(),
        panel.grid = element_blank())

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

no.of.draws <- 100

# Based of an fBm RKHS
mod <- optim(c(1, 1, 0), dev_FBMkern_iprior, method = "L-BFGS", y = y)
n <- length(y)
H <- fnH3(x, gamma = pnorm(mod$par[3])); class(H) <- NULL
alpha <- mean(y)
lambda <- exp(mod$par[1])
psi <- exp(mod$par[2])
Vy <- psi * (lambda * H) %*% (lambda * H) + diag(1 / psi, n)
w.hat <- psi * lambda * H %*% solve(Vy, y - alpha)
H.star <- fnH3(x = x, y = x.true, gamma = pnorm(mod$par[3]))
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

# Posterior predictive covariance matrix -------------------------------------
varyprior <- abs(diag(Vf.pri)) + 1 / psi
varystar <- abs(diag(Vf.pos)) + 1 / psi
dat.fit <- data.frame(x.true, y.fitted, sdev = sqrt(varystar),
                      type = "95% credible interval")
dat.f <- rbind(data.frame(x = x.true, y = mean(y), sdev = NA, type = "Prior"),
               data.frame(x = x.true, y = y.fitted, sdev = sqrt(varystar), type = "Posterior"))

# Prepare random draws for posterior predictive checks -----------------------
VarY.hat <- (lambda ^ 2) * H %*% solve(Vy, H) + diag(1 / psi, nrow(Vy))
ppc <- t(mvtnorm::rmvnorm(no.of.draws, mean = y.fitted2, sigma = VarY.hat))
melted.ppc <- melt(data.frame(x = x, ppc = ppc), id.vars = "x")
melted.ppc <- cbind(melted.ppc, type = "Posterior predictive check")

linlab <- "italic(f(x))"

## ---- prior ----
p <- ggplot() +
  geom_point(data = dat, aes(x = x, y = y), col = "grey60") +
  scale_x_continuous(breaks = NULL, name = expression(italic(x))) +
  scale_y_continuous(breaks = NULL, name = expression(italic(y))) +
  coord_cartesian(
    xlim = c(min(x.true) + 0.45, max(x.true) - 0.45),
    ylim = c(min(y) - 5, max(y) + 5)
  )

p.prior <- p +
  geom_line(data = subset(melted, type == "Prior"),
            aes(x = x, y = value, group = variable, linetype = "1",
                col = "1", size = "1"), alpha = 0.5) +
  geom_line(data = subset(dat.f, type == "Prior"),
            aes(x = x, y = y, linetype = "2", col = "2", size = "2")) +
  scale_linetype_manual(name = "", values = c("solid", "dashed"),
                        labels = c("Sample paths", "Mean path")) +
  scale_colour_manual(name = "", values = c("steelblue3", "grey20"),
                      labels = c("Sample paths", "Mean path")) +
  scale_size_manual(name = "", values = c(0.19, 0.8),
                    labels = c("Sample paths", "Mean path")) +
  labs(title = "Random draws from the fBm-0.5 I-prior") +
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

p.prior.th <- ggplot() +
  geom_point(data = dat, aes(x = x, y = y), col = "grey60", size = 1.1) +
  scale_x_continuous(
    breaks = NULL, name = expression(italic(x))
  ) +
  scale_y_continuous(
    breaks = NULL, name = expression(italic(y))
  ) +
  coord_cartesian(
    xlim = c(min(x.true) + 0.45, max(x.true) - 0.45),
    ylim = c(min(y) - 5, max(y) + 5)
  ) +
  geom_line(data = subset(melted, type == "Prior"),
            aes(x = x, y = value, group = variable),
            col = "steelblue3", size = 0.19, alpha = 0.5) +
  geom_line(data = subset(dat.f, type == "Prior"),
            aes(x = x, y = y), linetype = "dashed", col = "grey20", size = 0.8) +
  theme_void() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

ggsave("image/post_reg_prior_th.png", width = 4, height = 8 / 3, dpi = 150, scale = 0.8)
ggsave("image/post_reg_prior.png", p.prior, "png", width = 9, height = 6)
move_fig_to_site()

## ---- posterior ----
p.posterior <- p +
  geom_line(data = subset(melted, type == "Posterior"),
            aes(x = x, y = value, group = variable, linetype = "1",
                col = "1", size = "1"), alpha = 0.5) +
  geom_line(data = subset(dat.f, type == "Posterior"),
            aes(x = x, y = y, linetype = "2", col = "2", size = "2")) +
  scale_linetype_manual(name = "", values = c("solid", "dashed"),
                        labels = c("Sample paths", "Mean path")) +
  scale_colour_manual(name = "", values = c("steelblue3", "grey20"),
                      labels = c("Sample paths", "Mean path")) +
  scale_size_manual(name = "", values = c(0.19, 0.8),
                    labels = c("Sample paths", "Mean path")) +
  labs(title = "Random draws from the posterior of an fBm-0.5 RKHS under an I-prior") +
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

p.posterior.th <- ggplot() +
  geom_point(data = dat, aes(x = x, y = y), col = "grey60", size = 1.1) +
  scale_x_continuous(breaks = NULL, name = expression(italic(x))) +
  scale_y_continuous(breaks = NULL, name = expression(italic(y))) +
  coord_cartesian(
    xlim = c(min(x.true) + 0.45, max(x.true) - 0.45),
    ylim = c(min(y) - 5, max(y) + 5)
  ) +
  geom_line(data = subset(melted, type == "Posterior"),
            aes(x = x, y = value, group = variable),
            col = "steelblue3", size = 0.19, alpha = 0.5) +
  geom_line(data = subset(dat.f, type == "Posterior"),
            aes(x = x, y = y), linetype = "dashed", col = "grey20", size = 0.8) +
  theme_void() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

ggsave("image/post_reg_posterior_th.png", p.posterior.th, "png", width = 4, height = 8 / 3, dpi = 150, scale = 0.8)
ggsave("image/post_reg_posterior.png", p.posterior, "png", width = 9, height = 6)
move_fig_to_site()

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
  geom_point(data = dat, aes(x = x, y = y), col = "grey10") +
  geom_line(data = dat.fit, aes(x = x.true, y = y.fitted), col = "red3",
            size = 0.8) +
  labs(title = "Fitted regression line with 95% credibility interval for predicted values") +
  annotate(geom = "text", x = 7 - 0.3, y = dat.fit$y.fitted[1000] - 1,
           label = linlab, col = "red3", parse = TRUE) +
  theme_reg -> p.cred

ggplot() +
  scale_x_continuous(breaks = NULL, name = expression(italic(x))) +
  scale_y_continuous(breaks = NULL, name = expression(italic(y))) +
  coord_cartesian(
    xlim = c(min(x.true) + 0.45, max(x.true) - 0.45),
    ylim = c(min(y) - 5, max(y) + 5)
  ) +
  geom_ribbon(data = dat.fit, fill = "grey", alpha = 0.5,
              aes(x = x.true, ymin = y.fitted - 1.96 * sdev,
                  ymax = y.fitted + 1.96 * sdev)) +
  geom_point(data = dat, aes(x = x, y = y), col = "grey10", size = 1.1) +
  geom_line(data = dat.fit, aes(x = x.true, y = y.fitted), col = "red3",
            size = 0.8) +
  theme_void() -> p.cred.th

ggsave("image/post_reg_cred_th.png", p.cred.th, "png", width = 4, height = 8 / 3, dpi = 150, scale = 0.8)
ggsave("image/post_reg_cred.png", p.cred, "png", width = 9, height = 6)
move_fig_to_site()

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
  labs(y = "Density", title = "Posterior predictive density check") +
  theme_bw() +
  theme(legend.position = c(0.9, 0.5))

p.ppc.th <- ggplot() +
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
  theme_void() +
  theme(legend.position = "none")

ggsave("image/post_reg_ppc_th.png", p.ppc.th, "png", width = 4, height = 8 / 3, dpi = 150, scale = 0.8)
ggsave("image/post_reg_ppc.png", p.ppc, "png", width = 9, height = 6)
move_fig_to_site()












