source("site-00-prelim.R")

theme_reg <- theme_bw() +
  theme(axis.ticks = element_blank(), axis.text = element_blank(),
        panel.grid = element_blank(), legend.position = "none")

## ---- linear.regression ----
set.seed(789)
n <- 150
x <- seq(-1, 1, length = n)
y <- x * 2 + rt(n, df = 5)
f <- fitted(iprior(y, x))
plot.lin.df <- data.frame(x, y, f)
linlab <- "italic(f(x))"
p <- ggplot(plot.lin.df, aes(x, y)) +
  geom_point() +
  geom_line(aes(x, f), col = "red3", size = 0.9) +
  labs(x = expression(italic(x)), y = expression(italic(y)))
p.lin.th <- p +
  annotate(geom = "text", x = max(x) + 0.1, y = max(f) + 0.1,
           label = linlab, col = "red3", parse = TRUE) +
  theme_void(); p.lin.th
p.lin <- p +
  annotate(geom = "text", x = max(x) + 0.04, y = max(f) + 0.04,
           label = linlab, col = "red3", parse = TRUE) +
  theme_reg + labs(title = "Fitted regression line from the linear RKHS"); p.lin

ggsave("image/reg_lin_th.png", p.lin.th, "png", width = 4, height = 8 / 3, dpi = 150, scale = 1.1)
ggsave("image/reg_lin.png", p.lin, "png", width = 9, height = 6)
move_fig_to_site()

## ---- smooth.regression ----
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
mod <- iprior(y, x, model = list(kernel = "FBM"))
x.pred <- seq(min(x), max(x), length = 1000)
f <- predict(mod, list(matrix(x.pred)))
linlab <- "italic(f(x))"
p <- ggplot() +
  geom_point(aes(x = x, y = y)) +
  geom_line(aes(x = x.pred, y = f), col = "red3", size = 0.8) +
  labs(x = expression(italic(x)), y = expression(italic(y)))
p.fbm.th <- p +
  annotate(geom = "text", x = x[N] + 0.25, y = f[1000],
           label = linlab, col = "red3", parse = TRUE) +
  theme_void(); p.fbm.th
p.fbm <- p +
  annotate(geom = "text", x = x[N] + 0.15, y = f[1000],
           label = linlab, col = "red3", parse = TRUE) +
  theme_reg + labs(title = "Fitted regression line from the fBm-0.5 RKHS"); p.fbm

ggsave("image/reg_fbm_th.png", p.fbm.th, "png", width = 4, height = 8 / 3, dpi = 150, scale = 1.1)
ggsave("image/reg_fbm.png", p.fbm, "png", width = 9, height = 6)
move_fig_to_site()

## ---- multilevel.regression ----
set.seed(4567)
n <- 25
m <- 6
x <- seq(-1, 1, length = n)
beta <- mvtnorm::rmvnorm(m, c(0, 0), sigma = matrix(c(50, 2, 2, 8), nrow = 2))
beta[6, 2] <- -beta[6, 2]
dat <- as.data.frame(matrix(NA, nrow = n * m, ncol = 3))
for (j in seq_len(m)) {
  x <- seq(-1, 1, length = n) + rnorm(1)
  y <- beta[j, 1] + beta[j, 2] * x + rnorm(n, sd = 2)
  dat[(1:n) + n * (j - 1), ] <- cbind(y, x, j)
}
names(dat) <- c("y", "x", "grp")
dat$grp <- factor(dat$grp)
x.new <- c(min(dat$x), max(dat$x))
dat.new <- as.data.frame(matrix(NA, nrow = 2 * m, ncol = 3))
dat.new[, 2] <- rep(x.new, m)
dat.new[, 3] <- rep(1:m, each = 2)
dat.new[, 1] <- 1
names(dat.new) <- c("y", "x", "grp")
dat.new$grp <- factor(dat.new$grp)

# Random intercept
mod1 <- iprior(y ~ ., dat)
f <- predict(mod1, dat.new)
dat.new1 <- cbind(dat.new, f)
p <- ggplot(dat) +
  geom_point(aes(x, y, col = grp)) +
  geom_line(data = dat.new1, aes(x, f, col = grp), size = 0.8)
p.ri.th <- p +
  theme_void() +
  theme(legend.position = "none"); p.ri.th
p.ri <- p +
  theme_reg +
  theme(legend.position = "none") +
  labs(title = "Fitted regression line from the canonical + Pearson RKHS"); p.ri

# Random slope
mod2 <- iprior(y ~ . ^ 2, dat)
f <- predict(mod2, dat.new)
dat.new2 <- cbind(dat.new, f)
p <- ggplot(dat) +
  geom_point(aes(x, y, col = grp)) +
  geom_line(data = dat.new2, aes(x, f, col = grp), size = 0.8)
p.rs.th <- p +
  theme_void() +
  theme(legend.position = "none"); p.rs.th
p.rs <- p +
  theme_reg +
  theme(legend.position = "none") +
  labs(title = "Fitted regression line from the canonical x Pearson RKHS"); p.rs

ggsave("image/reg_multi_int_th.png", p.ri.th, "png", width = 4, height = 8 / 3, dpi = 150, scale = 1.1)
ggsave("image/reg_multi_int.png", p.ri, "png", width = 9, height = 6)
ggsave("image/reg_multi_slope_th.png", p.rs.th, "png", width = 4, height = 8 / 3, dpi = 150, scale = 1.1)
ggsave("image/reg_multi_slope.png", p.rs, "png", width = 9, height = 6)
move_fig_to_site()
