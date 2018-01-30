source("site-00-prelim.R")

theme_kernel_path <- theme_classic() +
  theme(legend.position = "none", axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank())

theme_kernel_path_th <- theme_void() + theme(legend.position = "none")

kernel_path_canonical <- function(n = 1000, seed = 789, th = FALSE) {
  x <- seq(-1, 1, length = n)
  m <- 5  # no. of paths
  f.prior <- matrix(NA, nrow = n, ncol = m)
  if (!is.null(seed)) set.seed(seed)
  for (i in seq_len(m)) f.prior[, i] <- iprior::fnH2(x) %*% rnorm(n)
  plot.df <- data.frame(x, f = f.prior)
  plot.df <- melt(plot.df, id = "x")
  p <- ggplot(plot.df, aes(x, value)) +
    geom_line(aes(col = variable))
  if (isTRUE(th)) p <- p + theme_kernel_path_th
  else p <- p +
    labs(x = expression(italic(x)), y = expression(italic(f(x))),
         title = "Sample linear I-prior paths") +
    theme_kernel_path
  p
}

kernel_path_fbm <- function(n = 1000, seed = 789, hurst = 0.1, th = FALSE) {
  the.title <- paste0("Sample fBm I-prior paths (Hurst = ", hurst, ")")
  x <- seq(-1, 1, length = n)
  m <- 5  # no. of paths
  f.prior <- matrix(NA, nrow = n, ncol = m)
  if (!is.null(seed)) set.seed(seed)
  for (i in seq_len(m)) f.prior[, i] <- iprior::fnH3(x, gamma = hurst) %*% rnorm(n)
  plot.df <- data.frame(x, f = f.prior)
  plot.df <- melt(plot.df, id = "x")
  p <- ggplot(plot.df, aes(x, value)) +
    geom_line(aes(col = variable))
  if (isTRUE(th)) p <- p + theme_kernel_path_th
  else p <- p +
    labs(x = expression(italic(x)), y = expression(italic(f(x))),
         title = the.title) +
    theme_kernel_path
  p
}

kernel_path_pearson <- function(n = 1000, seed = 789, th = FALSE) {
  x <- factor(sample(1:5, size = n, replace = TRUE))
  m <- 5  # no. of paths
  f.prior <- matrix(NA, nrow = n, ncol = m)
  if (!is.null(seed)) set.seed(seed)
  for (i in seq_len(m)) f.prior[, i] <- iprior::fnH1(x) %*% rnorm(n)
  plot.df <- data.frame(x, f = f.prior)
  plot.df <- plot.df[match(unique(plot.df$x), plot.df$x), ]
  plot.df <- melt(plot.df, id = "x")
  p <- ggplot(plot.df, aes(x, value)) +
    geom_point(aes(col = variable), size = 3, shape = 21, fill = NA, stroke = 1.2) +
    geom_point(aes(fill = variable, col = variable), size = 3, shape = 21, alpha = 0.5)
  if (isTRUE(th)) p <- p + theme_kernel_path_th
  else p <- p +
    labs(x = expression(italic(x)), y = expression(italic(f(x))),
         title = "Sample Pearson I-prior points") +
    theme_kernel_path
  p
}

kernel_path_constant <- function(n = 1000, seed = 789, th = FALSE) {
  x <- seq(-1, 1, length = n)
  m <- 5  # no. of paths
  f.prior <- matrix(NA, nrow = n, ncol = m)
  if (!is.null(seed)) set.seed(seed)
  for (i in seq_len(m)) f.prior[, i] <- matrix(1, nrow = n, ncol = n) %*% rnorm(n)
  plot.df <- data.frame(x, f = f.prior)
  plot.df <- melt(plot.df, id = "x")
  p <- ggplot(plot.df, aes(x, value)) +
    geom_line(aes(col = variable)) +
    coord_cartesian(ylim = c(min(f.prior) - 10, max(f.prior) + 10))
  if (isTRUE(th)) p <- p + theme_kernel_path_th
  else p <- p +
    labs(x = expression(italic(x)), y = expression(italic(f(x))),
         title = "Sample constant I-prior paths") +
    theme_kernel_path
  p
}

kernel_path_poly <- function(n = 1000, seed = 789, c = 0, d = 2, th = FALSE) {
  the.title <- paste0("Sample polynomial I-prior paths (degree = ", d, ")")
  x <- seq(-1, 1, length = n)
  m <- 5  # no. of paths
  f.prior <- matrix(NA, nrow = n, ncol = m)
  if (!is.null(seed)) set.seed(seed)
  for (i in seq_len(m)) f.prior[, i] <- fnH5(x, c = c, d = d) %*% rnorm(n)
  plot.df <- data.frame(x, f = f.prior)
  plot.df <- melt(plot.df, id = "x")
  p <- ggplot(plot.df, aes(x, value)) +
    geom_line(aes(col = variable))
  if (isTRUE(th)) p <- p + theme_kernel_path_th
  else p <- p +
    labs(x = expression(italic(x)), y = expression(italic(f(x))),
         title = the.title) +
    theme_kernel_path
  p
}


kernel_path_se <- function(n = 1000, seed = 789, l = 0.5, th = FALSE) {
  x <- seq(-1, 1, length = n)
  m <- 5  # no. of paths
  f.prior <- matrix(NA, nrow = n, ncol = m)
  if (!is.null(seed)) set.seed(seed)
  for (i in seq_len(m)) f.prior[, i] <- fnH4(x, l = l) %*% rnorm(n)
  plot.df <- data.frame(x, f = f.prior)
  plot.df <- melt(plot.df, id = "x")
  p <- ggplot(plot.df, aes(x, value)) +
    geom_line(aes(col = variable))
  if (isTRUE(th)) p <- p + theme_kernel_path_th
  else p <- p +
    labs(x = expression(italic(x)), y = expression(italic(f(x))),
         title = bquote(Sample~Gaussian~"I-prior"~paths~(sigma~"="~.(l)))) +
    theme_kernel_path
  p
}

## ---- kernel.path.plots ----
move_fig_to_site <- function() {
  files <- list.files("image/")
  file.copy(file.path("image", files),
            file.path("/Users/haziqjamil/Desktop/phd-poster-site/images", files),
            overwrite = TRUE)
}

# Canonical
ggsave("image/kernel_path_canonical_th.png", kernel_path_canonical(th = TRUE),
       "png", width = 3, height = 2)
ggsave("image/kernel_path_canonical.png", kernel_path_canonical(),
       "png", width = 9, height = 6)
move_fig_to_site()

# fBm
ggsave("image/kernel_path_fbm_01_th.png", kernel_path_fbm(hurst = 0.1, th = TRUE),
       "png", width = 3, height = 2)
ggsave("image/kernel_path_fbm_03_th.png", kernel_path_fbm(hurst = 0.3, th = TRUE),
       "png", width = 1.5, height = 1)
ggsave("image/kernel_path_fbm_05_th.png", kernel_path_fbm(hurst = 0.5, th = TRUE),
       "png", width = 1.5, height = 1)
ggsave("image/kernel_path_fbm_07_th.png", kernel_path_fbm(hurst = 0.7, th = TRUE),
       "png", width = 1.5, height = 1)
ggsave("image/kernel_path_fbm_09_th.png", kernel_path_fbm(hurst = 0.9, th = TRUE),
       "png", width = 1.5, height = 1)
ggsave("image/kernel_path_fbm_01.png", kernel_path_fbm(hurst = 0.1), "png", width = 9, height = 6)
ggsave("image/kernel_path_fbm_03.png", kernel_path_fbm(hurst = 0.3), "png", width = 9, height = 6)
ggsave("image/kernel_path_fbm_05.png", kernel_path_fbm(hurst = 0.5), "png", width = 9, height = 6)
ggsave("image/kernel_path_fbm_07.png", kernel_path_fbm(hurst = 0.7), "png", width = 9, height = 6)
ggsave("image/kernel_path_fbm_09.png", kernel_path_fbm(hurst = 0.9), "png", width = 9, height = 6)
move_fig_to_site()

# Pearson
ggsave("image/kernel_path_pearson_th.png", kernel_path_pearson(th = TRUE),
       "png", width = 6, height = 4)
ggsave("image/kernel_path_pearson.png", kernel_path_pearson(), "png",
       width = 9, height = 6)
move_fig_to_site()

# Constant
ggsave("image/kernel_path_const_th.png", kernel_path_constant(th = TRUE), "png",
       width = 1.5, height = 1)
ggsave("image/kernel_path_const.png", kernel_path_constant(), "png",
       width = 9, height = 6)
move_fig_to_site()

# Polynomial
ggsave("image/kernel_path_poly_2_th.png", kernel_path_poly(d = 2, th = TRUE),
       "png", width = 1.5, height = 1)
ggsave("image/kernel_path_poly_3_th.png", kernel_path_poly(d = 3, th = TRUE),
       "png", width = 1.5, height = 1)
ggsave("image/kernel_path_poly_4_th.png", kernel_path_poly(d = 4, th = TRUE),
       "png", width = 1.5, height = 1)
ggsave("image/kernel_path_poly_2.png", kernel_path_poly(d = 2),
       "png", width = 9, height = 6)
ggsave("image/kernel_path_poly_3.png", kernel_path_poly(d = 3),
       "png", width = 9, height = 6)
ggsave("image/kernel_path_poly_4.png", kernel_path_poly(d = 4),
       "png", width = 9, height = 6)
move_fig_to_site()

# SE
ggsave("image/kernel_path_se_01_th.png", kernel_path_se(l = 0.1, th = TRUE),
       "png", width = 1.5, height = 1)
ggsave("image/kernel_path_se_05_th.png", kernel_path_se(l = 0.5, th = TRUE),
       "png", width = 1.5, height = 1)
ggsave("image/kernel_path_se_1_th.png", kernel_path_se(l = 1, th = TRUE),
       "png", width = 1.5, height = 1)
ggsave("image/kernel_path_se_01.png", kernel_path_se(l = 0.1),
       "png", width = 9, height = 6)
ggsave("image/kernel_path_se_05.png", kernel_path_se(l = 0.5),
       "png", width = 9, height = 6)
ggsave("image/kernel_path_se_1.png", kernel_path_se(l = 1),
       "png", width = 9, height = 6)
move_fig_to_site()
