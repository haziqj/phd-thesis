# Chapter 4
# Estimated run-time and memory usage of iprior models
library(iprior)
library(tidyverse)
chapter.no <- "04"

# Memory usage ANOVA kernels (no need to fit)
n <- c(30, 50, 100, 200, 500, 1000)
p <- c(1, 2, 3, 4, 5)
tab <- matrix(NA, nrow = length(n), ncol = length(p))
rownames(tab) <- n
colnames(tab) <- p

for (i in seq_along(n)) {
  dat <- gen_smooth(n = n[i])
  X <- dat$X
  for (j in seq_along(p)) {
    if (j > 1) dat <- cbind(dat, X)
    stopifnot(p[j] == (ncol(dat) - 1))
    colnames(dat) <- c("y", seq_len(p[j]))
    mod <- kernL(y ~ . ^ 2, dat, kernel = "fbm")
    tab[i, j] <- as.numeric(get_size(mod)) / (1000 ^ 2)  # MB
    cat(i, j, " ")
    cat("\n")
  }
}

# prediction (n = 5000)
tab.pred <- seq_along(p)
for (i in seq_along(tab.pred)) {
  tab.pred[i] <- predict(lm(tab[, i] ~ n), data.frame(n = 5000))
}

tab.df <- reshape2::melt(tab)
tab.df$Var2 <- factor(tab.df$Var2)

tab.df2 <- reshape2::melt(rbind(tab, "5000" = tab.pred))
tab.df2$Var2 <- factor(tab.df2$Var2)

p1 <- ggplot(subset(tab.df2, tab.df2$Var1 > 500),
            aes(x = Var1, y = value, col = Var2, group = Var2)) +
  geom_line(linetype = "dashed") +
  directlabels::geom_dl(aes(label = paste0("p = ", Var2)), method = "last.points") +
  geom_line(data = tab.df ) +
  coord_trans(y = "log10") +
  scale_y_continuous(breaks = c(1, 10, 100, 1000)) +
  scale_x_continuous(limits = c(0, 5300)) +
  labs(x = "Sample size", y = "Object size (MB)") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("figure/04-iprior_size.pdf", p1, "pdf", width = 6.5, height = 3.5)

# I-prior runtime (no need to fit)
n <- c(30, 50, 100, 200, 500, 1000)
m <- 1:5
tab <- matrix(NA, nrow = length(n), ncol = length(m))
rownames(tab) <- n#paste("n =", n)
colnames(tab) <- m#paste("m =", m)
tab.em <- tab.di <- tab
tab <- tab[, 1:2]

for (i in seq_along(n)) {
  dat <- gen_smooth(n = n[i])
  for (j in seq_along(m)) {
    cat(i, j, " ")
    mod <- iprior(y ~ X, dat, control = list(maxit = 1000, silent = TRUE))
    tab.di[i, j] <- mod$time$time
    print(mod$time)
    mod <- iprior(y ~ X, dat, method = "em",
                  control = list(maxit = 1000, silent = TRUE))
    tab.em[i, j] <- mod$time$time
    print(mod$time)
    cat("\n")
  }
}
tab.em[5, 3] <- tab.em[5, 3] * 60  # minutes to seconds
tab.em[6, ] <- tab.em[6, ] * 60  # minutes to seconds
tab.di[6, ] <- tab.di[6, ] * 60  # minutes to seconds
tab[, 1] <- apply(tab.em, 1, mean, na.rm = TRUE)
tab[, 2] <- apply(tab.di, 1, mean, na.rm = TRUE)

# predictions
tab.pred <- seq_len(2)
for (i in seq_along(tab.pred)) {
  tab.pred[i] <- predict(lm(tab[, i] ~ n), data.frame(n = 5000))
}

tab.df <- reshape2::melt(tab)
tab.df$Var2 <- factor(tab.df$Var2)

tab.df2 <- reshape2::melt(rbind(tab, "5000" = tab.pred))
tab.df2$Var2 <- factor(tab.df2$Var2)

p2 <- ggplot(subset(tab.df2, tab.df2$Var1 > 500),
            aes(x = Var1, y = value, col = Var2, group = Var2)) +
  geom_line(linetype = "dashed") +
  directlabels::geom_dl(aes(label = ifelse(Var2 == 1, "EM", "Direct")),
                        method = "last.points") +
  geom_line(data = tab.df ) +
  coord_trans(y = "log10") +
  scale_y_continuous(breaks = c(1, 10, 100, 1000)) +
  scale_x_continuous(limits = c(0, 5300)) +
  labs(x = "Sample size", y = "Time (s)") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("figure/04-iprior_runtime.pdf", p2, "pdf", width = 6.5, height = 3)

move_fig_to_chapter()
