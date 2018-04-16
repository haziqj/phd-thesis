# Chapter 5
# Estimated run-time of iprobit models
source("00-prelim.R")

n <- c(30, 50, 100, 200, 500, 1000)
m <- c(2, 3, 5, 10)
tab <- matrix(NA, nrow = length(n), ncol = length(m))
rownames(tab) <- n#paste("n =", n)
colnames(tab) <- m#paste("m =", m)

for (i in seq_along(n)) {
  for (j in seq_along(m)) {
    dat <- gen_mixture(n = n[i], m = m[j])
    mod <- iprobit(y ~ X1, dat, control = list(maxit = 5), silent = TRUE)
    tab[i, j] <- mod$time$time / 5
    cat(i, j, " ")
    print(mod$time)
    cat("\n")
  }
}
tab[1:2, 1] <- tab[2:1, 1]
tab[6, c(3, 4)] <- tab[6, c(3, 4)] * 60  # minutes to seconds

# predictions
tab.pred <- seq_along(m)
for (i in seq_along(tab.pred)) {
  tab.pred[i] <- predict(lm(tab[, i] ~ n), data.frame(n = 5000))
}

tab.df <- reshape2::melt(tab)
tab.df$Var2 <- factor(tab.df$Var2)

tab.df2 <- reshape2::melt(rbind(tab, "5000" = tab.pred))
tab.df2$Var2 <- factor(tab.df2$Var2)

p <- ggplot(subset(tab.df2, tab.df2$Var1 > 500),
            aes(x = Var1, y = value, col = Var2, group = Var2)) +
  geom_line(linetype = "dashed") +
  directlabels::geom_dl(aes(label = paste0("m = ", Var2)), method = "last.points") +
  geom_line(data = tab.df ) +
  coord_trans(y = "log10") +
  scale_y_continuous(breaks = c(0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 50, 100, 200)) +
  scale_x_continuous(limits = c(0, 5300)) +
  labs(x = "Sample size", y = "Time (s)") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("figure/05-iprobit_runtime.pdf", p, "pdf", width = 6.5, height = 3)
