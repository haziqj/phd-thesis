# Chapter 5, examples
# Vowel recognition
library(iprobit)
library(tidyverse)
library(reshape2)
library(kernlab)

## ---- vowel.data ----
data("vowel.train", package = "ElemStatLearn")
data("vowel.test", package = "ElemStatLearn")
vowel.train$y <- as.factor(vowel.train$y); n.train <- nrow(vowel.train)
vowel.test$y <- as.factor(vowel.test$y); n.test <- nrow(vowel.test)
vowel.dat <- rbind(vowel.train, vowel.test)

## ---- vowel.lb ----
# Plot showing different starting values leads to multiple local optima. Best to
# restart many times and pick one with best lb/error rate
#
# tab.lb <- matrix(NA, nrow = 10, ncol = 10)
# for (i in 1:10) {
#   tmp <- iprobit(vowel.train$y, vowel.train[, -1], kernel = "fbm",
#                  control = list(maxit = 10))
#   tab.lb[, i] <- tmp$lower.bound
# }
# save(tab.lb, file = "data/vowel-tab-lb")
# load("data/vowel-tab-lb")
# as.tibble(tab.lb[2:10, ]) %>%
#   mutate(x = seq_len(9)) %>%
#   reshape2::melt(id = "x") %>%
#   ggplot(aes(x, value, group = variable, col = variable)) +
#   geom_line() +
#   theme_bw() +
#   guides(col = FALSE)

# For each of the models, want to see if different starting values yields
# different lbs. But it is ok, so no need to rerun for final model.

## ---- vowel.mod.can ----
# mod.can <- iprobit(y ~ ., vowel.dat, train.samp = 1:n.train, one.lam = TRUE,
#                    control = list(maxit = 10, restarts = TRUE))
# iplot_par_lb(mod.can)
# iplot_par_error(mod.can)
# iplot_par_error(mod.can, type = "test")
# mod.can <- iprobit(y ~ ., vowel.dat, train.samp = 1:n.train, one.lam = TRUE,
#                    control = list(maxit = 1000))
# save(mod.can, file = "data/vowel-mod-can")
load("data/vowel-mod-can")

## ---- vowel.mod.fbm ----
# mod.fbm <- iprobit(y ~ ., vowel.dat, train.samp = 1:n.train, one.lam = TRUE,
#                    kernel = "fbm", control = list(maxit = 10, restarts = TRUE))
# iplot_par_lb(mod.fbm)
# iplot_par_error(mod.fbm)
# iplot_par_error(mod.fbm, type = "test")
# mod.fbm <- iprobit(y ~ ., vowel.dat, train.samp = 1:n.train, one.lam = TRUE,
#                    kernel = "fbm", control = list(maxit = 1000))
# save(mod.fbm, file = "data/vowel-mod-fbm")
load("data/vowel-mod-fbm")

## ---- vowel.mod.se ----
# mod.se <- iprobit(y ~ ., vowel.dat, train.samp = 1:n.train, one.lam = TRUE,
#                   kernel = "se", control = list(maxit = 10, restarts = TRUE))
# iplot_par_lb(mod.se)
# iplot_par_error(mod.se)
# iplot_par_error(mod.se, type = "test")
# mod.se <- iprobit(y ~ ., vowel.dat, train.samp = 1:n.train, one.lam = TRUE,
#                   kernel = "se", control = list(maxit = 1000))
# save(mod.se, file = "data/vowel-mod-se")
load("data/vowel-mod-se")

## ---- vowel.confusion.matrix ----
conf_mat <- function(mod) {
  tibble(test = vowel.test$y, predicted = mod$test$y) %>%
    group_by(test, predicted) %>%
    tally() %>%
    mutate(prop = n / 42) %>%
    complete(predicted, fill = list(n = 0, prop = 0)) -> plot.df

  levels(plot.df$test) <- levels(plot.df$predicted) <- c(
    "heed", "hid", "head", "had", "hud", "hard", "hod", "hoard", "hood", "who'd", "heard"
  )

  ggplot(plot.df, aes(test, predicted)) +
    geom_tile(aes(fill = n)) +
    geom_text(aes(label = ifelse(n > 0, n, NA)), na.rm = TRUE) +
    scale_y_discrete(limits = rev(levels(plot.df$test))) +
    scale_x_discrete(position = "top") +
    scale_fill_continuous(low = "white", high = "slategray", limits = c(0, 42)) +
    labs(x = "Test data", y = "Predicted classes") +
    theme_bw() +
    theme(legend.position = "none", plot.title = element_text(hjust = 0))
  # ggtitle("Canonical kernel")
}

conf_mat(mod.can)
conf_mat(mod.fbm)
conf_mat(mod.se)

## ---- vowel.tab ----
# iprobit.fbm <- c(get_error_rate(mod), res$test.error)
# iprobit.can <- c(get_error_rate(mod.can), res.can$test.error)
# tab <- rbind(
#   c("Linear regression"  ,              48, 67),
#   c("Logistic regression",           22, 51),
#   c("Linear discriminant analysis"    , 32, 56),
#   c("Quadratic discriminant analysis" , 1, 53),
#   c("Decision trees"                  , 5, 54),
#   c("Neural networks"                 , "", 45),
#   c("k-Nearest neighbours"            , "", 44),
#   c("FDA/BRUTO"                       , 6, 44),
#   c("FDA/MARS"                        , 13, 39),
#   c("I-probit (fBm-0.5)"              , round(iprobit.fbm, 0)),
#     c("I-probit (linear)"             , round(iprobit.can, 0))
# )
# colnames(tab) <- c("Method", "Training", "Test")
# knitr::kable(tab, format = "latex", booktabs = TRUE, align = c("l", "r", "r"),
#              linesep = c("", "", "", "\\addlinespace", "", "", "", "", "\\addlinespace"),
#              caption = "Results of various classification methods for the vowel data set.") %>%
#   kableExtra::kable_styling(position = "center") %>%
#   kableExtra::add_header_above(c(" " = 1, "Error rates" = 2))

## ---- vowel.gpc ----
mod.gpc1 <- gausspr(y ~ ., vowel.train, kernel = "vanilladot")
mod.gpc2 <- gausspr(y ~ ., vowel.train, kernel = "rbfdot")
y.hat <- predict(mod.gpc2, vowel.test)
RMSE.train <- sum(fitted(mod.gpc2) != vowel.train$y) / length(vowel.train$y) * 100
RMSE.test <- sum(y.hat != vowel.test$y) / length(vowel.test$y) * 100
