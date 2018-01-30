# Chapter 7, Section 4
# Vowel recognition
source("00-prelim.R")

## ---- vowel.data ----
data("vowel.train")
data("vowel.test")
vowel.train$y <- as.factor(vowel.train$y)
vowel.test$y <- as.factor(vowel.test$y)

## ---- vowel.lb ----
# tab.lb <- matrix(NA, nrow = 10, ncol = 10)
# for (i in 1:10) {
#   tmp <- iprobit(vowel.train$y, vowel.train[, -1], kernel = "FBM",
#                  control = list(maxit = 10))
#   tab.lb[, i] <- tmp$lower.bound
# }
# save(tab.lb, file = "data/vowel-tab-lb")
load("data/vowel-tab-lb")
as.tibble(tab.lb[2:10, ]) %>%
  mutate(x = seq_len(9)) %>%
  reshape2::melt(id = "x") %>%
  ggplot(aes(x, value, group = variable, col = variable)) +
  geom_line() +
  theme_bw() +
  guides(col = FALSE)

## ---- vowel.mod.can ----
# (mod.can <- iprobit(vowel.train$y, vowel.train[, -1],
#                     control = list(restarts = 8, maxit = 2000,
#                                    restart.method = "error")))
# iplot_lb(mod.can)
# iplot_error(mod.can)
# save(mod.can, file = "data/vowel-mod-can")
load("data/vowel-mod-can")
(res.can <- predict(mod.can, list(vowel.test[, -1]), vowel.test$y))

## ---- vowel.mod.fbm ----
# (mod <- iprobit(vowel.train$y, vowel.train[, -1], kernel = "FBM",
#                 control = list(restarts = 8, maxit = 500, restart.method = "error")))
# iplot_lb(mod)
# iplot_error(mod)
# save(mod, file = "data/vowel-mod-fbm")
load("data/vowel-mod-fbm")
res <- predict(mod, list(vowel.test[, -1]), vowel.test$y)

## ---- vowel.confusion.matrix ----
tibble(test = vowel.test$y, predicted = res.can$y) %>%
  group_by(test, predicted) %>%
  tally() %>%
  mutate(prop = n / 42) %>%
  complete(predicted, fill = list(n = 0, prop = 0)) -> plot.df

levels(plot.df$test) <- levels(plot.df$predicted) <- c(
  "heed", "hid", "head", "had", "hud", "hard", "hod", "hoard", "hood", "who'd", "heard"
)

ggplot(plot.df, aes(test, predicted)) +
  geom_tile(aes(fill = n)) +
  geom_text(data = subset(plot.df, n > 0), aes(label = n)) +
  scale_y_discrete(limits = rev(levels(plot.df$test))) +
  scale_x_discrete(position = "top") +
  scale_fill_continuous(low = "white", high = "slategray", limits = c(0, 42)) +
  labs(x = "Test data", y = "Predicted classes") +
  theme_bw() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0))
  # ggtitle("Canonical kernel")

tibble(test = vowel.test$y, predicted = res$y) %>%
  group_by(test, predicted) %>%
  tally() %>%
  mutate(prop = n / 42) %>%
  complete(predicted, fill = list(n = 0, prop = 0)) -> plot.df

levels(plot.df$test) <- levels(plot.df$predicted) <- c(
  "heed", "hid", "head", "had", "hud", "hard", "hod", "hoard", "hood", "who'd", "heard"
)

ggplot(plot.df, aes(test, predicted)) +
  geom_tile(aes(fill = n)) +
  geom_text(data = subset(plot.df, n > 0), aes(label = n)) +
  scale_y_discrete(limits = rev(levels(plot.df$test))) +
  scale_x_discrete(position = "top") +
  scale_fill_continuous(low = "white", high = "slategray", limits = c(0, 42)) +
  labs(x = "Test data", y = "Predicted classes (fBm-0.5 kernel)") +
  theme_bw() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
  # ggtitle("fBm-0.5 kernel")

## ---- vowel.tab ----
iprobit.fbm <- c(get_error_rate(mod), res$test.error)
iprobit.can <- c(get_error_rate(mod.can), res.can$test.error)
tab <- rbind(
  c("Linear regression"  ,              48, 67),
  c("Logistic regression",           22, 51),
  c("Linear discriminant analysis"    , 32, 56),
  c("Quadratic discriminant analysis" , 1, 53),
  c("Decision trees"                  , 5, 54),
  c("Neural networks"                 , "", 45),
  c("k-Nearest neighbours"            , "", 44),
  c("FDA/BRUTO"                       , 6, 44),
  c("FDA/MARS"                        , 13, 39),
  c("I-probit (fBm-0.5)"              , round(iprobit.fbm, 0)),
    c("I-probit (linear)"             , round(iprobit.can, 0))
)
colnames(tab) <- c("Method", "Training", "Test")
knitr::kable(tab, format = "latex", booktabs = TRUE, align = c("l", "r", "r"),
             linesep = c("", "", "", "\\addlinespace", "", "", "", "", "\\addlinespace"),
             caption = "Results of various classification methods for the vowel data set.") %>%
  kableExtra::kable_styling(position = "center") %>%
  kableExtra::add_header_above(c(" " = 1, "Error rates" = 2))
