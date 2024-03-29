# Chapter 7, Section 3
# Meta-analysis of nicotine gum treatment for smoking cessation
source("00-prelim.R")

## ---- data.smoke ----
smoking <- read.table(
  "data/gum.txt",
  sep = "\t", header = TRUE
)
sum(smoking[, c(3, 5)])  # sample size
smoking.small <- smoking

levels(smoking[, 1]) <- c(levels(smoking[, 1]), "Summary measure")
as.tibble(
  rbind(smoking,
        c("Summary measure", apply(smoking[, -1], 2, sum)))
) %>%
  mutate(d1 = as.numeric(d1),  # number of patients quit in treatment group
         n1 = as.numeric(n1),  # number of patients in treatment group
         d0 = as.numeric(d0),  # number of patients quit in control group
         n0 = as.numeric(n0),  # number of patients in control group
         logodds = log((d1 / (n1 - d1)) / (d0 / (n0 - d0))),
         lower = NA,
         upper = NA,
         n = n0 + n1,
         method = "Raw odds") %>%
  select(Study, lower, logodds, upper, n, method) -> smoke.raw.odds
study.size <- smoke.raw.odds$n

makeSmoke <- function(x) {
  rbind(
    data.frame(y = rep(1, x[2]), Study = x[1], Group = "Treated"),  # Quit
    data.frame(y = rep(0, x[3] - x[2]), Study = x[1], Group = "Treated"),  # Remain
    data.frame(y = rep(1, x[4]), Study = x[1], Group = "Control"),  # Quit
    data.frame(y = rep(0, x[5] - x[4]), Study = x[1], Group = "Control")  # Remain
  )
}
dat.smoke <- NULL
for (i in 1:nrow(smoking.small)) {
  dat.smoke <- rbind(dat.smoke, makeSmoke(smoking.small[i, ]))
}  # ignore warning about row names discarded
dat.smoke$y <- as.factor(dat.smoke$y)
levels(dat.smoke$y) <- c("Remain", "Quit")
str(dat.smoke)

# Distribution among control and treatment group
tmp <- apply(smoking[, -1], 2, mean)[c(2, 4)]  # treatment vs control
(dist.mean <- tmp / sum(tmp) * 100)

# Raw odds ratio
(raw.odds.ratio <- exp(smoke.raw.odds$logodds[28]))

# Load Logistic RE model results
as.tibble(read.csv("data/smoking_re_results.csv")) %>%
  mutate(n = study.size, method = "Logistic RE",
         studynam = factor(studynam, levels = levels(smoke.raw.odds[[1]]))) %>%
  rename(Study = studynam, logodds = eff) -> smoke.re.res

## ---- plot.data.smoke ----
as.tibble(dat.smoke) %>%
  group_by(y, Group, Study) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(Group = factor(Group, levels = c("Control", "Treated")),
         Group = recode(Group, Treated = "Treatment")) %>%
  ggplot(aes(y, count, fill = y)) +
  facet_grid(. ~ Group) +
  geom_boxplot(notch = TRUE, outlier.alpha = 0) +
  geom_jitter(width = 0.12, height = 0, pch = 21) +
  # coord_flip() +
  # coord_trans(y = "log10") +
  scale_y_log10(breaks = c(1, 2, 5, 10, 25, 50, 100, 500)) +
  labs(x = NULL, y = "Count") +
  # coord_cartesian(ylim = c(0, 150)) +
  theme_bw() +
  theme(legend.position = "none")

## ---- mod.smoke ----
# This takes quite a long time, but saved object mod in data/smokingmod
# mod1 <- iprobit(y ~ Group, dat.smoke)
# save(mod1, file = "data/smoking-mod1")
# > mod1
# Training error rate: 23.65 %
# Lower bound value: -3210.76
#
# Class = 1 Class = 2
# Intercept  -0.72533  -0.72533
# res         1.00460   1.00460
# >
# mod2 <- iprobit(y ~ Group + Study, dat.smoke,
#                 control = list(theta0 = c(3.86e+03, 1.06e-03), alpha0 = -0.769))
# save(mod2, file = "data/smoking-mod2")
# > mod2
# Training error rate: 29.3 %
# Lower bound value: -3142.237
#
# Class = 1   Class = 2
# Intercept     -0.76267    -0.76267
# lambda[1,] 30593.69995 30593.69995
# lambda[2,]     0.00084     0.00084
# mod3 <- iprobit(y ~ Group * Study, dat.smoke,
#                 control = list(theta0 = c(0.014894, 0.000878), alpha0 = -0.767))
# save(mod3, file = "data/smoking-mod3")
# > mod3
# Training error rate: 23.48 %
# Lower bound value: -3091.204
#
# Class = 1 Class = 2
# Intercept   -0.76669  -0.76669
# lambda[1,]   0.01437   0.01437
# lambda[2,]   0.00087   0.00087

# ## ---- fit.smoke ----
# load("data/smoking-mod1"); l1 <- logLik(mod1); e1 <- get_error_rate(mod1)
# load("data/smoking-mod2"); l2 <- logLik(mod2); e2 <- get_error_rate(mod2)
# load("data/smoking-mod3"); l3 <- logLik(mod3); e3 <- get_error_rate(mod3)
# calc_odds <- function(n.samp = 100) {
#   studies <- levels(dat.smoke$Study)
#   as.tibble(dat.smoke) %>%
#     mutate(tmp = paste0(Study, Group)) -> tmp
#   unq.i <- match(unique(tmp$tmp), tmp$tmp)
#   quants <- predict(mod3, dat.smoke[unq.i, ], quantiles = TRUE, raw = TRUE,
#                     n.samp = n.samp)[[1]][[2]]  # prob quit
#
#   res <- quants[seq_len(length(studies) + 1), ]; res[] <- NA
#   for (i in seq_along(studies)) {
#     prob.treated <- quants[2 * i - 1, ]
#     prob.control <- quants[2 * i, ]
#     a <- prob.treated; b <- 1 - a  # probs in treatment groups
#     c <- prob.control; d <- 1 - c  # probs in control groups
#     # OR.quit = (a / b) / (c / d) = (a * d) / (c * b)
#     res[i, ] <- log(a * d) - log(b * c)
#   }
#
#   unq.i <- match(unique(dat.smoke$Group), dat.smoke$Group)
#   quants <- predict(mod1, dat.smoke[unq.i, ], quantiles = TRUE, raw = TRUE,
#                     n.samp = n.samp)[[1]][[2]]
#   prob.treated <- quants[1, ]
#   prob.control <- quants[2, ]
#   a <- prob.treated; b <- 1 - a
#   c <- prob.control; d <- 1 - c
#   res[nrow(res), ] <- log(a * d) - log(b * c)
#
#   tab <- data.frame(
#     Study = c(studies, "Summary measure"),
#     t(apply(res, 1, function(x) c(logodds = mean(x),
#                                   quantile(x, probs = c(0.05, 0.95))))),
#     n = study.size,
#     method = "I-probit"
#   )
#   names(tab)[3:4] <- c("lower", "upper")
#   tab
# }
# l <- c(l1, l2, l3)
# err <- c(e1, e2, e3)
# b <- c(get_brier_score(mod1), get_brier_score(mod2), get_brier_score(mod3))
# smoke.ip.res <- calc_odds(100)
# save(smoke.ip.res, file = "data/smoking-res")
# save(l, file = "data/smoking-lb")
# save(b, file = "data/smoking-brier")
# save(err, file = "data/smoking-error")

## ---- mod.compare.smoke ----
load("data/smoking-res")
load("data/smoking-lb")
load("data/smoking-brier")
load("data/smoking-error")
tab.compare <- as.data.frame(cbind(
  c("$f_1$", "$f_1 + f_2$",
    "$f_1 + f_2 + f_{12}$"),
  paste0("$-", dec_plac(-l, 2), "$"),

  dec_plac(b, 3),
  c(1, 2, 2)
))
colnames(tab.compare) <- c("Model", "Lower bound", "Brier score",
                           "No. of RKHS\\newline scale param.")
# kable(tab.compare, longtable = TRUE, booktabs = TRUE, escape = FALSE,
#       align = c("l", "r", "r", "r"), format = "latex",
#       caption = "Results of the I-prior model fit for three models.") %>%
#   column_spec(4, width = "2.4cm")

## ---- smoke.forest.plot ----
ggplot(rbind(smoke.raw.odds, smoke.re.res, smoke.ip.res), aes(Study, logodds)) +
  geom_hline(yintercept = 0, linetype = "dotted", col = "grey50") +
  geom_errorbar(aes(ymin = lower, ymax = upper, linetype = method, col = Study),
                width = 0.4, size = 0.5, position = position_dodge(width = 0.55)) +
  geom_point(aes(fill = Study, size = n, shape = method),
             position = position_dodge(width = 0.55), col = "white", alpha = 0.75) +
  scale_x_discrete(name = NULL, limits = rev(levels(smoke.re.res$Study))) +
  scale_y_continuous(name = NULL,
                     breaks = c(-0.5, 0, 0.5, 1, 1.5)
                     ) +
  scale_size_continuous(range = c(1.8, 8)) +
  coord_flip() +
  theme_bw() +
  scale_linetype_manual(name = NULL, values = c(1, 2, 0)) +
  scale_shape_manual(name = NULL, values = c(21, 22, 23)) +
  guides(col = FALSE, size = FALSE, fill = FALSE,
         linetype = guide_legend(override.aes = list(col = "black",
                                                     fill = "black",
                                                     # size = 2,
                                                     linetype = c(1, 2, 0)))) +
  # theme(legend.key.width = unit(2.5, "line"),
  #        legend.justification = c(1, 0), legend.position = c(0.99, 0),
  #       legend.background = element_rect(fill = scales::alpha("white", 0)))
  theme(legend.position = "top", legend.key.width = unit(2.5, "line"))
