# Chapter 4, Section 3
# Longitudinal modelling using I-priors (cows data set)
source("00-prelim.R")

## ---- data ----
data(cattle, package = "jmcm")
names(cattle) <- c("id", "time", "group", "weight")
cattle$id <- as.factor(cattle$id)  # convert to factors
# str(cattle)
as.tibble(cattle)

## ---- ipriorKernels ----
# Model 1: weight ~ f(time)
mod1 <- kernL(weight ~ time, data = cattle, model = list(kernel = "FBM"))

# Model 2: weight ~ f(time) + f(treatment) + f(time dependent treatment)
mod2 <- kernL(weight ~ group * time, cattle,
              model = list(kernel = "FBM"))

# Model 3: weight ~ f(time) + f(cow index) + f(time dependent cow index)
mod3 <- kernL(weight ~ id * time, cattle,
              model = list(kernel = "FBM"))

# Model 4: weight ~ f(time) + f(cow index) +  f(treatment)
#                   + f(time dependent cow index)
#                   + f(time dependent treatment)
mod4 <- kernL(weight ~ group * time +  id * time, cattle,
              model = list(kernel = "FBM"))

# Model 5: weight ~ f(time:cow:treatment)
mod5 <- kernL(weight ~ id * group * time, cattle,
              model = list(kernel = "FBM"))

## ---- model_fitting ----
mod1 <- kernL(weight ~ time, data = cattle, model = list(kernel = "FBM"))
mod1.fit <- ipriorOptim(mod1)
mod2.fit <- ipriorOptim(mod2, control = list(silent = TRUE))
mod3.fit <- ipriorOptim(mod3, control = list(silent = TRUE))
mod4.fit <- ipriorOptim(mod4, control = list(silent = TRUE))
mod5.fit <- ipriorOptim(mod5, control = list(silent = TRUE))

## ---- table ----
cow_table <- function(mod) {
  form <- capture.output(mod$ipriorKernel$formula)
  form <- substring(form, 10)
  tibble(formula = form, loglik = logLik(mod), error = mod$sigma,
         no_lambda = length(mod$lambda))
}
tab <- rbind(
  cow_table(mod1.fit), cow_table(mod2.fit), cow_table(mod3.fit),
  cow_table(mod4.fit), cow_table(mod5.fit)
)
tab <- cbind(model = 1:5, tab)
knitr::kable(tab, col.names = c(
  "Model", "Formula", "Log-likelihood", "Error S.D.", "No. of $\\lambda$"
))

## ---- cow_plot ----
# Create data for plot (new points)
# group_indices <- match(levels(cattle$id), cattle$id)
# cow_groups <- as.character(cattle$group[group_indices])
# nx <- 1000
# ncows <- length(cow_groups)
# new.time <- seq(0, 133, length = nx)
# plot.df <- data.frame(id = 1, time = new.time,
#                       group = cow_groups[1], weight = 1)
# for (i in seq_len(ncows)[-1]) {
#   plot.df <- rbind(plot.df,
#                    data.frame(id = i, time = new.time,
#                               group = cow_groups[i], weight = 1))
# }
# plot.df$id <- factor(plot.df$id)
# plot.df$weight <- predict(mod5.fit, plot.df)
# as.tibble(plot.df) %>%
#   mutate(group  = recode_factor(group, A = "Treatment A", B = "Treatment B")) -> plot.df

as.tibble(cattle) %>%
  select(id, time, group, weight) %>%
  mutate(iprior = fitted(mod5.fit),
         group  = recode_factor(group, A = "Treatment A", B = "Treatment B")) %>%
  ggplot(aes(x = time, y = weight, col = id)) +
  # geom_point(size = 0.5) +
  # geom_line(data = plot.df, aes(x = time, y = weight, col = id)) +
  geom_line(aes(y = iprior)) +
  # geom_dl(aes(label = id), method = list("last.bumpup", cex = 0.8)) +
  facet_grid(. ~ group) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = unique(cattle$time)[-11]) +
  labs(x = "Time (days)", y = "Weight (kg)", title = "Fitted growth curves")
