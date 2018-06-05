# Require iprior v0.6.4.9002 or higher for GPR support
library(iprior)
library(pushoverr)
library(RPEnsemble)
library(ggplot2)
library(foreach)
library(doSNOW)
no.cores <- parallel::detectCores()  # or set number of cores

# For push notifications using pushoverR, set here
# For more information, check out https://github.com/briandconnelly/pushoverr
userID <- "uyq2g37vnityt1b3yvpyicv6o9h456"
appToken <- "avxnrig1qppsgsw9woghwwmxsobo4a"

# Function to combine mean and se
mean_and_se <- function(x, y) paste0(dec_plac(x, 2), " (", dec_plac(y, 2), ")")

# Function to calculate ranks
tab_rank <- function(x, y) {
  # This is based on a weighted average. More weight given if low classification
  # error in small sammple size
  tmp <- apply(x + y, 1, function(x) sum(n * x) / sum(n))
  rank(tmp)
}

# Function to tabulate mean and se
tab_res <- function(...) {
  this <- list(...)
  K <- length(this)
  tab.mean <- tab.se <- tab <- NULL

  for (k in 1:K) {
    if (any(is.na(this[[k]]))) {
      tab.mean <- rbind(tab.mean, NA)
      tab.se <- rbind(tab.se, NA)
      tab <- rbind(tab, NA)
    } else {
      tab.mean.tmp <- apply(this[[k]], 2, mean)
      tab.se.tmp <- apply(this[[k]], 2, sd) / sqrt(nrow(this[[k]]))
      tab.mean.and.se <- mean_and_se(tab.mean.tmp, tab.se.tmp)
      tab.mean <- rbind(tab.mean, tab.mean.tmp)
      tab.se <- rbind(tab.se, tab.se.tmp)
      tab <- rbind(tab, tab.mean.and.se)
    }
  }

  rownames(tab.mean) <- rownames(tab.se) <- rownames(tab) <- names(this)
  colnames(tab.mean) <- colnames(tab.se) <- colnames(tab) <- colnames(this[[1]])
  list(
    tab = as.data.frame(tab),
    tab.mean = as.data.frame(tab.mean),
    tab.se = as.data.frame(tab.se)
  )
}

# I-prior probit simulation
my_iprobit_sim <- function(nsim = 100, kernel = c("canonical", "fbm", "se")) {
  # kernel <- match.arg(kernel, c("canonical", "fbm", "se"))
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  progress <- function(i) setTxtProgressBar(pb, i)

  cl <- makeCluster(no.cores)
  registerDoSNOW(cl)
  res <- foreach(i = 1:nsim, .combine = rbind,
                 .packages = c("iprior", "iprobit"),
                 .export = c("n", "y", "X"),
                 .options.snow = list(progress = progress)) %dopar% {
      res.tmp <- rep(NA, length(n))
      for (j in 1:length(n)) {
        tmp <- iprobit(y, X, kernel = kernel, control = list(maxit = 5),
                       train.samp = sample(seq_along(y), size = n[j]))
        res.tmp[j] <- predict(tmp)$error
      }
      res.tmp
    }
  close(pb)
  stopCluster(cl)
  # save.image(experiment.name)

  push.message <- paste0(
    experiment.name, ": ", kernel, " I-prior probit COMPLETED."
  )
  pushoverr::pushover(message = push.message, user = userID, app = appToken)

  colnames(res) <- paste0(c("n = "), n)
  res
}

# Function to plot
plot_res <- function() {
  plot.df <- cbind(tab.mean, id = rownames(tab.mean))
  # plot.df <- plot.df[names(sort(tab.ranks)), ]
  suppressMessages(
    plot.df <- reshape2::melt(plot.df)
  )
  id2 <- plot.df$id
  tmp <- levels(id2)
  tmp[grepl("I-probit", tmp)] <- "I-probit"
  levels(id2) <- tmp
  suppressMessages(
    plot.se <- reshape2::melt(tab.se)
  )
  plot.df <- cbind(plot.df, se = plot.se[, 2], id2 = id2)
  plot.df$id <- factor(plot.df$id,
                       levels = names(sort(tab.ranks, decreasing = TRUE)))
  ggplot(plot.df, aes(x = value, y = id, col = id2, label = dec_plac(value))) +
    geom_point() +
    # directlabels::geom_dl(method = "top.bumptwice") +
    geom_label(size = 3.8) +
    # geom_errorbarh(aes(xmin = value - 1.96 * se, xmax = value + 1.96 * se,
    #                    height = 0)) +
    facet_grid(. ~ variable) +
    coord_cartesian(xlim = c(18, 45)) +
    labs(x = "Misclassification rate", y = NULL) + guides(col = FALSE) +
    theme_bw()
}
