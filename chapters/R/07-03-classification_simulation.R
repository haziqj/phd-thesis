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

# Function to specify decimal places
decPlac <- function(x, k = 2) format(round(x, k), nsmall = k)

# Function to combine mean and se
meanAndSE <- function(x, y) paste0(decPlac(x, 2), " (", decPlac(y, 2), ")")

# Function to calculate ranks
tabRank <- function(x, y) {
  # This is based on a weighted average. More weight given if low classification
  # error in small sammple size
  tmp <- apply(x + y, 1, function(x) sum(n * x) / sum(n))
  rank(tmp)
}

# Function to create test and train set
testTrain <- function(n.testTrain, y.testTrain, X.testTrain, sex.testTrain,
                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  train.index <- sample(1:length(y.testTrain), size = n.testTrain,
                        replace = FALSE)
  test.index <- (1:length(y.testTrain))[-train.index]
  if (length(test.index) > 1000) {
    test.index <- sample(test.index, size = 1000, replace = FALSE)
  }

  X.train <- X.testTrain[train.index, ]
  sex.train <- sex.testTrain[train.index]
  y.train <- y.testTrain[train.index]
  X.test <- X.testTrain[test.index, ]
  sex.test <- sex.testTrain[test.index]
  y.test <- y.testTrain[test.index]

  list(X.train = X.train, sex.train = sex.train, y.train = y.train,
       X.test = X.test, sex.test = sex.test, y.test = y.test,
       N = length(y.testTrain))
}

# Function to tabulate mean and se
tabRes <- function(...) {
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
      tab.mean.and.se <- meanAndSE(tab.mean.tmp, tab.se.tmp)
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

# I-prior probit innersim
probitInnerSim <- function(y.innerSim, X.innerSim, sex.innerSim, n.innerSim,
                           kernel) {
  dat <- testTrain(n.innerSim, y.innerSim, X.innerSim, sex.innerSim)
  # Not letting the algorithm converge fully. Want a local optima for better
  # predictive performance.
  mod <- iprobit(dat$y.train, dat$X.train,
                 # dat$sex.train,
                 kernel = kernel, control = list(maxit = 5))
  y.test <- predict(mod, newdata = list(dat$X.test
                                        # dat$sex.test
                                        ))$y
  # coef.mod <- c(mod$alpha, mod$lambda)
  # names(coef.mod) <- c("alpha", "lambda")
  # print(coef.mod)
  mean(y.test != dat$y.test) * 100
}

# I-prior probit simulation
ipmySim <- function(y.mySim = y, X.mySim = X, sex.mySim = sex, nsim = 100,
                    n.mySim = n, kernel = c("Canonical", "FBM")) {
  kernel <- match.arg(kernel, c("Canonical", "FBM"))
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  progress <- function(i) setTxtProgressBar(pb, i)

  cl <- makeCluster(no.cores)
  registerDoSNOW(cl)
  res <- foreach(i = 1:nsim, .combine = rbind,
                 .packages = c("iprior", "iprobit"),
                 .export = c("testTrain", "probitInnerSim"),
                 .options.snow = list(progress = progress)) %dopar% {
      res.tmp <- rep(NA, length(n.mySim))
      for (j in 1:length(n.mySim)) {
        res.tmp[j]  <- probitInnerSim(y.mySim, X.mySim, sex.mySim, n.mySim[j],
                                      kernel = kernel)
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

  colnames(res) <- paste0(c("n = "), n.mySim)
  res
}

# Function to plot
plotRes <- function() {
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
  ggplot(plot.df, aes(x = value, y = id, col = id2, label = decPlac(value))) +
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
