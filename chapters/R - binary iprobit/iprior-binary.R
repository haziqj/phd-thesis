## ---- prelim ----
library(iprior)
library(iprobit)
library(ggplot2)
library(progress)
library(gridExtra)


## ---- pmf.pdf ----
# f(x)
fiprior <- function(theta, w) {
  alpha <- theta[1]
  lambda <- theta[2]
  as.numeric(alpha + lambda * H %*% w)
}
# f(y | w)
l <- function(theta, w) {
  sum(
    y * pnorm(fiprior(theta, w), log.p = TRUE) +
      (1 - y) * pnorm(-fiprior(theta, w), log.p = TRUE)
  )
}

## ---- variational.bayes ----
library(iprobit)

## ---- iris.data ----
data(iris)
str(iris, strict.width = "cut", width = 70)
y <- ifelse(iris$Species == "setosa", 1, 0)
X <- iris[, -(3:5)]
n <- length(y)
H <- fnH2(X)  # canonical kernel
setosa <- as.factor(y)
levels(setosa) <- c("Setosa", "Others")

## ---- iris.plot1 ----
ggplot(data = cbind(X, Class = setosa),
       aes(x = Sepal.Length, y = Sepal.Width, col = Class)) +
  geom_point(size = 3)

## ---- iris.res ----
system.time(mod <- iprobit(y, X, silent = TRUE, maxit = 100))
print(mod)

## ---- iris.plot2 ----
plot(mod, 30, levels = c("Setosa", "Others"))

## ---- ionosphere.data ----
# n = 350, p = 34
ion <- read.table("ionosphere.data.txt", sep = ",", header = TRUE)
summary(ion$g)
X <- as.matrix(ion[, -35])
y <- as.numeric(ion$g)
y[y == 2] <- 0  # convert good = 0
train.index <- sample(1:length(y), 200)
test.index <- (1:length(y))[-train.index]
X.train <- X[train.index, ]
y.train <- y[train.index]
X.test <- X[test.index, ]
y.test <- y[test.index]

## ---- ionosphere.res ----
mod <- iprobit(y.train, X.train, kernel = "FBM", silent = TRUE)
print(mod)
print(mod, X.test, y.test)  # Test error rate

## ---- ionosphere.plot ----
plot(mod, 15, levels = c("good", "bad"))

## ---- cardiac.data ----
# n = 451, p = 194
load("Arrh194.RData")
tmp <- as.factor(ArrhDataNew$y)
levels(tmp) <- c("Normal", "Arrhythmia")
summary(tmp)
train.index <- sample(1:length(ArrhDataNew$y), 300)
test.index <- (1:length(ArrhDataNew$y))[-train.index]
X.train <- ArrhDataNew$x[train.index, ]
y.train <- ArrhDataNew$y[train.index] - 1
X.test <- ArrhDataNew$x[test.index, ]
y.test <- ArrhDataNew$y[test.index] - 1

## ---- cardiac.res ----
mod <- iprobit(y.train, X.train, kernel = "FBM", silent = F, maxit = 100)
print(mod)
print(mod, X.test, y.test)  # Test error rate

## ---- cardiac.plot ----
plot(mod, 15, levels = c("Normal", "Arrhythmia"))

