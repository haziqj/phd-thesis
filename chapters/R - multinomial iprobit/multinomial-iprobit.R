## ---- prelim ----
library(iprobit)
library(iprior)
library(ggplot2)
library(kableExtra)

## ---- iris.data ----
data(iris)
y <- iris$Species
X <- iris[, -5]
ggplot(data = iris, aes(x = Sepal.Length, y = Sepal.Width, color = Species)) +
  geom_point(size = 3) + theme_bw()
ggplot(data = iris, aes(x = Petal.Length, y = Petal.Width, color = Species)) +
  geom_point(size = 3) + theme_bw()
ggplot(data = iris, aes(x = Sepal.Length, y = Petal.Length, color = Species)) +
  geom_point(size = 3) + theme_bw()
ggplot(data = iris, aes(x = Sepal.Width, y = Petal.Width, color = Species)) +
  geom_point(size = 3) + theme_bw()
ggplot(data = iris, aes(x = Sepal.Length, y = Petal.Width, color = Species)) +
  geom_point(size = 3) + theme_bw()
ggplot(data = iris, aes(x = Petal.Length, y = Sepal.Width, color = Species)) +
  geom_point(size = 3) + theme_bw()

## ---- iris.mod ----
(mod <- iprobit_mult(y, X, silent = TRUE))
plot(mod)
iplot_lb(mod)

## ---- toy.data ----
genSample <- function(n, noiseVar = 0) {
  ## Class 1 and 2 (x ~ U(0,1))
  u <- 4. * matrix(runif(2*n), nrow = n, ncol = 2) - 2.
  i <- which(((u[, 1]^2 + u[, 2]^2) > .1) & ((u[, 1]^2 + u[, 2]^2) < .5))
  j <- which(((u[, 1]^2 + u[, 2]^2) > .6) & ((u[, 1]^2 + u[, 2]^2) < 1))
  X <- u[c(i, j), ]
  t.class <- c(rep(1, length(i)), rep(2, length(j)))
  ## Class 3 (x ~ N(0,1))
  x <- 0.1 * matrix(rnorm(2 * length(i)), ncol = 2, nrow = length(i) )
  k <- which((x[, 1]^2 + x[, 2] ^ 2) < 0.1)
  X <- rbind(X, x[k, ])
  t.class <- c(t.class, rep(3, length(k)))
  ## Add random columns
  if (noiseVar > 0) X <- cbind(X, matrix(rnorm(noiseVar * nrow(X)),
                                         ncol = noiseVar, nrow = nrow(X)))
  list(class = as.factor(t.class), X = X)
}
set.seed(123)
toy.train <- genSample(2200)
toy.test  <- genSample(2200)
df.plot <- data.frame(X = toy.train$X, Class = toy.train$class)
colnames(df.plot) <- c("X1", "X2", "Class")
ggplot(df.plot, aes(x = X1, y = X2, col = Class)) +
  geom_point() + theme_bw()

## ---- toy.mod ----
mod <- iprobit_mult(toy.train$class, toy.train$X, kernel = "FBM", maxit = 10)
plot(mod)
print(mod, toy.test$X, toy.test$class)

## ---- vowel.data ----
library(ElemStatLearn)
data("vowel.test")
data("vowel.train")
vow.tr <- vowel.train
vow.ts <- vowel.test[sample(1:nrow(vowel.test), 100), ]
vow.tr$y <- as.factor(vow.tr$y)
vow.ts$y <- as.factor(vow.ts$y)
names(vow.tr)[1] <- names(vow.ts)[1] <- "class"
head(vow.tr)
vow.tr <- vow.tr[order(vow.tr[, 1]), ]

## ---- vowel.mod ----
set.seed(123)
(mod <- iprobit_mult(vow.tr$class, vow.tr[, -1], kernel = "FBM", silent = TRUE))
plot(mod)
predict(mod, X.test = vow.ts[, -1], y.test = vow.ts[, 1])

## ---- vowel.tab ----
# tmp <- predict(mod, X.test = vow.ts[, -1], y.test = vow.ts[, 1])
iprobit.fbm <- c(0, 41) # tmp$test
tab <- rbind(
  "k-Nearest neighbours" = c(NA, 44),
  "Linear regression" = c(48, 67),
  "Linear discriminant analysis" = c(32,56),
  "Neural network" = c(NA, 45),
  "FDA/BRUTO" = c(6, 44),
  "FDA/MARS" = c(13, 39),
  "I-probit (FBM-0.5)" = c(0, 41)
)
colnames(tab) <- c("Training", "Test")
knitr::kable(tab, format = "latex", booktabs = TRUE) %>%
  kableExtra::kable_styling(position = "center") %>%
  kableExtra::add_header_above(c(" " = 1, "Error rates" = 2))


