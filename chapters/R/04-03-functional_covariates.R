# Chapter 4, Section 4
# Models with functional covariates (Tecator data set)
source("00-prelim.R")

## ---- tecator_data ----
data(tecator)
endpoints <- as.data.frame(endpoints)
colnames(endpoints) <- c("water", "fat", "protein")

# Plot of 10 random spectra predicting fat content
# (reproduced from package caret, v6.0-68)
set.seed(123)
inSubset <- sample(1:dim(endpoints)[1], 10)

absorpSubset <- absorp[inSubset,]
endpointSubset <- endpoints$fat[inSubset]

newOrder <- order(absorpSubset[,1])
absorpSubset <- absorpSubset[newOrder,]
endpointSubset <- endpointSubset[newOrder]

as.tibble(data.frame(id = factor(1:10), fat = endpointSubset, wave = absorpSubset)) %>%
  melt(id.vars = c("id", "fat")) %>%
  as.tibble() %>%
  mutate(x = rep(1:100, each = 10)) %>%
  ggplot(aes(x, value, col = id)) +
  geom_line() +
  geom_dl(aes(label = fat), method = list("last.polygons", cex = 0.9)) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = c(1, 25, 50, 75, 100)) +
  coord_cartesian(xlim = c(0, 105)) +
  labs(x = "Wavelength index", y = "Absorbance Units",
       title = "Percentage fat content predictor profiles for 10 random samples")

# Prepare data for iprior
# n = 215, use first 160 for training
absorpTrain <- -t(diff(t(absorp)))	# this takes first differences using diff()
absorpTest <- absorpTrain[161:215,]
absorpTrain <- absorpTrain[1:160,]

# Other variables
fatTrain <- endpoints$fat[1:160]
fatTest <- endpoints$fat[161:215]
waterTrain <- endpoints$water[1:160]
waterTest <- endpoints$water[161:215]

## ---- tecator_1 ----
# Model 1: canonical RKHS (linear)
(mod1 <- kernL(y = fatTrain, absorpTrain))

## ---- tecator_2 ----
# Model 2: canonical RKHS (quadratic)
mod2 <- kernL(y = fatTrain, absorpTrain, absorpTrain ^ 2,
              model = list(order = c("1", "1^2")))
## ---- tecator_3 ----
# Model 3: canonical RKHS (cubic)
mod3 <- kernL(y = fatTrain, absorpTrain, absorpTrain ^ 2, absorpTrain ^ 3,
              model = list(order = c("1", "1^2", "1^3")))

## ---- tecator_4 ----
# Model 4: FBM RKHS (default Hurst = 0.5)
mod4 <- kernL(y = fatTrain, absorpTrain, model = list(kernel = "FBM"))

## ---- tecator_5 ----
# Model 5: FBM RKHS + extra covariate
(mod5 <- kernL(y = fatTrain, absorpTrain, waterTrain,
               model = list(kernel = c("FBM", "Canonical"))))

## ---- tecator_fit ----
mod1.fit <- ipriorOptim(mod1, control = list(silent = TRUE))  # linear
mod2.fit <- ipriorOptim(mod2, control = list(silent = TRUE))  # quadratic
mod3.fit <- ipriorOptim(mod3, control = list(silent = TRUE))  # cubic
mod4.fit <- ipriorOptim(mod4, control = list(silent = TRUE))  # smooth
mod5.fit <- fbmOptim(mod4, silent = TRUE)  # smooth, MLE
mod6.fit <- fbmOptim(mod5, silent = TRUE)  # smooth, MLE with extra covariate

RMSE.Train1 <- mod1.fit$sigma
fatTestPredicted1 <- predict(mod1.fit, list(absorpTest))
RMSE.Test1 <- sqrt(mean((fatTestPredicted1 - fatTest) ^ 2))

RMSE.Train2 <- mod2.fit$sigma
fatTestPredicted2 <- predict(mod2.fit, list(absorpTest, absorpTest ^ 2))
RMSE.Test2 <- sqrt(mean((fatTestPredicted2 - fatTest) ^ 2))

RMSE.Train3 <- mod3.fit$sigma
fatTestPredicted3 <- predict(mod3.fit,
                             list(absorpTest, absorpTest ^ 2, absorpTest ^ 3))
RMSE.Test3 <- sqrt(mean((fatTestPredicted3 - fatTest) ^ 2))

RMSE.Train4 <- mod4.fit$sigma
fatTestPredicted4 <- predict(mod4.fit, list(absorpTest))
RMSE.Test4 <- sqrt(mean((fatTestPredicted4 - fatTest) ^ 2))

RMSE.Train5 <- mod5.fit$sigma
fatTestPredicted5 <- predict(mod5.fit, list(absorpTest))
RMSE.Test5 <- sqrt(mean((fatTestPredicted5 - fatTest) ^ 2))

RMSE.Train6 <- mod6.fit$sigma
fatTestPredicted6 <- predict(mod6.fit, list(absorpTest, matrix(waterTest, ncol = 1)))
RMSE.Test6 <- sqrt(mean((fatTestPredicted6 - fatTest) ^ 2))

## ---- tecator_compare ----
tab <- c(mod1.fit$log.lik, RMSE.Train1, RMSE.Test1)
tab <- rbind(tab, c(mod2.fit$log.lik, RMSE.Train2, RMSE.Test2))
tab <- rbind(tab, c(mod3.fit$log.lik, RMSE.Train3, RMSE.Test3))
tab <- rbind(tab, c(mod4.fit$log.lik, RMSE.Train4, RMSE.Test4))
tab <- rbind(tab, c(mod5.fit$log.lik, RMSE.Train5, RMSE.Test5))
tab <- rbind(tab, c(mod6.fit$log.lik, RMSE.Train6, RMSE.Test6))
rownames(tab) <- c("Linear", "Quadratic", "Cubic", "FBM gam=0.50",
                   paste0("FBM gam=", mod5.fit$ipriorKernel$model$Hurst),
                   paste0("FBM gam=", mod6.fit$ipriorKernel$model$Hurst[1],
                          " + extra cov."))
colnames(tab) <- c("Log-lik", "Training RMSE", "Test RMSE")
tab

# Likelihood ratio test
# G <- 2 * (mod6.fit$log.lik - mod5.fit$log.lik)
# (pval <- pchisq(G, 1))
