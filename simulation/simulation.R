# Set working directory and random seed
setwd(code_path)
set.seed(seed)

# Load necessary packages
library(plyr)
library(dplyr)
library(glmnet)
library(hydroGOF)
library(future.apply)
library(MASS)
library(OpenMx)
library(gower)
library(progressr)
library(MASS)
library(mvtnorm)
library(statmod)
library(lme4)
library(aSPU)
library(mosaic)
library(daarem)
library(pROC)
library(tools)
library(e1071)
library(caret)
library(randomForest)
library(stats)
library(parallel)
library(kernelshap)
library(doFuture)

# Set plot device
pdf(file = NULL)

# Read data and extract parameter values
df <- readRDS(file = paste0(data_path, "/data.RData"))

# Extract training, test set, and parameter values from the given data
traindf <- df$train
testdf <- df$test
ta <- rbind(traindf, testdf)
outcome_type <- df$info$outcome_type
diffcut <- 1 / N
M <- B + 50
newton <- 1

# Create pair matrix
pairlist <- pair(outcome_type)  # matrix

# Load functions
pairfunlist <- c(
  Case_2_discrete = Case_2_discrete,
  Case_1_discrete = Case_1_discrete,
  Case_2_continuous = Case_2_continuous
)

# Generate training and test set
index <- sample(1:nrow(ta), N, replace = FALSE)
traindf <- ta[index, ]
testdf <- ta[-index, ]
X <- traindf[, 1:p]
Y <- traindf[, (p + 1):(p + q)]
test.X <- testdf[, 1:p]
test.Y <- testdf[, (p + 1):(p + q)]

# Bootstrap sampling
samplelist <- bootsample(N = N, boots = M, samplesize = 1.0, seed = seed)
featlist <- randomfeat(p = p, maxfeatsample = maxfeatsample, minfeatsample = minfeatsample, boots = M, seed = seed)

# Run the model
fitlist <- future.apply::future_lapply(1:M, function(x) {
  sampleindex <- samplelist[[x]]
  featureindex <- featlist[[x]]
  newtraindfX <- X[sampleindex, featureindex]
  newtraindfY <- Y[sampleindex, ]
  if (nrow(X[-sampleindex, ]) == 0) print("No valid test samples!")
  res_sub <- Univ(X = newtraindfX, newtraindfY, outcome_type)
  beta_0 <- res_sub[[1]]
  sigma_0 <- res_sub[[2]]
  rho_0 <- diag(q)
  fit <- MMRM(newtraindfX, newtraindfY, beta_0, sigma_0, rho_0, pairlist, outcome_type, diffcut, newton = newton)
  if (length(fit) == 0) return()
  beta_v <- rep(0, p)
  names(beta_v) <- colnames(beta_0)[-1]
  beta_v[featureindex] <- 1
  est_beta <- fit$beta
  newvalidfX <- X[-sampleindex, featureindex]
  newvalidfY <- Y[-sampleindex, ]
  g <- gdist(newvalidfX, newvalidfY, est_beta, outcome_type) / nrow(newvalidfX)
  
  return(list(fit = fit, gsample <- c(beta_v, g)))
}, future.seed = TRUE)

# Fit a model for feature performance and model performance
fitlist <- fitlist[!sapply(fitlist, is.null)]
pm <- sapply(1:B, function(x) fitlist[[x]][[2]])
pm <- t(pm)
colnames(pm) <- c(colnames(X), "g_dist")

# Supervised learning
if (method == "LA") {
  grid <- 10^seq(10, -2, length = 100)
  cv.out <- cv.glmnet(as.matrix(pm[, -ncol(pm)]), as.matrix(pm[, ncol(pm)]), alpha = 1, lambda = grid, nfolds = 10)
  lambda_selected <- cv.out$lambda.min
  lasso.mod <- glmnet(pm[, -ncol(pm)], pm[, ncol(pm)], alpha = 1, lambda = lambda_selected)
  coef <- coef(lasso.mod, s = "lambda.min")
  feature <- coef@Dimnames[[1]][which(coef < 0)]
  feature <- feature[feature != "(Intercept)"]
} else if (method == "EL") {
  elastic <- train(g_dist ~ ., data = pm, method = "glmnet", trControl = trainControl("cv", number = 5), tuneLength = 10)
  alpha <- elastic$bestTune[1]
  lambda <- elastic$bestTune[2]
  Matrix <- model.matrix(g_dist ~ ., as.data.frame(pm))[, -1]
  model <- glmnet(Matrix, pm[, ncol(pm)], lambda = lambda, alpha = alpha, standardize = FALSE)
  coef <- as.matrix(coef(model))
  feature <- rownames(coef)[coef < 0]
  feature <- feature[feature != "(Intercept)"]
} else if (method == "RF") {
  bestMtry <- tuneRF(pm[, -ncol(pm)], pm[, ncol(pm)], stepFactor = 1.5, improve = 1e-5, ntree = 1000)
  bmtry <- as.numeric(names(which(bestMtry[, 2] == min(bestMtry[, 2]))))
  rf.fit <- randomForest(g_dist ~ ., data = pm, ntree = 1000, keep.forest = TRUE, importance = TRUE, mtry = bmtry)
  imp <- scale(importance(rf.fit, typ = 1))
  feature <- colnames(pm)[which(imp > 0)]
  pm <- as.data.frame(pm)
  registerDoFuture()
  plan(multicore, workers = availableCores())
  s <- kernelshap(object = rf.fit, X = pm[, -ncol(pm)], bg_X = pm, feature_names = feature, parallel = TRUE)
  S <- s$S
  M <- pm[, which(imp > 0)]
  stat <- sapply(1:length(feature), wilcox, S = S, M = M, p = length(feature))
  feature <- feature[which(stat == TRUE)]
}
  
if (length(feature) != 0) {
  # Select features based on the model
  X.select <- X[, colnames(X) %in% feature, drop = FALSE]
  res_select <- Univ(X = X.select, Y, outcome_type)
  beta_select <- res_select[[1]]
  sigma_select <- res_select[[2]]
  rho_select <- diag(q)
  fit.select <- MMRM(X.select, Y, beta_select, sigma_select, rho_select, pairlist, outcome_type, diffcut, newton = newton)
  model.beta <- fit.select$beta
  test.X.select <- test.X[, colnames(X) %in% feature]
  pred.Y <- cbind(1, test.X.select) %*% t(model.beta)
  
  # Evaluate performance of the selected model
  performance <- c()
  for (a in 1:q) {
    if (outcome_type[a] == "continuous") {
      performance <- c(performance, mse(test.Y[, a], pred.Y[, a]))
    } else if (outcome_type[a] == "binary") {
      perf.roc <- roc(test.Y[, a], pnorm(pred.Y[, a]))
      performance <- c(performance, auc(perf.roc))
    }
  }
} else {
  performance <- "No performance"
}

# Evaluate feature selection performance
allvar <- paste("X", 1:p, sep = "")
truevar <- paste("X", 1:truep, sep = "")
res.true <- allvar %in% truevar
res.select <- allvar %in% feature

res.select <- factor(res.select, levels = c(FALSE, TRUE))

# Create a contingency table and calculate performance metrics
res.table <- table(res.select, res.true)
s.t <- res.table[2, 2]
s.tp <- res.table[2, 2] / truep
s.fp <- res.table[2, 1] / (p - truep)
s.tn <- res.table[1, 1] / (p - truep)
s.fn <- res.table[1, 2] / truep
s.f <- res.table[2, 1]

# Calculate rate metrics
rate <- c(s.tp, s.tn, s.fp, s.fn)

# Write results to files
write.table(feature, file = paste0(result_path, "/feature_rep", seed, ".txt"), row.names = FALSE, sep = ",")
write.table(rate, file = paste0(result_path, "/rate_rep", seed, ".txt"), row.names = FALSE, sep = ",")
write.table(performance, file = paste0(result_path, "/performance_rep", seed, ".txt"), row.names = FALSE, sep = ",")
