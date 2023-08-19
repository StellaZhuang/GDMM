#######################################################################
# This script contains functions for univariate analysis
# and feature selection for multiple outcomes
# Univ: Fits linear/generalized linear models 
#       for each outcome in the given data
# Univ.feature: Conducts feature selection based on 
#               feature significance results from univariate analysis
#######################################################################

# Function Univ: Performs univariate analysis and model fitting
Univ = function(X, Y, outcome_type = c("binary", "continuous")) {
  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)
  X = as.matrix(X)
  
  res_fit = lapply(1:q, function(x) {
    if (outcome_type[x] == "continuous") {
      model = lm(Y[, x] ~ X) 
      return(model)
    } else if (outcome_type[x] == "binary") {
      model = glm(Y[, x] ~ X, family = binomial(link = "probit"))
      return(model)
    }
  })
  
  # Extract beta coefficients
  beta = t(as.data.frame(sapply(res_fit, `[`, "coefficients")))
  rownames(beta) = colnames(Y)
  colnames(beta) = c("Intercept", colnames(X))
  
  # Extract sigma & rho
  residual = sapply(1:q, function(x) {
    if (outcome_type[x] == "continuous") {
      residuals.lm(res_fit[[x]])
    } else if (outcome_type[x] == "binary") {
      residuals.glm(res_fit[[x]])
    }
  })
  
  sigma = apply(residual, 2, sd)
  names(sigma) = colnames(Y)
  
  rho = cor(residual)
  if (is.null(rownames(rho))) rownames(rho) = paste0('Y', 1:nrow(rho))
  if (is.null(colnames(rho))) colnames(rho) = paste0('Y', 1:ncol(rho))
  
  return(list(beta = beta, sigma = sigma, rho = rho))
}

# Function Univ.feature: Performs feature selection based on univariate analysis results
Univ.feature = function(X, Y, outcome_type = c("binary", "continuous"), select_type = c("u", "i")) { 
  # u: union, i: intersection
  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)
  X = as.matrix(X)
  
  res_fit = lapply(1:q, function(x) {
    if (outcome_type[x] == "continuous") {
      model = lm(Y[, x] ~ X)
      significant_vars <- colnames(X)[summary(model)$coefficients[-1, 4] < 0.05]
      return(significant_vars)
    } else if (outcome_type[x] == "binary") {
      model = glm(Y[, x] ~ X, family = binomial(link = "probit"))
      significant_vars <- colnames(X)[summary(model)$coefficients[-1, 4] < 0.05]
      return(significant_vars)
    }
  })
  
  # Select features based on select_type
  if (select_type == "u") {
    feature = unique(unlist(res_fit))
  } else {
    feature = Reduce(intersect, res_fit)
  }
  
  return(feature)
}
