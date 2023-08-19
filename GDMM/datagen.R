##############################################################################
# This function generates simulated data for training and testing
# Input:
#   train_sample, test_sample: sample size for the training and test sets
#   correlation_var: number of features with correlations
#   correlation_val: correlation coefficient among variables
#   varnum: number of features
#   outcome_type_gen: a vector indicating outcome types (binary first)
#   q: number of outcomes
#   seed: random seed
#   setting: feature correlation setting (No_correlation/fixed)
#   var_con: continuous outcome variance
#   rho: correlation among outcomes
# Output:
#   A list with three elements:
#   train: the generated training set
#   test: the generated test set
#   info: a list containing all settings of the generated dataset
##############################################################################


Datagen <- function(train_sample, test_sample, correlation_var, correlation_val, varnum, outcome_type_gen, q, seed, setting, var_con, rho) {
  
  # Generate the correlation matrix for X
  if (setting == "No_Correlation") {
    Sigma = matrix(rep(0, varnum), nrow = varnum, ncol = varnum, byrow = FALSE)
    diag(Sigma) <- 10
  } else if (setting == "fixed") {
    Sigma = matrix(rep(0, varnum), nrow = varnum, ncol = varnum, byrow = FALSE)
    for (i in 1:varnum) {
      for (j in 1:varnum) {
        if (i == j) {
          Sigma[i, j] = 10
        } else if (i <= correlation_var & j <= correlation_var) {
          Sigma[i, j] = Sigma[j, i] = correlation_val
        } else {
          Sigma[i, j] = Sigma[j, i] = 0
        }
      }
    }
  } else {
    Sigma = randcorr(varnum = varnum, maxcorr = correlation_val, type = c("matrix"), seed = seed)
  }
  
  # Set random seed
  set.seed(seed)
  
  # Generate X 
  ta = data.frame(MASS::mvrnorm(n = train_sample + test_sample, rep(0, varnum), Sigma / 10))
  
  variablelist = list()
  for (i in 1:varnum) {
    variablelist[[i]] = gsub(" ", "", paste("X", i))
    ta[, i] = mosaic::zscore(ta[, i])
  }
  variablelist = unlist(variablelist)
  colnames(ta) = variablelist
  
  # Create the outcome matrix Y
  
  # Generate formula
  mar_var = paste(colnames(ta), collapse = "+")
  f = as.formula(paste("~", mar_var))
  
  main_mat = model.matrix(f, ta) # Design matrix
  
  # Generate variance-covariance matrix
  covmatrix = matrix(0, nrow = q, ncol = q)
  n.bin = sum(outcome_type_gen == "binary")
  diag(covmatrix) = c(rep(1, n.bin), var_con) # Variance
  covmatrix[col(covmatrix) != row(covmatrix)] = sqrt(var_con) * rho
  
  error = data.frame(MASS::mvrnorm(n = train_sample + test_sample, rep(0, q), covmatrix))
  
  Y = matrix(nrow = train_sample + test_sample, ncol = q)
  for (i in 1:q) {
    if (outcome_type_gen[i] == 'continuous') {
      Y[, i] = main_mat %*% var_effect_all[, i] + error[, i]
    } else {
      Y[, i] = ProbitSimulate(c(var_effect_all[, i], 1), cbind(main_mat, error[, i]))
    }
  }
  colnames(Y) = paste("Y", 1:q, sep = "")
  
  # Generate training and test sets
  ta = cbind(ta, Y)
  ta = as.matrix(ta)
  index = sample(1:nrow(ta), train_sample, replace = FALSE)
  traindf = ta[index, ]
  validationdf = ta[-index, ]
  
  info <- list(train_sample = train_sample, test_sample = test_sample, correlation_var = correlation_var,
               correlation_val = correlation_val, varnum = varnum, outcome_type_gen = outcome_type_gen,
               q = q, seed = seed, setting = setting, var_con = var_con, rho = rho)
  
  return(list(train = traindf, test = validationdf, info = info))
}