#############################################################
# This script contains all supporting functions
# bootsample: generate bootstrap sample datasets
# randomfeat: random sampling for feature space
# gdist: calculate Gower distance
# wilcox: perform Wilcoxon test for GDMM-random forest
# pairfun: Generate a sequence of functions for reference
# in the Newton-Raphson algorithm
#############################################################

# Function bootsample: Generates bootstrap sample datasets
bootsample = function(N = N, boots, samplesize = 1.0, seed = seed) { 
  # boots: number of bootstrap datasets
  # N: size of bootstrap datasets
  
  if (is.na(boots)) {
    boots = 100
  }
  
  samindexlist = lapply(1:boots, function(x) {
    set.seed(x * seed)
    sam = sample(1:N, N, replace = TRUE)
    sam = sort(sam)
  })
  
  return(samindexlist)
}

# Function randomfeat: Randomly samples features from the dataset
randomfeat = function(p = p, maxfeatsample = ncol(traindf) - length(outvar), minfeatsample = 2, boots = 1, seed = seed) {
  if (is.na(boots)) {
    boots = 100
  }
  
  featindexlist = lapply(1:boots, function(x) {
    set.seed(x * seed)
    sam = sample(1:p, sample(minfeatsample:maxfeatsample, 1), replace = FALSE)
    sam = sort(sam)
  })
  
  return(featindexlist)
}

# Function gdist: Calculates Gower distance between two datasets
gdist <- function(X, Y, est_beta, outcome_type) {
  m <- cbind(1, X)
  m = as.matrix(m)
  pred <- m %*% t(est_beta)
  
  for (i in 1:ncol(pred)) {
    if (outcome_type[i] == "binary") {
      pred[, i] = pnorm(0, mean = pred[, i], sd = 1, lower.tail = FALSE)
    }
  }
  
  colnames(pred) = colnames(Y)
  g = sum(gower_dist(as.data.frame(Y), as.data.frame(pred)))
  
  return(g)
}

# Function wilcox: Performs Wilcoxon test for GDMM-random forest
wilcox <- function(i, S, M, p) {
  shap1 = S[which(M[, i] == 1), i]
  shap0 = S[which(M[, i] == 0), i]
  w.test = wilcox.test(shap1, shap0, alternative = "less")
  
  return(ifelse(w.test$p.value < 0.05 / p, TRUE, FALSE)) 
}

# Function pair: Generates a matrix of pairwise outcome comparisons
pair = function(outcome_type) {
  q = length(outcome_type)
  pairfun = matrix(nrow = q, ncol = q)
  
  for (i in 1:q - 1) {
    for (j in (i + 1):q) {
      pairfun[i, j] = (outcome_type[i] == "binary") + (outcome_type[j] == "binary")
    }
  }
  
  pairfun = as.data.frame(pairfun)
  pairfun = pairfun %>%
    dplyr::mutate(across(everything(),
                         ~case_when(.x == 2 ~ "Case_2_discrete",
                                    .x == 1 ~ "Case_1_discrete",
                                    .x == 0 ~ "Case_2_continuous")))
  
  return(pairfun)
}
