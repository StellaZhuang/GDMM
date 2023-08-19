#################################################
# This function is used for multivariate 
# mixed response model likelihood estimation
#################################################
MMRM <- function(X_sub, Y_s, beta_0, sigma_0, rho_0, pairlist, outcome_type, diffcut = 1/N, newton = 1, firststep = 1e-5) {
  rhos_v <- rho_0[lower.tri(rho_0)]
  eta_0 <- as.numeric(c(t(beta_0), sigma_0, rhos_v))
  p <- ncol(X_sub)
  q <- ncol(Y_s)
  
  # Computation Process
  L_return = L_return0 = LL_ha(X_sub, Y_s, beta_0, sigma_0, rho_0, pairlist, outcome_type)
  if (length(L_return) == 0) return()
  
  eta_old = gradient_old = list()
  diff <- 1
  rua <- 1
  
  if (newton == 0) {
    minrua = 4
  } else {
    minrua = 3
  }
  
  while ((diff > diffcut | rua < minrua) & rua <= 100) {
    if (newton == 1) {
      eta_new <- eta_0 - ginv(L_return[[3]]) %*% L_return[[1]]
    } else {
      if (rua <= 2) {
        ga = firststep
        eta_new = eta_0 - ga * L_return0[[1]]
      } else {
        q1 = eta_old[[rua - 1]] - eta_old[[rua - 2]]
        q2 = gradient_old[[rua - 1]] - gradient_old[[rua - 2]]
        ga = t(q1) %*% q2 / (t(q2 %*% q2))
        eta_new = eta_0 - c(ga) * c(gradient_old[[rua - 1]])
      }
    }
    
    diff <- sum(abs(eta_new - eta_0))
    eta_0 <- eta_new
    
    beta_new = eta_new[1:length(beta_0)]
    beta_new = matrix(beta_new, byrow = TRUE, nrow = nrow(beta_0))
    
    sigma_new = eta_new[(length(beta_0) + 1):(length(beta_0) + q)]
    sigma_new = ifelse(sigma_new > 0, sigma_new, 0)
    sigma_new[sigma_new == 0] = 0.9 * sigma_0[sigma_new == 0]
    
    rho_new = matrix(nrow = q, ncol = q)
    diag(rho_new) = 1
    rho_new[lower.tri(rho_new)] = eta_new[(length(beta_0) + q + 1):length(eta_new)]
    rho_new[upper.tri(rho_new)] = t(rho_new)[upper.tri(rho_new)]
    rho_new[rho_new > 1] = 0.9
    rho_new[rho_new < (-1)] = -0.9
    
    L_return <- LL_ha(X_sub, Y_s, beta_new, sigma_new, rho_new, pairlist, outcome_type)
    if (length(L_return) == 0) return()
    
    eta_old[[rua]] = eta_new
    gradient_old[[rua]] = L_return[[1]]
    rua = rua + 1
  }
  
  est_beta = matrix(eta_new[1:((p + 1) * q)], nrow = q, byrow = TRUE)
  dimnames(est_beta) = dimnames(beta_0)
  
  sig_index = which(outcome_type != "binary")
  est_sigma = eta_new[(p + 1) * q + sig_index]
  names(est_sigma) = paste("sigma", sig_index, sep = "")
  
  est_rho = matrix(nrow = q, ncol = q)
  diag(est_rho) = 1
  est_rho[lower.tri(est_rho)] = eta_new[((p + 2) * q + 1):length(eta_new)]
  est_rho[upper.tri(est_rho)] = t(est_rho)[upper.tri(est_rho)]
  
  sd = L_return[[4]]
  est_list = list(beta = est_beta, sigma = est_sigma, rho = est_rho, sd = sd)
}
