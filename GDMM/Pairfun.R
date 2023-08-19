#####################################################################
# This script contains functions for pairwise likelihood estimation
# Case_2_discrete: Pair of two binary outcomes
# Case_2_continuous: Pair of two continuous outcomes
# Case_1_discrete: Pair of one binary and one continuous outcome
#####################################################################

# Function Case_2_discrete: Pair of two binary outcomes
Case_2_discrete <- function(data_1, data_2, beta_1, beta_2, sig_1, sig_2, rho_12 ) {
  N_1 <- length(beta_1) #beta_1: intercept + coefficient
  N_2 <- length(beta_2)
  N<-nrow(data_1)

  rho <- rho_12

  Y_1<- data_1[,ncol(data_1)] #outcome j column
  Y_2<- data_2[,ncol(data_2)] #outcome k column
  X_1<- as.matrix (data_1[,-ncol(data_1)]) # intercept and X column
  X_2<- as.matrix (data_2[,-ncol(data_2)])

  mu_1 <- X_1 %*% beta_1
  mu_2 <- X_2 %*% beta_2
  q_1 <- 2 * Y_1 - 1     #convert 0,1 to -1,1
  q_2 <- 2 * Y_2 - 1
  rho_tulda <- q_1 * q_2 * rho
  w_1 <- q_1 * mu_1
  w_2 <- q_2 * mu_2

  Q_1 <- w_1/sig_1
  Q_2 <- w_2/sig_2
  delta <- 1/sqrt(1- rho_tulda^2)
  V_1 <- delta * (Q_2 - rho_tulda * Q_1)
  V_2 <- delta * (Q_1 - rho_tulda * Q_2)
  G_1 <- dnorm(Q_1) * pnorm(V_1)
  G_2 <- dnorm(Q_2) * pnorm(V_2)
  sigma_tulda <- array(0, dim=c(N, 2, 2))

  phi_2 <- rep(0,N)
  Phi_2 <- rep(0,N)
  for (ids in c(1:N))
  {
    sigma_tulda[ids,,] <- matrix (c(1,rho_tulda[ids] ,rho_tulda[ids], 1 ),nrow=2,  ncol=2,  byrow = TRUE)
    phi_2[ids] <- mvtnorm::dmvnorm( c(Q_1[ids] , Q_2[ids] ), mean = rep(0, 2), sigma = sigma_tulda[ids,,] )
    Phi_2[ids] <- mvtnorm::pmvnorm(lower=-Inf, upper= c(Q_1[ids] , Q_2[ids] ) , mean=rep(0, 2),  sigma=sigma_tulda[ids,,] )
  }
   if (any(Phi_2==0)==T) return()
  A_beta_1 <- - ( Q_1 * G_1 / Phi_2 + rho_tulda * phi_2/Phi_2 + G_1^2 / Phi_2^2   ) / sig_1 ^2
  A_beta_2 <- - ( Q_2 * G_2 / Phi_2 + rho_tulda * phi_2/Phi_2 + G_2^2 / Phi_2^2   ) / sig_2 ^2
  A_beta_12 <- ( phi_2/Phi_2 - G_1 * G_2 / Phi_2^2 )/(sig_1 * sig_2 )

  Arr_beta_1 <- array(0, dim = c(N_1,N_1,N))
  Arr_beta_2 <- array(0, dim = c(N_2,N_2,N))
  Arr_beta_12 <- array(0, dim = c(N_1,N_2,N))


  d_l_d_beta_1 <- t(X_1) %*% ( q_1/sig_1 * G_1 /Phi_2 )
  d_l_d_beta_2 <- t(X_2) %*% ( q_2/sig_2 * G_2 /Phi_2 )
  d_l_d_rho <- sum (q_1 * q_2 * phi_2 / Phi_2 )

  d_l_d_beta_1_point <- X_1 * rep ( q_1/sig_1 * G_1 /Phi_2 , N_1 )
  d_l_d_beta_2_point <- X_2 * rep ( q_2/sig_2 * G_2 /Phi_2 , N_2 )
  d_l_d_rho_point <- q_1 * q_2 * phi_2 / Phi_2


  for (ids in c(1:N))
  {
    Arr_beta_1[,,ids] <- X_1[ids,] %*% t(X_1[ids,]) * A_beta_1 [ids]
    Arr_beta_2[,,ids] <- X_2[ids,] %*% t(X_2[ids,]) * A_beta_2 [ids]
    Arr_beta_12[,,ids] <- X_2[ids,] %*% t(X_1[ids,]) * A_beta_12 [ids]
  }
  d2_l_d_beta1_beta1 <- matrix(0,N_1,N_1)
  for(j in 1:N_1){
    for(i in 1:N_1){
      d2_l_d_beta1_beta1[i,j] <- sum(Arr_beta_1[j,i,])
    }
  }
  d2_l_d_beta2_beta2 <- matrix(0,N_2,N_2)
  for(j in 1:N_2){
    for(i in 1:N_2){
      d2_l_d_beta2_beta2[i,j] <- sum(Arr_beta_2[j,i,])
    }
  }
  d2_l_d_beta1_beta2 <- matrix(0,N_1,N_2)
  for(j in 1:N_2){
    for(i in 1:N_1){
      d2_l_d_beta1_beta2[i,j] <- sum(Arr_beta_12[j,i,])
    }
  }

  d2_l_d_beta1_rho <-   t(X_1) %*%  ( q_2 /sig_1 * phi_2/Phi_2* (rho_tulda * delta * V_1 - Q_1 -  G_1 / Phi_2))
  d2_l_d_beta2_rho <-   t(X_2) %*%  ( q_1 /sig_2 * phi_2/Phi_2* (rho_tulda * delta * V_2 - Q_2 -  G_2 / Phi_2))
  d2_l_d_rho_rho <- sum( phi_2/Phi_2 * ( delta^2 * rho_tulda * ( 1 - delta^2 * (Q_1^2 + Q_2^2 - 2*rho_tulda* Q_1 * Q_2 )) + delta^2*Q_1*Q_2 -phi_2/Phi_2 ))



  grad<- c(t(d_l_d_beta_1), t(d_l_d_beta_2), d_l_d_rho )
  grad_point<- cbind(d_l_d_beta_1_point, d_l_d_beta_2_point, d_l_d_rho_point )

  hessian<-matrix(0,(N_1 + N_2 + 1),(N_1 + N_2 + 1))


  hessian [1:N_1,1:N_1] <- d2_l_d_beta1_beta1
  hessian [(N_1+1):(N_1+N_2),(N_1+1):(N_1+N_2)] <- d2_l_d_beta2_beta2
  hessian [1:N_1,(N_1+1):(N_1+N_2)] <- d2_l_d_beta1_beta2
  hessian [(N_1+1):(N_1+N_2),1:N_1] <- d2_l_d_beta1_beta2
  hessian [(N_1+N_2+1),(N_1+N_2+1)] <- d2_l_d_rho_rho
  hessian [(1:N_1),(N_1+N_2+1)] <- d2_l_d_beta1_rho
  hessian [(N_1+N_2+1),1:N_1] <- d2_l_d_beta1_rho
  hessian [(N_1+1):(N_1+N_2),(N_1+N_2+1)] <- d2_l_d_beta2_rho
  hessian [(N_1+N_2+1),(N_1+1):(N_1+N_2)] <- d2_l_d_beta2_rho



  l_12<-0
  attr(l_12,"gradient")<- grad
  attr(l_12,"hessian")<- hessian
  attr(l_12,"grad_point")<- grad_point
  l_12

  return(l_12)
}

# Function Case_2_continuous: Pair of two continuous outcomes
Case_2_continuous<- function(data_1, data_2, beta_1, beta_2, sig_1, sig_2, rho_12) {

  N_3 <- length(beta_1)
  N_4<-length(beta_2)
  rho <- rho_12
  X_3 <- as.matrix (data_1[,-ncol(data_1)])
  X_4 <- as.matrix (data_2[,-ncol(data_2)])
  Y_3<-data_1[,ncol(data_1)]
  Y_4<-data_2[,ncol(data_2)]

  mu_3 = X_3 %*% beta_1
  mu_4 = X_4 %*% beta_2
  eta = 1/(1-rho^2) *( (Y_3 - mu_3 )^2/sig_1^2 + ( Y_4 - mu_4 )^2/sig_2^2  - 2 * rho*(Y_3 - mu_3) * (Y_4 - mu_4)/(sig_1 *sig_2)  )
  A_3=(Y_3-mu_3)/sig_1
  A_4=(Y_4-mu_4)/sig_2

  l_34.beta3<- t(X_3) %*% (1/(1-rho^2)*( ( Y_3 - mu_3 )/sig_1^2 - rho * (Y_4 - mu_4)/(sig_1 *sig_2) ))
  l_34.beta4<- t(X_4) %*% (1/(1-rho^2)*( ( Y_4 - mu_4 )/sig_2^2 - rho * (Y_3 - mu_3)/(sig_1 *sig_2) ))
  l_34.sig_1 = sum(-1/(sig_1)*( 1- 1/(1-rho^2)*(( Y_3 - mu_3)^2/sig_1^2 - rho*(Y_3 - mu_3)*(Y_4 - mu_4)/(sig_1 *sig_2) )) )
  l_34.sig_2 = sum(-1/(sig_2)*( 1- 1/(1-rho^2)*(( Y_4 - mu_4)^2/sig_2^2 - rho*(Y_3 - mu_3)*(Y_4 - mu_4)/(sig_1 *sig_2) )) )
  l_34.rho = sum( 1/(1-rho^2) * (rho - rho*eta + (Y_3 - mu_3)*(Y_4 - mu_4)/(sig_1 *sig_2) ) )

  l_34.beta3_point<- X_3 * rep ( (1/(1-rho^2)*( ( Y_3 - mu_3 )/sig_1^2 - rho * (Y_4 - mu_4)/(sig_1 *sig_2) )), N_3)
  l_34.beta4_point<- X_4 * rep ( (1/(1-rho^2)*( ( Y_4 - mu_4 )/sig_2^2 - rho * (Y_3 - mu_3)/(sig_1 *sig_2) )), N_4)
  l_34.sig_1_point = -1/(sig_1)*( 1- 1/(1-rho^2)*(( Y_3 - mu_3)^2/sig_1^2 - rho*(Y_3 - mu_3)*(Y_4 - mu_4)/(sig_1 *sig_2) ))
  l_34.sig_2_point = -1/(sig_2)*( 1- 1/(1-rho^2)*(( Y_4 - mu_4)^2/sig_2^2 - rho*(Y_3 - mu_3)*(Y_4 - mu_4)/(sig_1 *sig_2) ))
  l_34.rho_point =  1/(1-rho^2) * (rho - rho*eta + (Y_3 - mu_3)*(Y_4 - mu_4)/(sig_1 *sig_2) )


  # Hessian  ;


  l_34.beta_1beta_1<--1/((1-rho^2)*sig_1^2) * t(X_3) %*% X_3
  l_34.beta_2beta_2<--1/((1-rho^2)*sig_2^2) * t(X_4) %*% X_4
  l_34.beta_1beta_2<- rho/((1-rho^2)*sig_1* sig_2) * t(X_3) %*% X_4

  l_34.sig_1sig_1 <-  sum( 1/(sig_1^2) * ( ( 1- 1/(1-rho^2)*(( Y_3-mu_3)^2/sig_1^2 - rho*(Y_3 - mu_3)*(Y_4 - mu_4)/(sig_1 *sig_2) ) ) +
                                             1/(1-rho^2)*( - 2*(Y_3-mu_3)^2 /sig_1^2 + rho*(Y_3 - mu_3)*(Y_4 - mu_4)/(sig_1 *sig_2) )  )   )


  l_34.sig_2sig_2 <- sum( 1/(sig_2^2)* ( ( 1- 1/(1-rho^2)*(( Y_4 - mu_4)^2/sig_2^2 - rho*(Y_3 - mu_3)*(Y_4 - mu_4)/(sig_1 *sig_2) )) +
                                           1/(1-rho^2)*( - 2*(Y_4-mu_4)^2 /sig_2^2 + rho*(Y_3 - mu_3)*(Y_4 - mu_4)/(sig_1 *sig_2) )  )   )


  #l_34.rho_rho <- sum ( (1/(1-rho^2) + 2*rho^2/(1-rho^2)^2) * (1 - eta + (Y_3 - mu_3)*(Y_4 - mu_4)/(rho*sig_1 *sig_2) )  -
                         # rho/(1-rho^2)* (  2/(1-rho^2)* ( rho * eta - (Y_3 - mu_3)*(Y_4 - mu_4)/(sig_1 *sig_2) ) +
                                             # (Y_3 - mu_3)*(Y_4 - mu_4)/(rho^2 * sig_1 *sig_2) ) )
l_34.rho_rho<-sum((1+rho^2)/((1-rho^2)^2)+2*rho*A_3*A_4/((1-rho^2)^2)-(A_3^2+A_4^2)*(3*rho^2+1)/((1-rho^2)^3)+A_3*A_4*(4*rho*(rho^2+1))/((1-rho^2)^3))


  l_34.beta_1sig_1<- t(X_3) %*% (-1/(1-rho^2) * ( 2*(Y_3-mu_3)/sig_1^3 - rho * (Y_4-mu_4)/(sig_1^2 *sig_2) ) )
  l_34.beta_1sig_2<- t(X_3) %*% ( 1/(1-rho^2) * (  rho * (Y_4-mu_4)/(sig_1 *sig_2^2) ) )
  l_34.beta_1_rho<-  t(X_3) %*% ( ( 2*rho/(1-rho^2)^2 * ( (Y_3-mu_3)/sig_1^2 - rho * (Y_4-mu_4)/(sig_1 *sig_2) ) ) - 1/(1-rho^2) * (Y_4-mu_4)/(sig_1 *sig_2))

  l_34.beta_2sig_1<- t(X_4) %*% ( 1/(1-rho^2) * (  rho * (Y_3-mu_3)/(sig_1^2 *sig_2) ) )
  l_34.beta_2sig_2<- t(X_4) %*% (-1/(1-rho^2) * ( 2*(Y_4-mu_4)/sig_2^3 - rho * (Y_3-mu_3)/(sig_1 *sig_2^2) ) )
  l_34.beta_2_rho<-  t(X_4) %*% ( ( 2*rho/(1-rho^2)^2 * ( (Y_4-mu_4)/sig_2^2 - rho * (Y_3-mu_3)/(sig_1 *sig_2) ) ) - 1/(1-rho^2) * (Y_3-mu_3)/(sig_1 *sig_2))


  l_34.sig_1sig_2 <- sum ( 1/(1-rho^2) * rho *(Y_3 - mu_3)*(Y_4 - mu_4)/(sig_1^2 * sig_2^2) )


  l_34.sig_1_rho <- sum(  2 *rho /( (1-rho^2)^2 * sig_1) * ( (Y_3-mu_3)^2/sig_1^2 - rho * (Y_3-mu_3) * (Y_4-mu_4)/(sig_1 *sig_2) )-
                            1 /( (1-rho^2) * sig_1) *(Y_3 - mu_3)*(Y_4 - mu_4)/(sig_1 * sig_2)  )


  l_34.sig_2_rho <- sum(  2 *rho /( (1-rho^2)^2 * sig_2) * ( (Y_4-mu_4)^2/sig_2^2 - rho * (Y_3-mu_3) * (Y_4-mu_4)/(sig_1 *sig_2) )-
                            1 /( (1-rho^2) * sig_2) *(Y_3 - mu_3)*(Y_4 - mu_4)/(sig_1 * sig_2)  )


  hessian<-matrix(0, ( N_3 + N_4 +3) ,(N_3 + N_4 +3) )


  hessian[1:N_3,1:N_3] <-l_34.beta_1beta_1
  hessian[(N_3 + 1) : (N_3 + N_4) , (N_3 + 1) : (N_3 + N_4)  ] <-l_34.beta_2beta_2


  hessian[(N_3 + N_4 +1) , (N_3 + N_4 +1)] <-l_34.sig_1sig_1
  hessian[(N_3 + N_4 +2) , (N_3 + N_4 +2)] <-l_34.sig_2sig_2
  hessian[(N_3 + N_4 +3), (N_3 + N_4 +3)] <-l_34.rho_rho

  hessian[1:N_3, (N_3 +1):(N_3 + N_4)] <-l_34.beta_1beta_2
  hessian[(N_3 +1):(N_3 + N_4) , 1:N_3] <-t(l_34.beta_1beta_2)

  hessian[1:N_3,(N_3 + N_4 +1)] <- l_34.beta_1sig_1
  hessian[(N_3 + N_4 +1),1:N_3] <- t(l_34.beta_1sig_1 )
  hessian[1:N_3,(N_3 + N_4 +2)] <- l_34.beta_1sig_2
  hessian[(N_3 + N_4 +2),1:N_3] <- t(l_34.beta_1sig_2)
  hessian[1:N_3,(N_3 + N_4 +3)] <- l_34.beta_1_rho
  hessian[(N_3 + N_4 +3),1:N_3] <- t(l_34.beta_1_rho)

  hessian[(N_3 +1):(N_3 + N_4) ,(N_3 + N_4 +1)] <- l_34.beta_2sig_1
  hessian[(N_3 + N_4 +1),(N_3 +1):(N_3 + N_4 )] <- t(l_34.beta_2sig_1 )

  hessian[(N_3 +1):(N_3 + N_4) ,(N_3 + N_4 +2)] <- l_34.beta_2sig_2
  hessian[(N_3 + N_4 +2),(N_3 +1):(N_3 + N_4) ] <- t(l_34.beta_2sig_2)
  hessian[(N_3 +1):(N_3 + N_4) ,(N_3 + N_4 +3)] <- l_34.beta_2_rho
  hessian[(N_3 + N_4 +3),(N_3 +1):(N_3 + N_4) ] <- t(l_34.beta_2_rho)

  hessian[(N_3 + N_4 +1),(N_3 + N_4 +2)] <- l_34.sig_1sig_2
  hessian[(N_3 + N_4 +2),(N_3 + N_4 +1)] <- l_34.sig_1sig_2
  hessian[(N_3 + N_4 +1),(N_3 + N_4 +3)] <- l_34.sig_1_rho
  hessian[(N_3 + N_4 +3),(N_3 + N_4 +1)] <- l_34.sig_1_rho
  hessian[(N_3 + N_4 +2),(N_3 + N_4 +3)] <- l_34.sig_2_rho
  hessian[(N_3 + N_4 +3),(N_3 + N_4 +2)] <- l_34.sig_2_rho


  grad<- c(t(l_34.beta3), t(l_34.beta4), l_34.sig_1, l_34.sig_2 , l_34.rho)
  grad_point <- cbind(l_34.beta3_point, l_34.beta4_point, l_34.sig_1_point, l_34.sig_2_point , l_34.rho_point)

  l_34<-0
  attr(l_34,"gradient")<- grad
  attr(l_34,"hessian")<- hessian
  attr(l_34,"grad_point")<-grad_point

  return(l_34)
}


# Function Case_1_discrete: Pair of one binary and one continuous outcome
Case_1_discrete <- function( data_1,data_2, beta_1, beta_2, sig_1,sig_2, rho_12) {

  rho <- rho_12
  N<-nrow(data_1)
  N_j <-length(beta_1)
  N_k <- length(beta_2)
  Y_j<- data_1[,ncol(data_1)]
  Y_k<- data_2[,ncol(data_2)]
  X_j<- as.matrix (data_1[,-ncol(data_1)])
  X_k<- as.matrix (data_2[,-ncol(data_2)])


  s_k<-(Y_k- X_k %*% beta_2 )/sig_2
  o_j <- -( X_j %*% beta_1 + rho * s_k )/(sqrt(1-rho^2 ))
     
  phi_oj <- dnorm (o_j)
  
  pnormo_j<-pnorm(o_j)
  #pnormo_j<-ifelse((1-pnorm(o_j))**2<0.01,0.9,pnormo_j)
  #pnormo_j<-ifelse(pnorm(o_j)**2<0.01,0.1,pnormo_j)
  Phi_oj <- rep(0, length(o_j))
  Phi_oj_2 <- rep(0, length(o_j))
  for ( val in 1:length(o_j)  )   {
    Phi_oj[val] <- if (Y_j[val] == 0 )  1/pnormo_j[val] else - 1/ ( 1- pnormo_j[val])
    Phi_oj_2[val] <- if (Y_j[val] == 0 )  - 1/(pnormo_j[val] **2) else  - 1/(( 1- pnormo_j[val])**2 )
  }
    if (any(is.infinite(Phi_oj))==T | any(is.infinite(Phi_oj_2))==T) return()
 #Phi_oj[Phi_oj==Inf]= 1666667
  #Phi_oj[Phi_oj==-Inf]=-1666667
  #Phi_oj_2[Phi_oj_2==Inf]= 1666667
  #Phi_oj_2[Phi_oj_2==-Inf]=-1666667

  d_s_beta_2 <- - 1/sig_2 * X_k
  d_s_sig_2 <- - 1/sig_2 * s_k

  d_o_beta_2 <-  X_k * rho/(sig_2 * sqrt (1- rho^2))
  d_o_beta_1 <-  - X_j / (  sqrt (1- rho^2))
  d_o_sig_2 <-   rho * s_k / (sig_2 * sqrt (1- rho^2))
  d_o_d_rho <-  -s_k/sqrt (1- rho^2) - X_j %*% beta_1 * rho / (  (1- rho^2)^(3/2) ) - rho^2 * s_k /(1- rho^2)^(3/2)

  d2_s_sig_2_2 <-  2*s_k /sig_2^2
  d2_s_beta_2_sig_2 <-  1/sig_2^2 * X_k

  d2_o_rho_sig_2 <- s_k/sig_2 * (  1/ sqrt(1-rho^2 ) + rho^2 * (1- rho^2)^(-3/2)   )
  d2_o_beta_1_rho <- -  X_j * rho /  (   (1- rho^2)^(3/2) )
  d2_o_beta_2_rho <-  X_k /sig_2 * (  1/ sqrt(1-rho^2 ) + rho^2 * (1- rho^2)^(-3/2) )
  d2_o_sig_2_2 <- - 2 * rho * s_k  / ( sig_2^2 * sqrt(1-rho^2 )  )
  d2_o_beta_2_sig_2 <- - X_k * rho / ( sig_2^2 * sqrt(1-rho^2 )  )
  d2_o_rho_2 <- - rho * s_k /((1-rho^2)^(3/2)) -  X_j %*% beta_1 *  (  1/ ((1-rho^2)^(3/2)) + 3* rho ^2 / ((1-rho^2)^(5/2)) ) -
    s_k * (2*rho/ ((1-rho^2)^(3/2)) + 3* rho ^3 / ((1-rho^2)^(5/2)) )


  d_sig_2 <- matrix( 0,nrow=length(o_j),ncol= (N_k+N_j+2))
  d_s <- matrix( 0,nrow=length(o_j),ncol= (N_k+N_j+2))
  d_o <- matrix( 0,nrow=length(o_j),ncol= (N_k+N_j+2))

  d_sig_2[,(N_k+N_j+1)]<- 1
  d_s[,1:N_k] <- d_s_beta_2
  d_s[,(N_k+N_j+1)] <- d_s_sig_2
  d_o[,1:N_k] <- d_o_beta_2
  d_o[,(N_k + 1):( N_k + N_j )] <- d_o_beta_1
  d_o[,(N_k+N_j+1)] <- d_o_sig_2
  d_o[,(N_k+N_j+2)] <- d_o_d_rho



  s_k_M <- rep( s_k, (N_k+N_j+2))
  Phi_oj_M <- rep (Phi_oj, (N_k+N_j+2))
  phi_oj_M <- rep (phi_oj, (N_k+N_j+2))

  d_l_d_eta_N <- -1/sig_2 * d_sig_2 - s_k_M * d_s + Phi_oj_M * phi_oj_M * d_o
  d_l_d_eta<-colSums( d_l_d_eta_N  )
  d_l_d_eta


  # Hessian

  d2_l_eta_etat <- matrix(0,(N_k+N_j+2),(N_k+N_j+2))

  d2_s_M <- array(0, dim=c( (N_k+N_j+2), (N_k+N_j+2), N))
  for ( ids in 1:N)
  {
    d2_s_M[(N_k+N_j+1),(N_k+N_j+1), ids] <-  d2_s_sig_2_2[ids]
    d2_s_M[1:N_k,(N_k+N_j+1), ids ] <- d2_s_beta_2_sig_2[ids,]
    d2_s_M[(N_k+N_j+1),1:N_k, ids] <- d2_s_beta_2_sig_2[ids,]
  }
  d2_o_M <- array(0, dim=c((N_k+N_j+2), (N_k+N_j+2), N ))
  for ( ids in 1:N)
  {
    d2_o_M[(N_k+N_j+2),(N_k+N_j+2), ids ] <- d2_o_rho_2[ids]
    d2_o_M[(N_k+N_j+1),(N_k+N_j+1), ids ] <-  d2_o_sig_2_2[ids]
    d2_o_M[1:N_k,(N_k+N_j+1), ids ] <- d2_o_beta_2_sig_2[ids,]
    d2_o_M[(N_k+N_j+1),1:N_k, ids ] <- d2_o_beta_2_sig_2[ids,]
    d2_o_M[1:N_k,(N_k+N_j+2),ids ] <- d2_o_beta_2_rho[ids,]
    d2_o_M[(N_k+N_j+2),1:N_k,ids]<- d2_o_beta_2_rho[ids,]
    d2_o_M[(N_k+1):(N_k+N_j),(N_k+N_j+2),ids]<- d2_o_beta_1_rho[ids,]
    d2_o_M[(N_k+N_j+2),(N_k+1):(N_k+N_j),ids]<- d2_o_beta_1_rho[ids,]

    d2_o_M[(N_k+N_j+2),(N_k+N_j+1),ids]<- d2_o_rho_sig_2[ids]
    d2_o_M[(N_k+N_j+1),(N_k+N_j+2),ids]<- d2_o_rho_sig_2[ids]
  }






  for ( rows in c(1: (N_k+N_j+2)  ))   {
    for ( cols in c(1:(N_k+N_j+2))) {
      d2_l_eta_etat [rows, cols] <- sum(  1/sig_2^2 * d_sig_2[, rows] * d_sig_2[,cols] - d_s[,rows] * d_s[,cols] - s_k * d2_s_M[rows,cols, ]
                                          +  Phi_oj_2 * phi_oj^2  * d_o[,rows] * d_o[,cols]
                                          +  Phi_oj * ( - o_j * phi_oj * d_o[,rows] * d_o[,cols] + phi_oj * d2_o_M [rows,cols, ] )
      )


    }
  }

  l_jk <-0

  attr(l_jk,"gradient")<- d_l_d_eta
  attr(l_jk,"hessian")<- d2_l_eta_etat
  attr(l_jk,"grad_point")<-  d_l_d_eta_N
  l_jk


  return(l_jk)
}
