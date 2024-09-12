cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499",
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
# set the commonly used color
blue <- cbPalette[6]
red <- cbPalette[7]
green <- cbPalette[4]
orange <- cbPalette[2]
grey <- cbPalette[1]
yellow <- cbPalette[5]
purple <- cbPalette[8]
skyblue <- cbPalette[3]


##################################################
# Expression-wide regression functions
##################################################


#########################################################
# model 2: xi_1 + xi_2 -> X, xi -> Y, xi_1 + xi_2 -> Z_j
#########################################################
###################################################
# construct model with two latent variables
###################################################

# M-step
two_factor_EIV_Mstep <- function(para){
  para$beta <- ( t(para$E_xi[,1]) %*% para$Y_0asNA ) / ( t(para$E2_xi[,1]) %*% para$Y_NAidx )
  para$sum_xi_xiT <- t(para$Z_NAidx) %*% para$E2_xi
  Zn_xiT <- t( (t(para$E_xi) %*% para$Z_0asNA) )
  para$eta_j <- sapply(1:para$J, function(j) Zn_xiT[j,,drop =F] %*% solve(matrix(para$sum_xi_xiT[j,], nrow = 2, ncol = 2)) ) %>% t()

  # estimate the variance across variables (Phi)
  para$sigma_X <- mean(rowSums(para$E2_xi) - 2*para$X*(rowSums(para$E_xi)) + para$X^2, na.rm = T)
  para$sigma_Y <- mean(para$E2_xi[,1]*para$beta[1]^2 - 2*para$Y*para$E_xi[,1]*para$beta[1] + para$Y^2, na.rm = T)
  para$Phi_hat_diag <- ( c( t( colSums(para$Z_0asNA^2)) ) -
                           2* rowSums(para$eta_j * (t(para$Z_0asNA) %*% para$E_xi))  +
                           sapply(1:para$J, function(j) t(para$eta_j[j,]) %*% matrix(para$sum_xi_xiT[j,], nrow = 2, ncol = 2) %*% para$eta_j[j,] ) ) /para$N_zj

  para$sigma_xi_1 <- mean(para$E2_xi[,1], na.rm = T)
  para$sigma_xi_2 <- mean(para$E2_xi[,4], na.rm = T)
  return(para)
}
# E-step:
two_factor_EIV_Estep <- function(para){
  # note: to handle missing data, we use only those non-missing variables in the likelihood function; therefore we have different M_n for each n
  HT_H <- cbind(para$eta_j[,1]^2, para$eta_j[,1]*para$eta_j[,2], para$eta_j[,2]*para$eta_j[,1], para$eta_j[,2]^2) # note here we used 4-element vector to represent the 2*2 matrix
  HT_phi_H <- para$Z_NAidx %*% (diag(1/para$Phi_hat_diag) %*% HT_H) +
    1/para$sigma_X * ( para$X_NAidx %*% t(c(1,1,1,1))) + # X and Y need to be added
    1/para$sigma_Y * ( para$Y_NAidx %*% t(c(para$beta^2,0,0,0)))

  para$M_n <- matrix(rep(c(1/para$sigma_xi_1, 0,0, 1/para$sigma_xi_2),dim(HT_phi_H)[1]), byrow = T, nrow = dim(HT_phi_H)[1])  +
    HT_phi_H
  HT_phi_Xn <- (cbind(c(1,1), c(para$beta, 0), t(para$eta_j)) %*%
                  diag(c(1/para$sigma_X, 1/para$sigma_Y, 1/para$Phi_hat_diag)) ) %*%
    t( cbind(para$X_0asNA, para$Y_0asNA, para$Z_0asNA) )

  # compute the inverse of M_n, again expressed 2x2 as 4x1 for computational efficiency
  M_n_inv <- sapply(1:dim(HT_phi_Xn)[2], function(j)
    c(solve(matrix(para$M_n[j,], nrow = 2))) ) %>%
    t()

  para$E_xi <- sapply(1:dim(HT_phi_Xn)[2], function(j)
    matrix(M_n_inv[j,], nrow = 2) %*%  HT_phi_Xn[,j]) %>%
    t()

  para$E2_xi <- cbind(para$E_xi[,1]^2, para$E_xi[,1]*para$E_xi[,2], para$E_xi[,2]*para$E_xi[,1], para$E_xi[,2]^2) +
    M_n_inv
  return(para)
}

# lower bound
two_factor_EIV_lb <- function(para){
  # only for Zs compute HT_Phi_T
  HT_H <- cbind(para$eta_j[,1]^2, para$eta_j[,1]*para$eta_j[,2], para$eta_j[,2]*para$eta_j[,1], para$eta_j[,2]^2) # note here we used 4-element vector to represent the 2*2 matrix
  etaT_phi_eta <- para$Z_NAidx %*% (diag(1/para$Phi_hat_diag) %*% HT_H)

  ELBO_factor_EIV_univariate <- -1/2 * log(para$sigma_X ) * para$N_X -1/2 * log(para$sigma_Y) * para$N_Y - 1/2 * (log(para$sigma_xi_1) + log(para$sigma_xi_2))* length(para$X) -1/2 * (log(para$Phi_hat_diag) %*% para$N_zj) -
    1/2 * ( (sum(para$E2_xi[,1], na.rm = T))/para$sigma_xi_1 + (sum(para$E2_xi[,4], na.rm = T))/para$sigma_xi_2) -
    1/2* sum( (rowSums(para$E2_xi) - 2*para$X*rowSums(para$E_xi) + para$X^2)/para$sigma_X, na.rm = T) -
    1/2* sum( (para$E2_xi[,1] * as.numeric(para$beta)^2 - 2*para$Y* para$E_xi[,1] * as.numeric(para$beta) + para$Y^2)/para$sigma_Y, na.rm = T) -
    # Vector format of the quadratic term
    1/2* ( sum(para$E2_xi * etaT_phi_eta)-
             2 * sum( ( ( para$E_xi %*% t(para$eta_j) ) * para$Z_0asNA ) %*% (1/para$Phi_hat_diag) ) +
             sum(para$Z_0asNA^2 %*% (1/para$Phi_hat_diag)) ) +
    # entropy term (ignore constant): note this is a positive as the entropy is defined with a negative term
    1/2* sum( sapply(1:dim(etaT_phi_eta)[1], function(j)
      log(abs(det(matrix(para$M_n[j,], nrow = 2)) ) ) ) )
  return(ELBO_factor_EIV_univariate)
}

wrapper_two_factor_EIV <- function(X, Y, Z, default_update_max = 2000, converge_ratio = 10^(-6)){
  para <- list()

  para$X <- scale(X)
  para$Y <- scale(Y)
  para$Z <- apply(Z, 2, scale)
  para$max_itr <- default_update_max
  para$converge_thre <- converge_ratio
  para$lb <- data.frame("Iteration" = as.numeric(),"Lower_bound" = as.numeric())

  para$X_0asNA <- para$X
  para$X_0asNA[is.na(para$X_0asNA)] <- 0
  para$X_NAidx <- 1-is.na(para$X) # if NA use 0, otherwise 1
  para$Y_0asNA <- para$Y
  para$Y_0asNA[is.na(para$Y_0asNA)] <- 0
  para$Y_NAidx <- 1-is.na(para$Y) # if NA use 0, otherwise 1
  para$Z_0asNA <- para$Z
  para$Z_0asNA[is.na(para$Z_0asNA)] <- 0
  para$Z_NAidx <- 1-is.na(para$Z)

  para$J <- dim(para$Z)[2] # number of helper proteins

  para$N_zj <- sapply(1:para$J, function(j) sum(!is.na(para$Z[,j]))) # how many measurement for each helper protein?
  para$N_X <- sum(!is.na(para$X))
  para$N_Y <- sum(!is.na(para$Y))

  # initialise expected value and M-step to speed up the estimation
  # para$E_xi <-  cbind(rnorm(length(para$X_0asNA)) + para$X_0asNA, rnorm(length(para$X_0asNA)) + para$X_0asNA)
  Y_X_lm <- lm(para$X_0asNA ~ para$Y_0asNA)
  para$E_xi <-  cbind(c(Y_X_lm$fitted.values) + 0.01 * rnorm(length(para$X_0asNA)), c(Y_X_lm$residuals) + 0.01 * rnorm(length(para$X_0asNA)))
  para$E2_xi <- cbind(para$E_xi[,1]^2,  para$E_xi[,1]*para$E_xi[,2], para$E_xi[,2]*para$E_xi[,1], para$E_xi[,2]^2) + 0.1  # note: for computational efficiency, we will use a single 4 element vector to represent the 2*2 covariance matrix

  # updates
  for(itr in 1:para$max_itr){
    para <- two_factor_EIV_Mstep(para)
    para <- two_factor_EIV_Estep(para)
    # compute ELBO
    ELBO_factor_EIV_univariate <- two_factor_EIV_lb(para)
    print(paste0("Current Lower bound ", ELBO_factor_EIV_univariate, " at iteration: ",itr))
    para$lb[nrow(para$lb) + 1,] <- c(itr, ELBO_factor_EIV_univariate)
    if(itr > 1){
      curr_lb <- pull(filter(para$lb, Iteration == itr), Lower_bound)
      prev_lb <- pull(filter(para$lb, Iteration == (itr -1 )), Lower_bound)

      try({
        if(is.finite((curr_lb - prev_lb)) & (curr_lb - prev_lb)/abs(prev_lb) < para$converge_thre ){
          print(paste0("Optimization converged at step ", itr))
          break
        }
      })
    }
  }
  return(para)
}

bioProp_EWR_factor_adjustment <- function(X, Y, Z, factor_model_update_max = 2000){ # X is Nx1, Y is Nx1, Z is NxJ
  J <- dim(Z)[2]

  PxPy <- summary(lm(scale(X) ~ scale(Y)))$coefficients[2,1]^2
  # beta_EWR
  OLS_x2z <- c()
  OLS_y2z <- c()
  for(j in 1:J){
    x2z_coef <- summary(lm(scale(Z[,j]) ~ scale(X)))$coefficients[2,1]
    y2z_coef <- summary(lm(scale(Z[,j]) ~ scale(Y)))$coefficients[2,1]

    OLS_x2z <- c(OLS_x2z, x2z_coef)
    OLS_y2z <- c(OLS_y2z, y2z_coef)
  }
  beta_EWR <- summary(lm(OLS_y2z ~ OLS_x2z))$coefficients[2,1]

  approx_Bioprop_xi1_xi2 <- sqrt(PxPy/beta_EWR^2)

  suppressMessages({
    para_simulation <- wrapper_two_factor_EIV(X = X, Y= Y, Z = Z, default_update_max = factor_model_update_max)
  })

  print(paste0("Factor model converged at update ", dim(para_simulation$lb)[1]))
  # compute two ratios separately (use (eta_1^2 * xi1^2) / (eta_1^2 * xi1^2) as it is invariant of scaling of eta_1 and xi_1)
  coef_1 <- summary(lm(para_simulation$X_0asNA ~ para_simulation$E_xi[,1]))$coefficients[2,1]
  coef_2 <- summary(lm(para_simulation$X_0asNA ~ para_simulation$E_xi[,2]))$coefficients[2,1]
  ratio_xi1_xi2 <- (para_simulation$sigma_xi_1 * coef_1^2)/(para_simulation$sigma_xi_2 * coef_2^2)
  ratio_eta_j_xi <- (var(para_simulation$eta_j[,1])/var(para_simulation$eta_j[,2])) *
    (para_simulation$sigma_xi_1/para_simulation$sigma_xi_2)
  # use the ratio to recover biological proportion
  Bioprop_xi1_xi2 <- (1+ratio_xi1_xi2) * ratio_eta_j_xi / (1+ratio_xi1_xi2* ratio_eta_j_xi) * approx_Bioprop_xi1_xi2

  results_biprop <- list()
  results_biprop$finish_itr <- dim(para_simulation$lb)[1]
  results_biprop$Bioprop_xi1_xi2 <- Bioprop_xi1_xi2
  results_biprop$PxPy <- PxPy
  results_biprop$beta_EWR <- beta_EWR
  results_biprop$approx_Bioprop_xi1_xi2 <- approx_Bioprop_xi1_xi2

  return(results_biprop)
}

######################################################
# 20240711: five different methods of latent factor analysis
######################################################
padding_X_Y_cov <- function(X, Y, Z, X_padding_num, Y_padding_num){
  X <- X %>%
    scale()
  Y <- Y %>%
    scale()
  Z <- apply(Z, 2, scale)
  ZtZ_cov <- cov(Z, use = "pairwise.complete.obs")

  # first extend the covariance matrix for paddings of Y
  Y_Z_cov <- cov(x = Z, y = Y, use = "pairwise.complete.obs")
  a <- matrix(var(Y, na.rm = T), nrow = Y_padding_num, ncol = Y_padding_num)
  b <- matrix(rep(Y_Z_cov, Y_padding_num) , nrow = length(Y_Z_cov), ncol = Y_padding_num)
  c <- ZtZ_cov
  # combine covariance abc
  ab <- rbind(a, b)
  bc <- rbind(t(b), c)
  Y_Z_padded <- cbind(ab, bc)

  # second extend the paddings to X
  X_Z_cov <- cov(x = Z, y = X, use = "pairwise.complete.obs")
  a <- matrix(var(X, na.rm = T), nrow = X_padding_num, ncol = X_padding_num)
  # b is exteded with covariance between Y and X
  b <- matrix(rep(c(rep(cov(x = X, y = Y, use = "pairwise.complete.obs"), Y_padding_num),X_Z_cov), X_padding_num),
              nrow = length(X_Z_cov) + Y_padding_num, ncol = X_padding_num)
  c <- Y_Z_padded
  ab <- rbind(a, b)
  bc <- rbind(t(b), c)
  XY_Z_padded <- cbind(ab, bc)
  return(XY_Z_padded)
}
# upweight X and Y
padding_XY_Z_cov <- function(X, Y, Z, num_padding_XY){
  XY_Z_padded <- padding_X_Y_cov(X, Y, Z, X_padding_num=num_padding_XY, Y_padding_num=num_padding_XY)
  return(XY_Z_padded)
}


######## 4. ECME padded factor analysis; assuming latent factor having variance = 1
# E-step:for the covariance only estimation, assuming missing at random.
two_factor_padding_COV_ECME_Estep <- function(para){

  ############ below needs to be changed -- M need to consider the padded matrix.
  Phi_inv_H <- matrix(rep(1/para$Phi_hat_diag_padded, dim(para$H_padded)[2]), ncol = dim(para$H_padded)[2]) * para$H_padded
  HT_phi_H <-  t(para$H_padded) %*% Phi_inv_H

  para$M_n <- diag(rep(1, dim(para$H_padded)[2])) + HT_phi_H
  para$delta <- solve(para$M_n) %*% t(Phi_inv_H)

  # for E step, no need to consider the padding as padding doesn't affect E and M-step of H.
  para$CX_xi <- para$XY_Z_padded %*% t(para$delta)
  para$Cxi_xi <- para$delta %*% para$XY_Z_padded %*% t(para$delta) + solve(para$M_n)
  return(para)
}
two_factor_padding_COV_ECME_Mstep_H <- function(para){
  para$alpha <- para$CX_xi[1,,drop=F] %*% solve(para$Cxi_xi)
  para$beta <- para$CX_xi[para$XY_padding_num + 1, 1] / para$Cxi_xi[1,1]
  para$eta_j <- para$CX_xi[(para$XY_padding_num*2 + 1):(dim(para$H_padded)[1]),,drop=F] %*% solve(para$Cxi_xi)

  # assigning the number to the original matrix
  para$H <- rbind(para$alpha, c(para$beta, 0), para$eta_j)
  para$H_padded <- rbind(matrix(rep(para$alpha, each = para$XY_padding_num), ncol = 2),
                         matrix(rep(c(para$beta, 0), each = para$XY_padding_num), ncol = 2), para$eta_j)
  return(para)
}

Marginal_likelihood_cov_ECME <- function(x, para){
  phi_padded <- c(rep(x[1], para$XY_padding_num),
                  rep(x[2], para$XY_padding_num),
                  x[3:length(x)])

  #####################################
  ### fast inverse matrix of woodbury identity
  dim_H <- dim(para$H_padded)[2]
  Phi_inv_H <- matrix(rep(1/phi_padded , dim_H), ncol = dim_H) * para$H_padded
  inverse_woodbury_fast <- diag(1/phi_padded) - Phi_inv_H %*% solve(diag(rep(1,dim_H)) + t(para$H_padded) %*% Phi_inv_H) %*% t(Phi_inv_H)

  # max(abs(inverse_woodbury_fast - solve((HHT + phi_matrix)) ))
  #####################################
  ### fast inverse matrix of determinant
  logdet_fast <- log(det(diag(rep(1,dim_H)) + t(para$H_padded) %*% Phi_inv_H)) +  sum(log(phi_padded))
  # logdet <- log(det(diag(phi_padded) + para$H_padded %*% t(para$H_padded)))
  # ingoring a factor of N/2
  return( -logdet_fast-sum(inverse_woodbury_fast * para$XY_Z_padded) )
}

gradient_likelihood_cov_ECME <-  function(x, para){
  phi_padded <- c(rep(x[1], para$XY_padding_num),
                  rep(x[2], para$XY_padding_num),
                  x[3:length(x)])

  #####################################
  ### fast inverse matrix of woodbury identity
  dim_H <- dim(para$H_padded)[2]
  Phi_inv_H <- matrix(rep(1/phi_padded , dim_H), ncol = dim_H) * para$H_padded
  inverse_woodbury_fast <- diag(1/phi_padded) - Phi_inv_H %*% solve(diag(rep(1,dim_H)) + t(para$H_padded) %*% Phi_inv_H) %*% t(Phi_inv_H)

  deriv_L <- diag(-inverse_woodbury_fast) + rowSums(inverse_woodbury_fast * t(para$XY_Z_padded %*% inverse_woodbury_fast))
  # diag(-inverse_woodbury_fast + inverse_woodbury_fast %*% para$XY_Z_padded %*% inverse_woodbury_fast)
  # ingnoring a factor of N/2
  gradient <- c(sum(deriv_L[1:para$XY_padding_num]),
                sum(deriv_L[(1 + para$XY_padding_num):(2*para$XY_padding_num)]),
                deriv_L[(2*para$XY_padding_num+1):length(deriv_L)])
  return( gradient )
}


two_factor_padding_COV_ECME_CMstep_uniqueness <- function(para){
  estimate_uniqueness <- optim(par = c(para$sigma_X, para$sigma_Y, para$Phi_hat_diag),
                               fn = function(x) -Marginal_likelihood_cov_ECME(x, para),
                               gr = function(x) -gradient_likelihood_cov_ECME(x, para),
                               method = "L-BFGS-B",
                               lower = 10^-6, upper = 1-10^-6)$par # optim is bounded for the uniqueness (variance term)
  para$sigma_X <- estimate_uniqueness[1]
  para$sigma_Y <- estimate_uniqueness[2]
  para$Phi_hat_diag <- estimate_uniqueness[3:length(estimate_uniqueness)]
  para$Phi_hat_diag_padded <- c(rep(para$sigma_X, para$XY_padding_num),
                                rep(para$sigma_Y, para$XY_padding_num),
                                para$Phi_hat_diag)
  return(para)
}

# evidence lower bound
two_factor_padding_COV_ECME_lb <- function(para){
  dim_H <- dim(para$H_padded)[2]
  Phi_inv_H <- matrix(rep(1/para$Phi_hat_diag_padded, dim_H), ncol = dim_H) * para$H_padded
  # ignoring a factor of N/2
  term1 <- -sum(log(para$Phi_hat_diag_padded))
  term2 <- -sum(diag(para$Cxi_xi))
  term3 <- -sum(diag(para$XY_Z_padded) * 1/para$Phi_hat_diag_padded) +
    2*sum(diag(t(Phi_inv_H) %*% para$CX_xi)) -
    sum(diag(t(Phi_inv_H) %*% para$H_padded %*% para$Cxi_xi))
  term4 <- log(det(para$M_n))
  return(term1 + term2 + term3 + term4)
}

wrapper_two_factor_padding_COV_ECME <- function(X, Y, Z, XY_padding_num = 1, default_update_max = 2000, converge_ratio = 10^(-6)){
  X <- scale(X)
  Y <- scale(Y)
  Z <- apply(Z, 2, scale)

  XY_Z_padded <- padding_XY_Z_cov(X, Y, Z, num_padding_XY = XY_padding_num)
  # get the covariance matrix without padding
  XY_Z <- padding_XY_Z_cov(X, Y, Z, num_padding_XY = 1)

  para <- list()

  para$XY_Z_padded <- XY_Z_padded
  para$XY_Z <- XY_Z
  para$XY_padding_num <- XY_padding_num
  para$Z <- apply(Z, 2, scale)
  para$max_itr <- default_update_max
  para$converge_thre <- converge_ratio
  para$lb <- data.frame("Iteration" = as.numeric(),"Log_likelihood" = as.numeric(), "ELBO" = as.numeric())

  # initialize the parameters
  para$alpha <- c(1,1) # initial alpha as equal to xi1 and xi2
  para$beta <- para$XY_Z_padded[para$XY_padding_num + 1,1]
  para$eta_j <- para$XY_Z_padded[(2*para$XY_padding_num + 1):dim(XY_Z_padded)[1],c(1:2)]
  para$H <- rbind(para$alpha, c(para$beta, 0), para$eta_j)
  para$H_padded <- rbind(matrix(rep(para$alpha, each = para$XY_padding_num), ncol = 2),
                         matrix(rep(c(para$beta, 0), each = para$XY_padding_num), ncol = 2), para$eta_j)
  para$sigma_X <- 1- para$beta^2
  para$sigma_Y <- 1- para$beta^2
  para$Phi_hat_diag <- 1 - para$eta_j[,1]^2
  para$Phi_hat_diag_padded <- c(rep(para$sigma_X, para$XY_padding_num),
                                rep(para$sigma_Y, para$XY_padding_num),
                                para$Phi_hat_diag)

  # updates
  for(itr in 1:para$max_itr){
    para <- two_factor_padding_COV_ECME_Estep(para)
    para <- two_factor_padding_COV_ECME_Mstep_H(para)
    para <- two_factor_padding_COV_ECME_CMstep_uniqueness(para)
    # compute ELBO
    marginal_ll_ECME <- Marginal_likelihood_cov_ECME(c(para$sigma_X, para$sigma_Y, para$Phi_hat_diag), para)
    ELBO_ECME <- two_factor_padding_COV_ECME_lb(para)
    print(paste0("Current log likelihood ", marginal_ll_ECME, " at iteration: ",itr))
    para$lb[nrow(para$lb) + 1,] <- c(itr, marginal_ll_ECME, ELBO_ECME)
    if(itr > 1){
      curr_lb <- pull(filter(para$lb, Iteration == itr), Log_likelihood)
      prev_lb <- pull(filter(para$lb, Iteration == (itr -1 )), Log_likelihood)

      try({
        if(is.finite((curr_lb - prev_lb)) & (curr_lb - prev_lb)/abs(prev_lb) < para$converge_thre ){
          print(paste0("Optimization converged at step ", itr))
          break
        }
      })
    }
  }
  return(para)
}

######  6. factor model with multiple factors: this allows more factors in the factor model; extended from a simple two factor model
Marginal_likelihood_cov_multifactor <- function(x, para){
  phi_padded <- c(rep(x[1], para$X_padding_num),
                  rep(x[2], para$Y_padding_num),
                  x[3:(length(para$Phi_hat_diag) + 2)])

  H_padded <- rbind(matrix(rep(c(x[(length(para$Phi_hat_diag) + 3):(length(para$Phi_hat_diag) + para$XZ_factor + 3)], rep(0,para$Z_factor)), each = para$X_padding_num), nrow = para$X_padding_num),
                    matrix(rep(c(x[(length(para$Phi_hat_diag) + length(para$alpha) + 3)], rep(0, para$XZ_factor + para$Z_factor) ), each = para$Y_padding_num), nrow = para$Y_padding_num),
                    matrix(x[(length(para$Phi_hat_diag) + length(para$alpha) + 4): length(x)], ncol = 1+para$XZ_factor + para$Z_factor))
  # identical(H_padded, para$H_padded)

  #####################################
  ### fast inverse matrix of woodbury identity
  dim_H <- dim(H_padded)[2]
  Phi_inv_H <- matrix(rep(1/phi_padded , dim_H), ncol = dim_H) * H_padded
  inverse_woodbury_fast <- diag(1/phi_padded) - Phi_inv_H %*% solve(diag(rep(1,dim_H)) + t(H_padded) %*% Phi_inv_H) %*% t(Phi_inv_H)

  # max(abs(inverse_woodbury_fast - solve((HHT + phi_matrix)) ))
  #####################################
  ### fast inverse matrix of determinant
  logdet_fast <- log(det(diag(rep(1,dim_H)) + t(H_padded) %*% Phi_inv_H)) +  sum(log(phi_padded))
  # logdet <- log(det(diag(phi_padded) + H_padded %*% t(H_padded)))
  # ingoring a factor of N/2
  return( -logdet_fast-sum(inverse_woodbury_fast * para$XY_Z_padded) )
}

gradient_likelihood_cov_multifactor <-  function(x, para){
  phi_padded <- c(rep(x[1], para$X_padding_num),
                  rep(x[2], para$Y_padding_num),
                  x[3:(length(para$Phi_hat_diag) + 2)])

  H_padded <- rbind(matrix(rep(c(x[(length(para$Phi_hat_diag) + 3):(length(para$Phi_hat_diag) + para$XZ_factor + 3)], rep(0,para$Z_factor)), each = para$X_padding_num), nrow = para$X_padding_num),
                    matrix(rep(c(x[(length(para$Phi_hat_diag) + length(para$alpha) + 3)], rep(0, para$XZ_factor + para$Z_factor) ), each = para$Y_padding_num), nrow = para$Y_padding_num),
                    matrix(x[(length(para$Phi_hat_diag) + length(para$alpha) + 4): length(x)], ncol = 1+para$XZ_factor + para$Z_factor))
  #####################################
  ### fast inverse matrix of woodbury identity
  dim_H <- dim(H_padded)[2]
  Phi_inv_H <- matrix(rep(1/phi_padded , dim_H), ncol = dim_H) * H_padded
  inverse_woodbury_fast <- diag(1/phi_padded) - Phi_inv_H %*% solve(diag(rep(1,dim_H)) + t(H_padded) %*% Phi_inv_H) %*% t(Phi_inv_H)
  gradient_matrix <- -inverse_woodbury_fast + inverse_woodbury_fast %*%  para$XY_Z_padded %*% inverse_woodbury_fast
  deriv_Phi <- diag(gradient_matrix)
  gradient_Phi <- c(sum(deriv_Phi[1:para$X_padding_num]),
                    sum(deriv_Phi[(1 + para$X_padding_num):(para$X_padding_num + para$Y_padding_num)]),
                    deriv_Phi[(para$X_padding_num + para$Y_padding_num+1):length(deriv_Phi)])
  deriv_H <- 2* gradient_matrix %*% H_padded
  gradient_H <- c(colSums(deriv_H[1:para$X_padding_num, 1:(1+para$XZ_factor), drop = F]),
                  sum(deriv_H[(1+para$X_padding_num):(para$X_padding_num + para$Y_padding_num), 1]),
                  deriv_H[(1+para$X_padding_num + para$Y_padding_num):(dim(deriv_H)[1]), , drop = F])
  return( c(gradient_Phi, gradient_H))
}

multi_factor_padding_COV_directML <- function(para){

  lower_bound_LBFGS <- c(rep(10^-6, 2+length(para$Phi_hat_diag)),
                         rep(-Inf, length(para$alpha)+1+length(c(para$eta_j))))
  upper_bound_LBFGS <- c(rep(1-10^-6, 2+length(para$Phi_hat_diag)),
                         rep(Inf, length(para$alpha)+length(c(para$eta_j))))

  # convergence is determined by b2 = 0 (keep updating it until b2 < 0.0001)
  para_b1b2 <- list()
  para_b1b2$b1_b2 <- rep(1,1+para$XZ_factor)
  # max number of iteration
  max_itr <- 10
  idx <- 1
  while(sum(abs(para_b1b2$b1_b2[2:(para$XZ_factor+1)]) > 10^(-4)) & idx <= max_itr){
    idx <- idx + 1
    estimate_allpar <- optim(par = c(para$sigma_X, para$sigma_Y, para$Phi_hat_diag, para$alpha, para$beta, para$eta_j),
                             fn = function(x) -Marginal_likelihood_cov_multifactor(x, para),
                             gr = function(x) -gradient_likelihood_cov_multifactor(x, para),
                             method = "L-BFGS-B",
                             lower = lower_bound_LBFGS, upper = upper_bound_LBFGS)$par # optim is bounded for the uniqueness (variance term)
    # assign parameters
    para$sigma_X <- estimate_allpar[1]
    para$sigma_Y <- estimate_allpar[2]
    para$Phi_hat_diag <- estimate_allpar[3:(length(para$Phi_hat_diag) + 2)]
    para$Phi_hat_diag_padded <- c(rep(para$sigma_X, para$X_padding_num),
                                  rep(para$sigma_Y, para$Y_padding_num),
                                  para$Phi_hat_diag)
    para$alpha <- estimate_allpar[(length(para$Phi_hat_diag) + 3):(length(para$Phi_hat_diag) + para$XZ_factor + 3)]
    para$beta <- estimate_allpar[(length(para$Phi_hat_diag) + length(para$alpha) + 3)]
    para$eta_j <- matrix(estimate_allpar[(length(para$Phi_hat_diag) + length(para$alpha) + 4): length(estimate_allpar)], ncol = 1+para$XZ_factor + para$Z_factor)

    # assigning the number to the original matrix
    para$H <- rbind(c(para$alpha, rep(0, para$Z_factor)), c(para$beta, rep(0, times = para$XZ_factor + para$Z_factor)), para$eta_j)
    para$H_padded <- rbind(matrix(rep(c(para$alpha, rep(0, para$Z_factor)), each = para$X_padding_num), nrow = para$X_padding_num),
                           matrix(rep(c(para$beta, rep(0, times = para$XZ_factor + para$Z_factor)), each = para$Y_padding_num), nrow = para$Y_padding_num),
                           para$eta_j)
    # check if b1b2 is to the direction (1,0,0...,0,0)
    para_b1b2 <- two_factor_padding_COV_ECME_Estep(para)
    para_b1b2$b1_b2 <- para_b1b2$CX_xi[para$X_padding_num + 1,,drop=F] %*% solve(para_b1b2$Cxi_xi)
    ### line below is helpful for checking the convergence
    # print(para_b1b2$b1_b2)
    # perform the rotation along each hyperplane in turn
    for(f in 1:para$XZ_factor){
      para_rotation <- two_factor_padding_COV_ECME_Estep(para)
      # para_rotation$a1_a2 <- para_rotation$CX_xi[1,,drop=F] %*% solve(para_rotation$Cxi_xi)
      para_rotation$b1_b2 <- para_rotation$CX_xi[para$X_padding_num + 1,,drop=F] %*% solve(para_rotation$Cxi_xi)
      rotation_angle <- atan(-para_rotation$b1_b2[f+1]/para_rotation$b1_b2[1])
      rotation_matrix <- diag(x = 1, nrow = para$XZ_factor + 1)
      rotation_matrix[1,1] <- cos(rotation_angle)
      rotation_matrix[f+1,1] <- sin(rotation_angle)
      rotation_matrix[1,f+1] <- -sin(rotation_angle)
      rotation_matrix[f+1,f+1] <- cos(rotation_angle)
      para$beta <- (rotation_matrix %*% para_rotation$b1_b2[1:(1+para$XZ_factor)])[1]
      para$alpha <- rotation_matrix %*% para$alpha
      para$eta_j[,1:(1+para$XZ_factor)] <- t(rotation_matrix %*% t(para$eta_j[,1:(1+para$XZ_factor)]))

      # assigning the number to the original matrix
      para$H <- rbind(c(para$alpha, rep(0, para$Z_factor)), c(para$beta, rep(0, times = para$XZ_factor + para$Z_factor)), para$eta_j)
      para$H_padded <- rbind(matrix(rep(c(para$alpha, rep(0, para$Z_factor)), each = para$X_padding_num), nrow = para$X_padding_num),
                             matrix(rep(c(para$beta, rep(0, times = para$XZ_factor + para$Z_factor)), each = para$Y_padding_num), nrow = para$Y_padding_num),
                             para$eta_j)
    }
  }
  para$itr_quasi_newton <- idx - 1 # number of rotation and quasi-newton updates
  return(para)
}

# XZ_factor is the number of factors shared by X and Z (but not Y); Z_factor are the factors that are not shared by X, Y
wrapper_multi_factor_padding_COV_directML <- function(X, Y, Z, X_padding_num = 1, Y_padding_num = 1, XZ_factor = 1, Z_factor = 0){
  X <- scale(X)
  Y <- scale(Y)
  Z <- apply(Z, 2, scale)

  # padding X,Y with different padding
  XY_Z_padded <- padding_X_Y_cov(X, Y, Z, X_padding_num = X_padding_num, Y_padding_num = Y_padding_num)

  para <- init_para_COV_factor(XY_Z_padded, X_padding_num, Y_padding_num, XZ_factor, Z_factor)
  para <- multi_factor_padding_COV_directML(para)
  return(para)
}

# bootstrapping samples from the covariance matrix
##########
####### use jackknife instead of bootstrapping (bootstrapping is likely to cause numeric issue as it would have repeated rows)
jk_cov <- function(XY_Z_padded, X_padding_num, Y_padding_num, num_jk_block = 10){
  # divide the ids into num_jk_block numbers: num_jk_block > 1
  size_block <- floor((dim(XY_Z_padded)[1] - X_padding_num - Y_padding_num)/num_jk_block)
  jk_XY_Z_COV <- list()
  for(jk_idx in 1:num_jk_block){
    if(jk_idx == num_jk_block){
      ids_exclude <-  (1+X_padding_num + Y_padding_num + (jk_idx - 1)*size_block ):(dim(XY_Z_padded)[1])
    }else{
      ids_exclude <-  (1+X_padding_num + Y_padding_num + (jk_idx - 1)*size_block ):(X_padding_num + Y_padding_num + jk_idx *size_block )
    }
    ids_include <- setdiff(1:dim(XY_Z_padded)[1], ids_exclude)
    jk_XY_Z_COV[[jk_idx]] <- XY_Z_padded[ids_include,ids_include]
  }
  return(jk_XY_Z_COV)
}
# create jackknife samples for OLS analysis
jk_ols_x2z <- function(OLS_x2z, X_padding_num, Y_padding_num, num_jk_block = 10){
  # divide the ids into num_jk_block numbers
  size_block <- floor((length(OLS_x2z) - X_padding_num - Y_padding_num)/num_jk_block)
  jk_ols_effects <- list()
  for(jk_idx in 1:num_jk_block){
    if(jk_idx == num_jk_block){
      ids_exclude <-  (1+X_padding_num + Y_padding_num + (jk_idx - 1)*size_block ):(length(OLS_x2z))
    }else{
      ids_exclude <-  (1+X_padding_num + Y_padding_num + (jk_idx - 1)*size_block ):(X_padding_num + Y_padding_num + jk_idx *size_block )
    }
    ids_include <- setdiff(1:length(OLS_x2z), ids_exclude)
    jk_ols_effects[[jk_idx]] <- OLS_x2z[ids_include]
  }
  return(jk_ols_effects)
}

bootstrapping_cov <- function(XY_Z_padded,X_padding_num, Y_padding_num){
  bs_idx <- sample((1+X_padding_num + Y_padding_num):dim(XY_Z_padded)[1], replace = T)
  return(XY_Z_padded[c(1:(X_padding_num + Y_padding_num), bs_idx ),
                     c(1:(X_padding_num + Y_padding_num), bs_idx )])
}
# initialise the parameters for the factor analysis
init_para_COV_factor <- function(XY_Z_padded, X_padding_num, Y_padding_num, XZ_factor, Z_factor){
  para <- list()

  para$XY_Z_padded <- XY_Z_padded
  para$X_padding_num <- X_padding_num
  para$Y_padding_num <- Y_padding_num
  para$XZ_factor <- XZ_factor
  para$Z_factor <- Z_factor

  # initialize the parameters: need to initialise all the factors (1+XZ_factor + Z_factor)
  para$alpha <- rep(1,1+para$XZ_factor) # initial alpha as equal to xi1 and xi2
  para$beta <- para$XY_Z_padded[para$X_padding_num + 1,1]
  para$eta_j <- para$XY_Z_padded[(X_padding_num + Y_padding_num + 1):dim(XY_Z_padded)[1],c(1, (X_padding_num + Y_padding_num + 1):(X_padding_num + Y_padding_num + para$XZ_factor + para$Z_factor))]
  para$H <- rbind(c(para$alpha, rep(0, para$Z_factor)), c(para$beta, rep(0, times = para$XZ_factor + para$Z_factor)), para$eta_j)
  para$H_padded <- rbind(matrix(rep(c(para$alpha, rep(0, para$Z_factor)), each = para$X_padding_num), nrow = para$X_padding_num),
                         matrix(rep(c(para$beta, rep(0, times = para$XZ_factor + para$Z_factor)), each = para$Y_padding_num), nrow = para$Y_padding_num),
                         para$eta_j)
  para$sigma_X <- 1- para$beta^2
  para$sigma_Y <- 1- para$beta^2
  para$Phi_hat_diag <- 1 - para$eta_j[,1]^2
  para$Phi_hat_diag_padded <- c(rep(para$sigma_X, para$X_padding_num),
                                rep(para$sigma_Y, para$Y_padding_num),
                                para$Phi_hat_diag)
  return(para)
}

# model selection

#' Factor model selection.
#'
#' `Factor_modelselection_COV_directML()` run different factor model with different number of latent factors.
#' The model select the model structure that maximise the variance explained in target protein. We use a
#' jackknifing procedure to test if adding a new latent factor explains more variance in the target protein.
#'
#' @param X a vector which denote the levels of the target protein.
#' @param Y a vector (same length of X) which denote the levels of the helper protein, should be significantly correlated with X,
#' @param Z a matrix (each row correspond to one element in X and Y) of auxiliary proteins that tags biological variance in X. We recommend at least 100 columns.
#' @param X_padding_num default 1: increasing this number will putting more weights on target protein using a padding procedure. Usually not required as factor model is already controlling latent factors in Z.
#' @param Y_padding_num default 1: increasing this number will putting more weights on helper protein using a padding procedure. Usually not required as factor model is already controlling latent factors in Z.
#' @param max_XZ_factor default 5: grid search up to the maximum number of latent factors of X that are independent of Y.
#' @param max_Z_factor default 5: grid search up to the maximum number of latent factors of Z that are independent of X and Y.
#' @param num_jk_block default 10: jackknife blocks used to test which model explain more variance in X.
#'
#' @return a list containing objects:
#' \describe{
#'   \item best_XZ: The choice of optimal latent factors of X that are independent of Y.
#'   \item best_Z: The choice of optimal latent factors of Z that are independent of X and Y.
#'   \item para_z_factors: comparison of metrics that are used to compare different models, we use sigma_X which denote the residual variance of X.
#'   \item para_z_factors_jk: jackknife results for choosing best_Z.
#'    \item  para_XZ_factors_jk: jackknife results for choosing best_XZ.
#'    \item  summary_jackknife_testing_Z_factor: summary across jackknife blocks for choosing best_Z.
#'    \item  summary_jackknife_testing_XZ: summary across jackknife blocks for choosing best_XZ.
#' }
#'
#' @export
#'
#' @examples
#' library(ExpressionWideRegression)
#' factor_model_selection <- Factor_modelselection_COV_directML(X_EXAMPLE, Y_EXAMPLE, Z_EXAMPLE,
#'  X_padding_num = 1,
#'  Y_padding_num = 1,
#'  max_XZ_factor = 3,
#'  max_Z_factor = 3,
#'  num_jk_block = 5)
#'
Factor_modelselection_COV_directML <- function(X, Y, Z, X_padding_num = 1, Y_padding_num = 1, max_XZ_factor = 5, max_Z_factor = 5, num_jk_block = 10){
  X <- scale(X)
  Y <- scale(Y)
  Z <- apply(Z, 2, scale)

  # padding X,Y with different padding
  XY_Z_padded <- padding_X_Y_cov(X, Y, Z, X_padding_num = X_padding_num, Y_padding_num = Y_padding_num)

  para_z_factors <- data.frame(Z_factor = as.numeric(),
                               XZ_factor = as.numeric(),
                               sigma_X = as.numeric(),
                               sigma_Y = as.numeric(),
                               AIC = as.numeric(),
                               BIC = as.numeric(),
                               Log_likelihood = as.numeric(),
                               ratio_1 = as.numeric(),
                               ratio_2 = as.numeric(),
                               adjusting_ratio = as.numeric())
  XZ_factor <- 1
  # first select best Z_factors using bisection method
  for(Z_factor in 0:max_Z_factor){
    print(paste0("Testing covariate factors = ", Z_factor))
    para <- init_para_COV_factor(XY_Z_padded, X_padding_num, Y_padding_num, XZ_factor, Z_factor)
    para <- multi_factor_padding_COV_directML(para)

    # compuate all factors
    par_H0 <- c(para$sigma_X,
                para$sigma_Y,
                para$Phi_hat_diag,
                para$alpha,
                para$beta,
                para$eta_j)
    LL_H0 <- dim(Z)[1]/2 * Marginal_likelihood_cov_multifactor(par_H0, para = para)

    ratio_1 <- apply(para$eta_j[,2:(1+para$XZ_factor), drop=F], 2, var)/var(para$eta_j[,1])
    ratio_2 <- para$alpha[2:(1+para$XZ_factor)]^2/para$alpha[1]^2
    adjusting_ratio <- (1+sum(ratio_2))/(1+sum(ratio_1 * ratio_2))

    para_z_factors <- para_z_factors %>%
      add_row(Z_factor = Z_factor,
              XZ_factor = XZ_factor,
              sigma_X = para$sigma_X,
              sigma_Y = para$sigma_Y,
              AIC = 2 *  (length(par_H0) - LL_H0),
              BIC = 2 *  (length(par_H0) * log(dim(Z)[1]) - LL_H0),
              Log_likelihood = LL_H0,
              ratio_1 = ratio_1,
              ratio_2 = ratio_2,
              adjusting_ratio = adjusting_ratio)
  }
  # using sigma_X to detect if the model is a better fit
  # perform bootstrapping to test if there are significant difference in sigma_X
  XZ_factor <- 1
  para_z_factors_jk <- data.frame(jk_idx = as.numeric(),
            Z_factor = as.numeric(),
            XZ_factor = as.numeric(),
            sigma_X = as.numeric())

  jk_cov_list <- jk_cov(XY_Z_padded,X_padding_num, Y_padding_num, num_jk_block = num_jk_block)
  for(jk_idx in 1:num_jk_block){
    print(paste0("Test if the factor model explains more variance in X, Jackknife index ", jk_idx))
    # XY_Z_padded_bt_idx <- bootstrapping_cov(XY_Z_padded,X_padding_num, Y_padding_num) # use the same boostrapping covariance matrix
    XY_Z_padded_jk_idx <- jk_cov_list[[jk_idx]]
    for(Z_factor in 0:max_Z_factor){
      para <- init_para_COV_factor(XY_Z_padded_jk_idx, X_padding_num, Y_padding_num, XZ_factor, Z_factor)
      para <- multi_factor_padding_COV_directML(para)
      para_z_factors_jk <- para_z_factors_jk %>%
        add_row(jk_idx = jk_idx,
                Z_factor = Z_factor,
                XZ_factor = XZ_factor,
                sigma_X = para$sigma_X)
    }
  }
  # compute the change in sigma_X
  subtraction_term <- para_z_factors_jk %>%
    mutate(matched_Z_factor_number = Z_factor + 1, sigma_X_baseline = sigma_X)
  summary_jackknife_testing_Z_factor <- para_z_factors_jk %>%
    left_join(select(subtraction_term, jk_idx, matched_Z_factor_number, sigma_X_baseline),
              by = c("Z_factor" = "matched_Z_factor_number", "jk_idx" = "jk_idx")) %>%
    mutate(diff_sigma_X = sigma_X - sigma_X_baseline) %>% # if diff_sigma_X < 0, means the model is capture more variance in X
    group_by(Z_factor, XZ_factor) %>%
    summarise(mean_diff_sigma_X = mean(diff_sigma_X), se_mean_diff_sigma_X =  sqrt(num_jk_block) * sd(diff_sigma_X)) %>%
    mutate(zscore_diff_sigma_X = mean_diff_sigma_X/se_mean_diff_sigma_X) %>%
    mutate(improved_fitting = zscore_diff_sigma_X < -1.645) %>% # one side test 95%
    mutate(improved_fitting = if_else(is.na(improved_fitting), 1, improved_fitting))

  # select the best Z_factors based on jackknife
  summary_jackknife_testing_Z_factor$improved_fitting <- cumprod(summary_jackknife_testing_Z_factor$improved_fitting)
  Z_factor <- summary_jackknife_testing_Z_factor$Z_factor[which.max(pull(filter(summary_jackknife_testing_Z_factor,improved_fitting == 1), Z_factor)) ]

  #select best XZ factor
  for(XZ_factor in 1:max_XZ_factor){
    print(paste0("Testing total X latent factors = ", XZ_factor + 1)) # note minimal = 2 as there is one factor that is shared with Y
    # initialize the parameters: need to initialise all the factors (1+XZ_factor + Z_factor)
    para <- init_para_COV_factor(XY_Z_padded, X_padding_num, Y_padding_num, XZ_factor, Z_factor)
    para <- multi_factor_padding_COV_directML(para)

    # compuate all factors
    par_H0 <- c(para$sigma_X,
                para$sigma_Y,
                para$Phi_hat_diag,
                para$alpha,
                para$beta,
                para$eta_j)
    LL_H0 <- dim(Z)[1]/2 * Marginal_likelihood_cov_multifactor(par_H0, para = para)

    ratio_1 <- apply(para$eta_j[,2:(1+para$XZ_factor), drop=F], 2, var)/var(para$eta_j[,1])
    ratio_2 <- para$alpha[2:(1+para$XZ_factor)]^2/para$alpha[1]^2
    adjusting_ratio <- (1+sum(ratio_2))/(1+sum(ratio_1 * ratio_2))

    para_z_factors <- para_z_factors %>%
      add_row(Z_factor = Z_factor,
              XZ_factor = XZ_factor,
              sigma_X = para$sigma_X,
              sigma_Y = para$sigma_Y,
              AIC = 2 *  (length(par_H0) - LL_H0),
              BIC = 2 *  (length(par_H0) * log(dim(Z)[1]) - LL_H0),
              Log_likelihood = LL_H0,
              ratio_1 = sum(ratio_1),
              ratio_2 = sum(ratio_2),
              adjusting_ratio = adjusting_ratio)

  }

  # using jacknife for selecting best z-number
  para_XZ_factors_jk <- data.frame(jk_idx = as.numeric(),
                                  Z_factor = as.numeric(),
                                  XZ_factor = as.numeric(),
                                  sigma_X = as.numeric())
  for(jk_idx in 1:num_jk_block){
    print(paste0("Test if there are multiple latent factors in X, Jackknife index ", jk_idx))
    # XY_Z_padded_bt_idx <- bootstrapping_cov(XY_Z_padded,X_padding_num, Y_padding_num) # use the same boostrapping covariance matrix
    XY_Z_padded_jk_idx <- jk_cov_list[[jk_idx]]
    for(XZ_factor in 1:max_XZ_factor){
      para <- init_para_COV_factor(XY_Z_padded_jk_idx, X_padding_num, Y_padding_num, XZ_factor, Z_factor)
      para <- multi_factor_padding_COV_directML(para)
      para_XZ_factors_jk <- para_XZ_factors_jk %>%
        add_row(jk_idx = jk_idx,
                Z_factor = Z_factor,
                XZ_factor = XZ_factor,
                sigma_X = para$sigma_X)
    }
  }
  subtraction_term <- para_XZ_factors_jk %>%
    mutate(matched_XZ_factor_number = XZ_factor + 1, sigma_X_baseline = sigma_X)
  summary_jackknife_testing_XZ <- para_XZ_factors_jk %>%
    left_join(select(subtraction_term, jk_idx, matched_XZ_factor_number, sigma_X_baseline),
              by = c("XZ_factor" = "matched_XZ_factor_number", "jk_idx" = "jk_idx")) %>%
    mutate(diff_sigma_X = sigma_X - sigma_X_baseline) %>% # if diff_sigma_X < 0, means the model is capture more variance in X
    group_by(Z_factor,XZ_factor) %>%
    summarise(mean_diff_sigma_X = mean(diff_sigma_X), se_mean_diff_sigma_X =  sqrt(num_jk_block) * sd(diff_sigma_X)) %>%
    mutate(zscore_diff_sigma_X = mean_diff_sigma_X/se_mean_diff_sigma_X) %>%
    mutate(improved_fitting = zscore_diff_sigma_X < -1.645) %>%
    mutate(improved_fitting = if_else(is.na(improved_fitting), 1, improved_fitting))
  # select the best Z_factors based on jackknife
  summary_jackknife_testing_XZ$improved_fitting <- cumprod(summary_jackknife_testing_XZ$improved_fitting)
  best_XZ <- summary_jackknife_testing_XZ$XZ_factor[which.max(pull(filter(summary_jackknife_testing_XZ,improved_fitting == 1), XZ_factor)) ]

  factor_model_selection <- list()
  factor_model_selection$best_XZ <- best_XZ
  factor_model_selection$best_Z <- Z_factor
  factor_model_selection$para_z_factors <- para_z_factors
  factor_model_selection$para_z_factors_jk <- para_z_factors_jk
  factor_model_selection$para_XZ_factors_jk <- para_XZ_factors_jk
  factor_model_selection$summary_jackknife_testing_Z_factor <- summary_jackknife_testing_Z_factor
  factor_model_selection$summary_jackknife_testing_XZ <- summary_jackknife_testing_XZ
  return(factor_model_selection)
}

# this function only uses the summary level data to perform fast jackknifing
jackknife_acrossZ_inputCOV <- function(XY_Z_padded, OLS_y2z, OLS_x2z, PxPy,
                                       XZ_factor,
                                       Z_factor,
                                       num_jk_block = 20,
                                       X_padding_num = 1,
                                       Y_padding_num = 1
                                       ){

  jk_cov_list <- jk_cov(XY_Z_padded,X_padding_num, Y_padding_num, num_jk_block = num_jk_block)
  jk_OLS_y2z <- jk_ols_x2z(OLS_y2z,X_padding_num, Y_padding_num, num_jk_block = num_jk_block)
  jk_OLS_x2z <- jk_ols_x2z(OLS_x2z,X_padding_num, Y_padding_num, num_jk_block = num_jk_block)

  jk_results <- data.frame(jk_idx = as.numeric(),
                           beta_EWR = as.numeric(),
             approx_Bioprop_xi1_xi2 = as.numeric(),
             ratio_1 = as.numeric(),
             ratio_2 = as.numeric(),
             adjusting_ratio = as.numeric(),
             Bioprop_xi1_xi2= as.numeric())
  for(jk_idx in 1:num_jk_block){
    print(paste("Jackknife index: ", jk_idx))
    # beta_EWR
    OLS_x2z_jk <- jk_OLS_x2z[[jk_idx]]
    OLS_y2z_jk <- jk_OLS_y2z[[jk_idx]]
    XY_Z_padded_jk <- jk_cov_list[[jk_idx]]
    beta_EWR <- summary(lm(OLS_y2z_jk ~ OLS_x2z_jk))$coefficients[2,1]
    approx_Bioprop_xi1_xi2 <- sqrt(PxPy/beta_EWR^2)

    para_jk <- init_para_COV_factor(XY_Z_padded_jk, X_padding_num, Y_padding_num, XZ_factor, Z_factor)
    para_jk <- multi_factor_padding_COV_directML(para_jk)

    ratio_1 <- apply(para_jk$eta_j[,2:(1+para_jk$XZ_factor), drop=F], 2, var)/var(para_jk$eta_j[,1])
    ratio_2 <- para_jk$alpha[2:(1+para_jk$XZ_factor)]^2/para_jk$alpha[1]^2
    adjusting_ratio <- (1+sum(ratio_2))/(1+sum(ratio_1 * ratio_2))
    Bioprop_xi1_xi2 <- adjusting_ratio * approx_Bioprop_xi1_xi2

    jk_results <- jk_results %>%
      add_row(jk_idx = jk_idx,
              beta_EWR = beta_EWR,
              approx_Bioprop_xi1_xi2 = approx_Bioprop_xi1_xi2,
              ratio_1 = sum(ratio_1),
              ratio_2 = sum(ratio_2),
              adjusting_ratio = adjusting_ratio,
              Bioprop_xi1_xi2= Bioprop_xi1_xi2)
  }
  return(jk_results)
}

# boot strapping step uses individual level data: repeat teach step with new X,Y,Z
Bootstrap_step <- function(X, Y, Z,  X_padding, number_XZ_factor_model, covariate_Z_factor){
  J <- dim(Z)[2]
  # two step procedure
  PxPy <- summary(lm(scale(X) ~ scale(Y)))$coefficients[2,1]^2
  # beta_EWR
  OLS_x2z <- c()
  OLS_y2z <- c()
  for(j in 1:J){
    x2z_coef <- summary(lm(scale(Z[,j]) ~ scale(X)))$coefficients[2,1]
    y2z_coef <- summary(lm(scale(Z[,j]) ~ scale(Y)))$coefficients[2,1]

    OLS_x2z <- c(OLS_x2z, x2z_coef)
    OLS_y2z <- c(OLS_y2z, y2z_coef)
  }
  beta_EWR <- summary(lm(OLS_y2z ~ OLS_x2z))$coefficients[2,1]
  approx_Bioprop_xi1_xi2 <- sqrt(PxPy/beta_EWR^2)

  para_directML <- wrapper_multi_factor_padding_COV_directML(X, Y, Z,
                                                             X_padding_num = X_padding,
                                                             Y_padding_num = 1,
                                                             XZ_factor = number_XZ_factor_model,
                                                             Z_factor = covariate_Z_factor)
  ratio_1 <- var(para_directML$eta_j[,2])/var(para_directML$eta_j[,1])
  ratio_2 <- para_directML$alpha[2]^2/para_directML$alpha[1]^2
  Bioprop_xi1_xi2 <- (1+ratio_2)/(1+ratio_1 * ratio_2) * approx_Bioprop_xi1_xi2
  return(Bioprop_xi1_xi2)
}


# functions to compute the bio-proportion

#' Biological proportion varaince computation
#'
#' @param X a vector which denote the levels of the target protein.
#' @param Y a vector (same length of X) which denote the levels of the helper protein, should be significantly correlated with X.
#' @param Z a matrix (each row correspond to one element in X and Y) of auxiliary proteins that tags biological variance in X. We recommend at least 100 columns.
#' @param X_padding default 1: increasing this number will putting more weights on target protein using a padding procedure. Usually not required as factor model is already controlling latent factors in Z.
#' @param PC_control_num default 10: the number of PCs in Z that would be regressed out from Z.
#' @param bootstrap_num Default 20. Number of Bootstrapping or Jackknifing blocks used to compute the confidence interval of the estimate.
#' @param bootstrap_samples_flg Default F: when TRUE, bootstrapping across individuals; recommend FALSE where sampling is done across the Z columns, this is more calibrated as EWR is estimated across columns of Zs.
#' @param number_XZ_factor_model Default 1: the choice of latent factors of X that are independent of Y.
#' @param covariate_Z_factor Default 1: The choice of optimal latent factors of Z that are independent of X and Y; use one factor to control the latent structure in Z.
#'
#' @return a list:
#' Bioprop_xi1_xi2: estimate of biological proportion (i.e. 1-noise) variance in protein level X.
#' sampling_se: standard error of the estimate Bioprop_xi1_xi2; estimated through sampling procedure.
#' sampling_mean: mean across the sample estimate (could be used as a bias detection tool from Jackknifing).
#' PC_variance: the proportion of X variance tagged by PCs of Z, which are regressed out from Z and then added back to the final estimate.
#' PxPy
#' PxPy_se: standard error of PxPy estimate.
#' beta_EWR: Slope of covariance y2z on x2z.
#' beta_EWR_se: standard error of beta_EWR estimate.
#' approx_Bioprop_xi1_xi2: term A in the paper.
#' ratio_2: term R2 in the paper.
#' ratio_1: term R2=1 in the paper.
#' sigma_X: residual variance in X that is not captured by the latent factors.
#' sigma_Y: residual variance in Y that is not captured by the latent factors.
#' itr_quasi_newton: the number of iteration between the (i) quasi-newton MLE and (ii) rotation of factor loadings. See the paper supplementary for details.
#'
#'
#' @export
#'
#' @examples
#' library(ExpressionWideRegression)
#' full_run_results_big_wrapper <- EWR(X_EXAMPLE,Z_EXAMPLE)
#'
bioProp_EWR_COVfactor_adjustment <- function(X, Y, Z,  X_padding = 1, PC_control_num = 10, bootstrap_num = 20, bootstrap_samples_flg = F,
                                             number_XZ_factor_model = 1, covariate_Z_factor = 1){ # X is Nx1, Y is Nx1, Z is NxJ
  J <- dim(Z)[2]
  print(paste0("Number of auxiliary proteins: ", J, "; Recommended at auxiliary protein > 100"))

  if(PC_control_num > 0){
    print(paste0("Compute PCA of auxiliary proteins to detect latent structure."))
    X <- scale(X)
    Y <- scale(Y)
    Z <- apply(Z, 2, scale)
    Z_PCA <- Z
    Z_PCA[is.na(Z)] <- 0
    pca_Z <- prcomp(Z_PCA, scale = F, rank = PC_control_num)
    top_PCs <- pca_Z$x[,1:PC_control_num]
    # only regress out PCs with low correlation with X
    # we regress out only the variable that explain relative small amount of variance in X (i.e. R2 of PC_X < 0.01)
    # ratio_Variance_PC_vs_X <- 0.1 # we use a ratio min(0.1/#PCs, 0.01) to denote relatively small proportion
    # R2_X_topPCs <- which(cor(X, top_PCs, use = "pairwise.complete.obs")^2 < min(ratio_Variance_PC_vs_X/PC_control_num, 0.01) )
    R2_X_topPCs <- which(cor(X, top_PCs, use = "pairwise.complete.obs")^2 < 0.01 )
    # R2_X_topPCs <- which(cor(X, top_PCs, use = "pairwise.complete.obs")^2 < (ratio_Variance_PC_vs_X* (pca_Z$sdev^2)[1:PC_control_num]/dim(Z)[2]) )
    print(paste("Regressing out PC", paste(R2_X_topPCs, collapse = ", "), " from the measurements"))
    lm_md_PC <- lm(X~top_PCs[,R2_X_topPCs],na.action="na.exclude")
    PC_variance <- summary(lm_md_PC)$r.square
    ####################################################################################################
    # we shouldn't regress out PCs from X,Y, only regress out from Z
    ####################################################################################################
    for(j in 1:J){
      Z[,j] <- resid(lm(Z[,j]~top_PCs[,R2_X_topPCs],na.action="na.exclude"))
    }
    print(paste("Variance in PC:", PC_variance))
  }

  # two step procedure
  print(paste("Perform covariance regression."))
  PxPy_model <- summary(lm(scale(X) ~ scale(Y)))
  PxPy <- PxPy_model$coefficients[2,1]^2
  PxPy_se <- PxPy_model$coefficients[2,2]
  # beta_EWR
  OLS_x2z <- c()
  OLS_y2z <- c()
  for(j in 1:J){
    x2z_coef <- summary(lm(scale(Z[,j]) ~ scale(X)))$coefficients[2,1]
    y2z_coef <- summary(lm(scale(Z[,j]) ~ scale(Y)))$coefficients[2,1]

    OLS_x2z <- c(OLS_x2z, x2z_coef)
    OLS_y2z <- c(OLS_y2z, y2z_coef)
  }
  beta_EWR <- summary(lm(OLS_y2z ~ OLS_x2z))$coefficients[2,1]
  beta_EWR_se <- summary(lm(OLS_y2z ~ OLS_x2z))$coefficients[2,2]
  print(paste("Slope of covariance y2z on x2z is: ", beta_EWR, "(", beta_EWR_se, ")", " with z-score ", beta_EWR/beta_EWR_se, ". |Z-score| > 5 required for the estimate to be approximately unbiased."))

  approx_Bioprop_xi1_xi2 <- sqrt(PxPy/beta_EWR^2)

  #####################
  # provide an option to perform factor analysis model selection.
  #####################
  # factor analysis
  print(paste("Use factor anlaysis to estimate the correction term."))
  # number_XZ_factor_model <- 1 # modeling a single latent factor that are not shared with Y
  # covariate_Z_factor <- 1 # number of factors used to control the structure in Z that are not related to X (similar to PCA)
  para_directML <- wrapper_multi_factor_padding_COV_directML(X, Y, Z,
                                                             X_padding_num = X_padding,
                                                             Y_padding_num = 1,
                                                             XZ_factor = number_XZ_factor_model,
                                                             Z_factor = covariate_Z_factor)
  # checking if the latent factors of Z also explain significant variance in X
  para_rotation <- two_factor_padding_COV_ECME_Estep(para_directML)
  para_rotation$a1_a2 <- para_rotation$CX_xi[1,,drop=F] %*% solve(para_rotation$Cxi_xi)
  print(paste("Covariate latent factors explain ", sum(para_rotation$a1_a2[(2+para_rotation$XZ_factor):(1+para_rotation$XZ_factor + para_rotation$Z_factor)]^2), " variance of X"))

  ratio_1 <- apply(para_directML$eta_j[,2:(1+para_directML$XZ_factor), drop=F], 2, var)/var(para_directML$eta_j[,1])
  ratio_2 <- para_directML$alpha[2:(1+para_directML$XZ_factor)]^2/para_directML$alpha[1]^2
  Bioprop_xi1_xi2 <- (1+sum(ratio_2))/(1+sum(ratio_1 * ratio_2)) * approx_Bioprop_xi1_xi2

  # adjusting the ratio for the PCs
  Bioprop_xi1_xi2 <- Bioprop_xi1_xi2 + PC_variance

  # bootstrap across individuals
  if(bootstrap_samples_flg){
    print(paste("Total bootstrapping samples: ", bootstrap_num))
    print(paste("Bootstrapping across samples, recomputing the covariance structure (slower)."))
    bt_estimat_list <- c()
    if(bootstrap_num > 0){
      for(bt in 1:bootstrap_num){
        print(paste("Bootsample: ", bt))
        boot_sample_idx <- sample(1:length(X), size = length(X), replace = T)
        bt_estimate <- Bootstrap_step(X = X[boot_sample_idx], Y = Y[boot_sample_idx], Z=Z[boot_sample_idx, ],  X_padding = X_padding,
                                      number_XZ_factor_model = number_XZ_factor_model, covariate_Z_factor = covariate_Z_factor)
        bt_estimat_list <- c(bt_estimat_list, bt_estimate)
      }
    }
    sampling_se <- sd(bt_estimat_list)
    sampling_mean <- mean(bt_estimat_list) + PC_variance
  }else{
    if(bootstrap_num > 1){
      print(paste("Jackknifing auxiliary proteins Z; directly reuse sufficient statistics (covariance). "))
      jk_results <- jackknife_acrossZ_inputCOV(para_directML$XY_Z_padded, OLS_y2z, OLS_x2z, PxPy,
                                               X_padding_num = 1,
                                               Y_padding_num = 1,
                                               num_jk_block = bootstrap_num,
                                               XZ_factor = number_XZ_factor_model,
                                               Z_factor = covariate_Z_factor)

      sampling_se <- sqrt(bootstrap_num) * sd(jk_results$Bioprop_xi1_xi2)
      sampling_mean <- mean(jk_results$Bioprop_xi1_xi2) + PC_variance
    }else{
      print("Jackknife blocks could not be less than 2.")
      sampling_se <- NULL
      sampling_mean <- NULL
    }

  }

  results_biprop <- list()
  results_biprop$Bioprop_xi1_xi2 <- Bioprop_xi1_xi2
  results_biprop$sampling_se <- sampling_se # jackknife (or boostraping) se
  results_biprop$sampling_mean <- sampling_mean # jackknife (or boostraping) mean, use to identify bias
  results_biprop$PC_variance <- PC_variance
  results_biprop$PxPy <- PxPy
  results_biprop$PxPy_se <- PxPy_se
  results_biprop$beta_EWR <- beta_EWR
  results_biprop$beta_EWR_se <- beta_EWR_se
  results_biprop$approx_Bioprop_xi1_xi2 <- approx_Bioprop_xi1_xi2
  results_biprop$ratio_2 <- sum(ratio_2)
  results_biprop$ratio_1 <- sum(ratio_1)
  results_biprop$sigma_X <- para_directML$sigma_X
  results_biprop$sigma_Y <- para_directML$sigma_Y
  results_biprop$itr_quasi_newton <- para_directML$itr_quasi_newton

  return(results_biprop)
}

# wrapper function of EWR: input is X and Z. will select Y and iterate through it based on correlation with the target protein.

#' Estimate non-noise proportion of variance in the protein measurement, using a set of auxiliary proteins.
#'
#' @param target_protein: a vector which denote the levels of the target protein.
#' @param auxiliary_proteins: matrix (each row correspond to one element in target proteins) of auxiliary proteins that tags biological variance in target proteins. We recommend at least 100 columns.
#' @param use_as_helper: default 10: we select top #use_as_helper strongly associated columns in auxiliary_proteins as helper proteins for the estimation. Each estimate should be estimate the same value.
#' @param Factor_Model_selection: default TRUE: true means performing a model selection to choose the model structure. If FALSE, default or user specified values will be used.
#' @param XZ_factor: Default 1: the choice of latent factors of X that are independent of Y.
#' @param Z_factor: Default 0: the choice of latent factors of Z that are independent of X and Y.
#'
#' @return list:
#' Bioprop_xi1_xi2: estimate of biological proportion (i.e. 1-noise) variance in protein level X.
#' sampling_se: standard error of the estimate Bioprop_xi1_xi2; estimated through sampling procedure.
#' sampling_mean: mean across the sample estimate (could be used as a bias detection tool from Jackknifing).
#' PC_variance: the proportion of X variance tagged by PCs of Z, which are regressed out from Z and then added back to the final estimate.
#' PxPy
#' PxPy_se: standard error of PxPy estimate.
#' beta_EWR: Slope of covariance y2z on x2z.
#' beta_EWR_se: standard error of beta_EWR estimate.
#' approx_Bioprop_xi1_xi2: term A in the paper.
#' ratio_2: term R2 in the paper.
#' ratio_1: term R2=1 in the paper.
#' sigma_X: residual variance in X that is not captured by the latent factors.
#' sigma_Y: residual variance in Y that is not captured by the latent factors.
#' itr_quasi_newton: the number of iteration between the (i) quasi-newton MLE and (ii) rotation of factor loadings. See the paper supplementary for details.
#' Y_idx: the column index of auxiliary_proteins that are chosen as the helper protein.
#' XZ_factor: the choice of latent factors of X that are independent of Y.
#' Z_factor: the choice of latent factors of Z that are independent of X and Y.
#'
#' @export
#'
#' @examples
#' full_run_results_big_wrapper <- EWR(X_EXAMPLE,Z_EXAMPLE, use_as_helper = 2)
EWR <- function(target_protein, auxiliary_proteins, use_as_helper = 10, Factor_Model_selection = T, XZ_factor = 1, Z_factor = 0){
  target_protein <- c(target_protein)
  if(length(target_protein) != dim(auxiliary_proteins)[1]){
    stop("dimension of target protein and auxiliary proteins is not matching.")
  }

  correlation_target_aux <- cor(target_protein, auxiliary_proteins, use = "pairwise.complete.obs")

  if(max(abs(correlation_target_aux)) > 0.99){
    stop("Really high correlation between target protein and on auxiliary protein. Did you accidentally put the target protein in the helper proteins?")
  }

  if(max(abs(correlation_target_aux)) < 0.05){
    stop("target protein doesn't seem to be correlated with any of auxiliary proteins. Have you checked the ordering is consistent between target and auxiliary proteins?")
  }

  Y_idx_list <- order(abs(correlation_target_aux), decreasing = T)[1:use_as_helper]

  X <- target_protein
  if(Factor_Model_selection){
    Y_idx <- Y_idx_list[1]
    Y <- auxiliary_proteins[, Y_idx]
    Z <- auxiliary_proteins[, -Y_idx]
    print(paste0("performing model selection."))
    factor_model_selection <- Factor_modelselection_COV_directML(X, Y, Z,
                                                                 X_padding_num = 1,
                                                                 Y_padding_num = 1,
                                                                 max_XZ_factor = 3,
                                                                 max_Z_factor = 3,
                                                                 num_jk_block = 5)
    XZ_factor <- factor_model_selection$best_XZ
    Z_factor <- factor_model_selection$best_Z
  }else{ # much faster, no model selection
    XZ_factor <- XZ_factor
    Z_factor <- Z_factor
  }

  bioProp_results <- list()
  for(idx in 1:length(Y_idx_list)){
    print(paste0("Helper protein idx: ", idx))
    Y_idx <- Y_idx_list[idx]
    Y <- auxiliary_proteins[, Y_idx]
    Z <- auxiliary_proteins[, -Y_idx]
    bioProp_results[[idx]] <- dplyr::bind_rows(bioProp_EWR_COVfactor_adjustment(X, Y, Z,  X_padding = 1, bootstrap_num = 20,
                                                                         number_XZ_factor_model = XZ_factor, covariate_Z_factor = Z_factor)) %>%
      mutate(Y_idx = Y_idx, XZ_factor = XZ_factor, Z_factor = Z_factor)
  }
  bioProp_results <- dplyr::bind_rows(bioProp_results)
  return(bioProp_results)
}

