test_that("Factor model selection: number of optimal factors for Z and Y", {
  factor_model_selection <- Factor_modelselection_COV_directML(X_EXAMPLE, Y_EXAMPLE, Z_EXAMPLE,
                                                               X_padding_num = 1,
                                                               Y_padding_num = 1,
                                                               max_XZ_factor = 3,
                                                               max_Z_factor = 3,
                                                               num_jk_block = 5)

})

test_that("Run through the inner wrapper", {
  set.seed(19940110)
  test_results <- bioProp_EWR_COVfactor_adjustment(X_EXAMPLE, Y_EXAMPLE, Z_EXAMPLE,  X_padding = 1, PC_control_num = 10, bootstrap_num = 2, bootstrap_samples_flg = F,
                                                   number_XZ_factor_model = 1, covariate_Z_factor = 0)
  expect_true(test_results$Bioprop_xi1_xi2 < 0.31)
  expect_true(test_results$Bioprop_xi1_xi2 > 0.3)
  test_random_X <- bioProp_EWR_COVfactor_adjustment(rnorm(dim(Z_EXAMPLE)[1]), Y_EXAMPLE, Z_EXAMPLE,  X_padding = 1, PC_control_num = 10, bootstrap_num = 2, bootstrap_samples_flg = F,
                                                   number_XZ_factor_model = 1, covariate_Z_factor = 0)
  expect_true(test_random_X$Bioprop_xi1_xi2 < 0.05)
})

test_that("Test the high-level wrapper. Do we got the expected errors?", {
  expect_error(EWR(X_EXAMPLE, cbind(X_EXAMPLE, Z_EXAMPLE)),
               "Really high correlation between target protein and on auxiliary protein. Did you accidentally put the target protein in the helper proteins?")
  expect_error(EWR(X_EXAMPLE[1:(length(X_EXAMPLE)-1)],  Z_EXAMPLE),
               "dimension of target protein and auxiliary proteins is not matching.")
  expect_error(EWR(rnorm(dim(Z_EXAMPLE)[1]),  Z_EXAMPLE[,1:10]),
               "target protein doesn't seem to be correlated with any of auxiliary proteins. Have you checked the ordering is consistent between target and auxiliary proteins?")
})

test_that("Test the high-level wrapper, both model should produce same results", {
  set.seed(19940110)
  full_run_results_big_wrapper <- EWR(X_EXAMPLE,Z_EXAMPLE, use_as_helper = 2)
  full_run_results_big_wrapper_noModel <- EWR(X_EXAMPLE,Z_EXAMPLE, use_as_helper = 2, Factor_Model_selection = F)
  expect_equal(full_run_results_big_wrapper, full_run_results_big_wrapper_noModel)
})
test_that("Test the fast jackknifing procedure", {
  set.seed(19940110)
  bootstrap_num <- 20
  X <- scale(X_EXAMPLE)
  Y <- scale(Y_EXAMPLE)
  Z <- apply(Z_EXAMPLE, 2, scale)
  J <- dim(Z)[2]
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

  print(paste("Use factor anlaysis to estimate the correction term."))
  para_directML <- wrapper_multi_factor_padding_COV_directML(X, Y, Z,
                                                             X_padding_num = 1,
                                                             Y_padding_num = 1,
                                                             XZ_factor = 1,
                                                             Z_factor = 0)
  # checking if the latent factors of Z also explain significant variance in X
  para_rotation <- two_factor_padding_COV_ECME_Estep(para_directML)
  para_rotation$a1_a2 <- para_rotation$CX_xi[1,,drop=F] %*% solve(para_rotation$Cxi_xi)
  print(paste("Covariate latent factors explain ", sum(para_rotation$a1_a2[(2+para_rotation$XZ_factor):(1+para_rotation$XZ_factor + para_rotation$Z_factor)]^2), " variance of X"))

  ratio_1 <- apply(para_directML$eta_j[,2:(1+para_directML$XZ_factor), drop=F], 2, var)/var(para_directML$eta_j[,1])
  ratio_2 <- para_directML$alpha[2:(1+para_directML$XZ_factor)]^2/para_directML$alpha[1]^2
  Bioprop_xi1_xi2 <- (1+sum(ratio_2))/(1+sum(ratio_1 * ratio_2)) * approx_Bioprop_xi1_xi2

  # adjusting the ratio for the PCs
  Bioprop_xi1_xi2 <- Bioprop_xi1_xi2

  print(paste("Jackknifing auxiliary proteins Z; directly reuse sufficient statistics (covariance). "))
  jk_cov_list <- jk_cov(para_directML$XY_Z_padded, 1, 1, num_jk_block = 20)
  jk_OLS_y2z <- jk_ols_x2z(OLS_y2z,num_jk_block = 20)
  jk_OLS_x2z <- jk_ols_x2z(OLS_x2z, num_jk_block = 20)

  expect_equal( jk_cov_list[[1]][3:dim(jk_cov_list[[1]])[1], 1], para_directML$XY_Z_padded[13:dim(para_directML$XY_Z_padded)[1], 1])
  expect_equal( jk_OLS_x2z[[1]], OLS_x2z[11:length(OLS_x2z)])

  # test max diff
  max_diff <- max(abs(para_directML$XY_Z_padded[3:dim(para_directML$XY_Z_padded)[1], 1]-OLS_x2z))
  block_id <- 1
  expect_equal(1,  prod(abs(jk_cov_list[[block_id]][3:dim(jk_cov_list[[block_id]])[1], 1]- jk_OLS_x2z[[block_id]]) <= max_diff) )

  # jk_results <- jackknife_acrossZ_inputCOV(para_directML$XY_Z_padded, OLS_y2z, OLS_x2z, PxPy,
  #                                          X_padding_num = 1,
  #                                          Y_padding_num = 1,
  #                                          num_jk_block = 20,
  #                                          XZ_factor = 1,
  #                                          Z_factor = 0)
  #
  # sampling_se <- sqrt(bootstrap_num) * sd(jk_results$Bioprop_xi1_xi2)
  # sampling_mean <- mean(jk_results$Bioprop_xi1_xi2)
})

