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

