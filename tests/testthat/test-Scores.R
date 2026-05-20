
set.seed(111)
X <- Iscores:::random_mcar_data(50, 3, 0.2)
methods_list <- list(exp = Iscores:::exp_imputation,
                     norm = Iscores:::norm_imputation)


test_that("We can compare IScores", {
  res <- compare_Iscores(X, methods_list = methods_list)

  expect_identical(round(res[, 1], 5), c(0.45178, 0.56466, 1.29138, 2.14579))

  expect_identical(res[, 2], c("energy_IScore", "energy_IScore",
                               "DR_IScore", "DR_IScore"))

  expect_identical(res[, 3], c("exp", "norm", "exp", "norm"))

})

