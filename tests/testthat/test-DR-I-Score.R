
set.seed(111)
X <- Iscores:::random_mcar_data(100, 3, 0.2)
imputation_func <- Iscores:::exp_imputation

test_that("DR I Score works", {

  set.seed(123)

  res <- DR_IScore(X, imputation_func)
  expect_equal(round(res, 4), 2.5068)

  set.seed(1)

  X_imp <- lapply(1:5, function(i) { Iscores:::norm_imputation(X) })

  set.seed(123)

  res <- DR_IScore(X, X_imp = X_imp)
  expect_equal(round(res, 4), 3.5289)


  set.seed(111)

  X <- Iscores:::random_mcar_data(20, 4, 0.2)

  imputation_func <- Iscores:::exp_imputation

  res <- DR_IScore(X, imputation_func)

  expect_equal(round(res, 4), -1.8155)

  set.seed(111)

  X <- Iscores:::random_mcar_data(20, 6, 0.1)
  X[2, 1] <- 0.1

  imputation_func <- Iscores:::exp_imputation
  res <- DR_IScore(X, imputation_func)

  expect_equal(round(res, 4), 3.1781)

})


test_that("DR I Score throws an error when thee is no imputed data and imputation function", {
  expect_error(DR_IScore(X), "You must provide one of imputation_func or X_imp!")
})


test_that("DR I Score throws an error when imputation function fails", {
  imp_fun <- stop
  expect_error(DR_IScore(X, imp_fun), "Errored imputing X using provided imputation_func!")
})

