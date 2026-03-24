
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
})


