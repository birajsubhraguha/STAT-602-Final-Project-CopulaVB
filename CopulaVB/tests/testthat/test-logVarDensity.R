test_that("derivative wrt rho is zero when means and arguments are zero", {
  out = log_var_density(z = c(0, 0), mu_mat = matrix(0, ncol = 4, nrow = 2),
                        std_mat = matrix(rexp(8), ncol = 4, nrow = 2), rho = 0)
  expect_equal(as.numeric(out$del_rho), 0)
  expect_equal(out$del_mu, matrix(0, ncol = 4, nrow = 2))
  expect_equal(out$del_std, matrix(0, ncol = 4, nrow = 2))
})

test_that("log density matches with that of standard normal of correct dimension with suitable parameters", {
  out = log_var_density(z = rep(0, 3), mu_mat = matrix(0, ncol = 2, nrow = 3),
                        std_mat = matrix(rep(1, 6), ncol = 2, nrow = 3), rho = 0)
  expect_equal(as.numeric(out$log_density), -1.5 * log(2 * pi))
})

