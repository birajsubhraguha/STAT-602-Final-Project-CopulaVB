test_that("marginal means are recovered for large samples", {
  out = generate_var_density(n_sample = 10000, mu_mat = matrix(1:8, nrow = 2, byrow = T),
                       std_mat = matrix(1, ncol = 4, nrow = 2), rho = .5, seed = 123)
  means = colMeans(out)
  expect_equal(max(abs(means - c(2.5, 6.5))), 0, tolerance = 1e-2)
})

test_that("sample correaltion tends to zero if copula correlation rho is zero", {
  out = generate_var_density(n_sample = 10000, mu_mat = matrix(rnorm(10, 3, 2), nrow = 2, ncol = 5),
                             std_mat = matrix(rexp(10), ncol = 5, nrow = 2), rho = 0, seed = 321)
  sample_correlation = cor(out[, 1], out[, 2])
  expect_equal(sample_correlation, 0, tolerance = 1e-2)
})
