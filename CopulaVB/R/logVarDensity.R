## this defines some useful quantities to be used later
trans_values = function(z, mu_mat, std_mat) {
  p = nrow(mu_mat)
  K = ncol(mu_mat)
  z_mat = matrix(z, nrow = p, ncol = K)
  arg_mat = (z_mat - mu_mat)/std_mat

  f_vec = (rowSums(dnorm(arg_mat)))/K
  F_vec = (rowSums(pnorm(arg_mat)))/K
  Y = qnorm(F_vec)
  return(list(arg_mat = arg_mat, f_vec = f_vec, F_vec = F_vec, Y = Y))
}

## this defines the log variational density
## INPUT
# z = p-dimensional argument, rho = equicorrelation coefficient
# mu_mat, std_mat = matrix of the Gaussian means/standard deviations of size p x K
## OUTPUT
# log density at z, derivatives of log density wrt all the parameters at z
log_var_density = function(z, mu_mat, std_mat, rho) {
  if(sum(dim(mu_mat) == dim(std_mat)) != 2) {
    stop("matrices mu_mat and std_mat should have the same dimensions")
  }
  p = nrow(mu_mat)
  K = ncol(mu_mat)
  val = trans_values(z, mu_mat, std_mat)
  arg_mat = val$arg_mat
  f_vec = val$f_vec
  F_vec = val$F_vec
  Y = val$Y

  f_mat = matrix(f_vec, ncol = K, nrow = p)

  log_det = ((p - 1) * log(1 - rho)) + log(1 + ((p - 1) * rho))

  log_marginals = log(f_vec)
  Y_term = (crossprod(Y) - ((sum(Y)^2)/(1 + ((p - 1) * rho))))
  quad_term = (rho/(1 - rho)) * Y_term

  # returning log of density
  log_density = (-0.5 * (log_det + quad_term)) + sum(log_marginals)

  const_vec = (rho / (rho - 1)) * Y_term * (1/(dnorm(Y))) / K
  const_mat = matrix(const_vec, ncol = K, nrow = p)
  del_1 = const_mat * dnorm(arg_mat) / std_mat
  del_2 = (arg_mat * dnorm(arg_mat) / (std_mat * f_mat))/K

  # derivative wrt mu
  del_mu_mat = del_1 + del_2

  # derivative wrt to std
  del_std_mat = arg_mat * (del_2 - del_1)

  del_rho_quad = crossprod(Y)/((1 - rho)^2) +
    ((sum(Y)^2 * (1 + (p - 1) * (rho^2))) / (((1 - rho)^2) * ((1 + (p - 1) * rho)^2)))
  del_rho_det = - (p * (p - 1) * rho) / ((1 - rho) * (1 + (p - 1) * rho))

  # derivative wrt rho
  del_rho = -.5 * (del_rho_det + del_rho_quad)

  return(list(log_density = log_density, del_mu_mat = del_mu_mat, del_std_mat = del_std_mat, del_rho = del_rho))
}






