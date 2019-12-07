## this defines the log variational density as well as the derivatives
## INPUT
# z = p-dimensional argument, rho = equicorrelation coefficient
# mu_mat, std_mat = matrix of the Gaussian means/standard deviations of size p x K
## OUTPUT
# log density at z, derivatives of log density wrt all the parameters at z

#' log density and derivatives
#'
#' \code{log_var_density} Calculates the log value of the Gaussian Copula model density at the provided input, using mixture Gaussian marginals and equicorrelation coefficient rho. Also outputs derivatives of the log density wrt to all the mean and sd parameters, as well as rho
#'
#' @param z : input
#' @param mu_mat : matrix of Gaussian means, components along rows and marginals along columns
#' @param std_mat : matrix of Gaussian sds, components along rows and marginals along columns
#' @param rho : equicorrealtion coefficient, lying betwwen -1 and 1
#'
#' @return : A list containing the log density and derivatives wrt mu, std and rho
#' @export
#'
#' @examples log_var_density(z = c(0, 0), mu_mat = matrix(0, ncol = 4, nrow = 2),
#' std_mat = matrix(rexp(8), ncol = 4, nrow = 2), rho = 0)
log_var_density = function(z, mu_mat, std_mat, rho) {
  if(sum(dim(mu_mat) == dim(std_mat)) != 2) {
    stop("matrices mu_mat and std_mat should have the same dimensions")
  }
  if(length(z) != nrow(mu_mat)) {
    stop("Argument z should have same length as row size of mu_mat")
  }
  if((rho > 1) | (rho < -1)) {
    stop("rho must lie between -1 and 1")
  }
  p = nrow(mu_mat)
  K = ncol(mu_mat)

  ## this defines some useful quantities to be used later
  trans_values = function(z, mu_mat, std_mat) {
    p = nrow(mu_mat)
    K = ncol(mu_mat)
    z_mat = matrix(z, nrow = p, ncol = K)
    arg_mat = (z_mat - mu_mat)/std_mat

    f_vec = (rowSums(stats::dnorm(arg_mat)))/K
    F_vec = (rowSums(stats::pnorm(arg_mat)))/K
    Y = stats::qnorm(F_vec)
    return(list(arg_mat = arg_mat, f_vec = f_vec, F_vec = F_vec, Y = Y))
  }

  val = trans_values(z, mu_mat, std_mat)
  arg_mat = val$arg_mat
  f_vec = val$f_vec
  F_vec = val$F_vec
  Y = val$Y

  f_mat = matrix(f_vec, ncol = K, nrow = p)

  log_det = ((p - 1) * log(1 - rho)) + log(1 + ((p - 1) * rho))

  log_marginals = log(f_vec)
  Y_term = (as.numeric(crossprod(Y)) - ((sum(Y)^2)/(1 + ((p - 1) * rho))))
  quad_term = (rho/(1 - rho)) * Y_term

  # returning log of density
  log_density = (-0.5 * (log_det + quad_term)) + sum(log_marginals)

  const_vec = as.vector(1/(stats::dnorm(Y)))
  const_vec = const_vec * ((rho / (rho - 1)) * Y_term) / K
  const_mat = matrix(const_vec, ncol = K, nrow = p)
  del_1 = const_mat * stats::dnorm(arg_mat) / std_mat
  del_2 = (arg_mat * stats::dnorm(arg_mat) / (std_mat * f_mat))/K

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






