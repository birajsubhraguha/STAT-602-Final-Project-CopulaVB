require(MASS)
require(stats)

## INPUT
# n_sample = no. of samples to be drawn
# mu_mat, std_mat = matrix of the Gaussian means/standard deviations of size p x K
# p = dimension of parameter space, K = number of mixture components, rho = equicorrelation coefficient
## OUTPUT
# n_samples from the variational density

#' generate_var_density
#'
#' \code{generate_var_density} Generates samples from Gaussian Copula model with marginals specified by finite Gaussian mixtures and having equicorrelation coefficient rho
#'
#' @param n_sample : number of samples
#' @param mu_mat : p x K - dimensional matrix of Gaussian means, components along rows and marginals along columns
#' @param std_mat : p x K - dimensional matrix of Gaussian sds, components along rows and marginals along columns
#' @param rho : scalar, equicorrealtion coefficient, lying betwwen -1 and 1
#' @param seed : seed
#'
#' @return n_sample x p - dimensional matrix, rows are samples from the variational density
#'
#' @export
#'
#' @examples
#' ## generating 2D data from the variational density
#' out = generate_var_density(n_sample = 10000, mu_mat = matrix(1:8, nrow = 2, byrow = TRUE),
#' std_mat = matrix(1, ncol = 4, nrow = 2), rho = .5, seed = 123)
#' out = as.data.frame(out)
#' colnames(out) = c("x", "y")
#'
#' ## 2D count plot of above data
#' require(ggplot2)
#' plot = ggplot2::ggplot(out, aes(x=x, y=y) ) + geom_bin2d(bins = 70) +
#' scale_fill_continuous(type = "viridis") + theme_bw()
#' print(plot + ggtitle("2D count plot"))
#'
generate_var_density = function(n_sample, mu_mat, std_mat, rho, seed = 1234) {
  if(sum(dim(mu_mat) == dim(std_mat)) != 2) {
    stop("matrices mu_mat and std_mat should have the same dimensions")
  }
  p = nrow(mu_mat)
  K = ncol(mu_mat)

  ## defining mixture of Gaussians
  ## INPUT
  # z = argument at which cdf is evaluated, mu_vec = vector of means, std_vec = vector of standard deviations
  # vectors mu_vec and std_vec must be of same length
  ## OUTPUT
  # value of the mixture cdf at z
  mix_pnorm = function(z, mu_vec, std_vec) {
    if(length(mu_vec) != length(std_vec)) {
      # checking compatibility
      stop("mu_vec and std_vec should be vectors of the same length")
    }

    K = length(mu_vec)
    out = numeric(0)
    out = (sum(stats::pnorm(z, mu_vec, std_vec)))/K
    return(out)
  }

  ## inverting cdf of mixture of Gaussians
  ## INPUT
  # a_vec = vector of arguments at which inverse-cdf is evaluated, mu_vec = vector of means, std_vec = vector of standard deviations
  # members of a_vec must lie between 0 and 1
  # vectors mu_vec and std_vec must be of same length
  ## OUTPUT
  # a_vec length vector of values of the inverse of the mixture cdf
  mix_qnorm = function(a_vec, mu_vec, std_vec) {
    # checking compatibility
    if(length(mu_vec) != length(std_vec)) {
      stop("mu_vec and std_vec should be vectors of the same length")
    }
    if((max(a_vec) > 1) | (min(a_vec) < 0)) {
      stop("all members of vector a must lie between 0 and 1")
    }

    out = numeric(0)
    for(i in 1:length(a_vec)) {
      temp_mix_pnorm = function(y) {
        out = mix_pnorm(y, mu_vec, std_vec) - a_vec[i]
        return(out)
      }
      val = stats::uniroot(temp_mix_pnorm, interval = c(-100, 100))
      out = c(out, val$root)
    }
    return(out)
  }

  ## The following codes the drawing of samples from the variational density
  # depends on pacakge "MASS"

  ## defining equicorrelation matrix with correlation rho
  ## INPUT
  # rho = scalar lying between -1 and 1, p = size of the equicorrelation matrix
  ## OUTPUT
  # equicorrelation matrix of size p and correlation rho
  rho_mat = function(rho, p) {
    if((rho > 1) | (rho < -1)) {
      stop("rho must be a scalar lying between -1 and 1")
    }
    out = ((1 - rho) * diag(p)) + (rho * tcrossprod(rep(1, p)))
    return(out)
  }

  # generating mean zero multivariate normal with covariance rho_mat
  set.seed(seed)
  mvt = MASS::mvrnorm(n_sample, mu = rep(0, p), Sigma = rho_mat(rho, p))

  # transforming marginally to uniforms
  uni_mvt = stats::pnorm(mvt)

  # applying inverse mixture cdf transform to get samples
  samples = numeric(0)
  for(i in 1:p) {
    temp = mix_qnorm(a_vec = uni_mvt[, i], mu_vec = mu_mat[i, ], std_vec = std_mat[i, ])
    samples = cbind(samples, temp)
  }
  return(samples)
}
