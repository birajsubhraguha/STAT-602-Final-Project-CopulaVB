# this function takes as input a theta vector of size p(p - 1)/2 and outputs a correlation matrix
cor_mat = function(theta) {
  size_theta = length(theta)

  # checking compatibility of length of input
  temp = sqrt(1 + (8 * size_theta))
  check = (temp - floor(temp) != 0)
  if(check) {
    stop("input theta vector must be of size p(p - 1)/2 for integer p greater than 1")
  }

  size_cor = as.integer((1 + temp)/2)

  # checking compatibility of range of thetas
  if((min(theta) < -pi) | (max(theta) > pi)) {
    stop("all thetas must lie between -pi and pi")
  }

  # defining starting index, starts from second l = 2 !!!
  l = 2:size_cor
  start_ind = (l^2 - 3*l + 4)/2

  # defining the cos and sin vectors
  cos_vec = cos(theta)
  sin_vec = sin(theta)

  # creating the cos and sin matrices, NOT SPACE EFFICIENT, as almost half the entries are zero !!
  cos_mat = sin_mat = matrix(c(1, rep(0, size_cor - 1)), nrow = 1)
  for (i in 1:start_ind) {
    k = start_ind[i]
    cos_mat = rbind(cos_mat1, c(1, cumprod(cos_vec[k : (k + i - 1)]), rep(0, size_cor - i - 1)))
    sin_mat = rbind(sin_mat, c(sin_vec[k : k + i - 1]), 1, rep(0, size_cor - i - 1))
  }

  # creating the Cholesky lower triangular matrix
  lower_chol = cos_mat * sin_mat

  # creating the correaltion matrix
  return(tcrossprod(lower_chol))
}


