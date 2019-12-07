# STAT-602-Final-Project-CopulaVB
Copula based Variational Bayes
# Package Name
**CopulaVB**
# Installation Instructions
1) Go to https://github.com/birajsubhraguha/STAT-602-Final-Project-CopulaVB and download the .tar.gz file.
2) Open RStudio, go to the 'Packages' pane and select install.
3) Select source of the package as .tar.gz type of files and then select the downloaded package from 'Browse'.
4) Hit 'Install'.
5) Use **library(CopulaVB)** in the RStudio console to use the package.
# General Description
The package currently contains initial functionality to implement Copula based Variational Bayes methods. The variational family used has finite Gaussian mixtures (equal weight, for simplicity) as marginals. It is equipped with a Gaussian Copula, whose correlation matrix is an equicorrelation matrix. Thus the set of variational parameters include all the Gaussian component means and standard deviations, and one single rho as the correlation parameter. The package helps to generate from this copula-defined joint variational density, which is a bettered version of mean-field approximations commonly found in literature of Variational Bayes procedures. It helps to capture, to some extent, the co-dependence of the posterior variables. The package currently provides the functionality to find the variational derivatives, that is almost universally used for any Variational algorithm one can choose to employ.
# Anticipated Updates
I aim to add the basic Black-Box Variational algorithm in future, as well as model some examples from where the log-joint density will be used to approximate the posterior.
# Functions available
## generate_var_density
This function draws samples from the specified variational density. These draws are essential in algorithms like the Black-Box Variational algorithm.
## log_var_density
This function calculates log-density as well as the variational derivatives.
# Tests
Two tests for each function were added:
## tests for generate_var_density
The first test checked whether the marginal means are correct, while the second checked whether setting the equicorrelation coefficient to zero produces independent samples. The tests use large samples, hence invoke Strong Law of Large Numbers. 
## tests for log_var_density
The first test checked for derivaties vanishing under zero correlation and zero means. The second test checked for density value at a specified point.
