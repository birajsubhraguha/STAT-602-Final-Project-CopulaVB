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
The package currently contains functionality to implement Copula based Variational Bayes methods. The variational family used has finite Gaussian mixtures (equal weight, for simplicity) as marginals. It is equipped with a Gaussian Copula, whose correlation matrix is an equicorrelation matrix. Thus the set of variational parameters include all the Gaussian component means and standard deviations, and one single $\rho$ as the correlation parameter. The package helps to generate from this copula-defined joint variational density, which is a bettered version of mean-field approximations commonly found in literature of Variational Bayes procedures. It helps to capture, to some extent, the co-dependence of the posterior variables. The package currently provides the functionality to find the variational derivatives, that is almost universally used for any Variational algorithm one can choose to employ.
# Anticipated Updates
# Functions available
## generate_var_density
## log_var_density
# Tests

