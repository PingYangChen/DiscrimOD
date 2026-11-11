# File Name: Discrim_Install.r
# Description:
#   To use the **DiscrimOD** package, one need to have R packages 'remotes',
#   'Rcpp' and 'RcppArmadillo' installed in advance.  This file is written for
#   users to install all the required R packages and our **DiscrimOD** package.
#		Please copy the codes below and run them in R.
# ----------------------------------------------------------------------------

# Install required R packages
pkgRequired <- c("remotes", "Rcpp", "RcppArmadillo")
install.packages(pkgRequired)
# Install the **DiscrimOD** package from GitHub
remotes::install_github("PingYangChen/DiscrimOD")
# Try if it is installed successfully
library(DiscrimOD)
?DiscrimOD
