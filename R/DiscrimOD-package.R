#' Finding Optimal Discrimination Designs by Hybirdizing Particle Swarm and L-BFGS Algorithms
#'
#' This package adopts a hybrid algorithm
#' to search for optimal discrimination designs
#' when there are two or more than two competing models
#' under normal or non-normal error assumption.
#' This hybrid algorithm is chosen to efficiently solve the maximin design criteria
#' in the optimal discrimination design problem which is usually a challenging task.
#' It combines the particle swarm optimization (PSO) algorithm and the L-BFGS algorithm
#' to tackle the outer and inner objectives of the maximin design criterion, respectively.
#' The equivalence theorems for various discriminaiton criteria are also available for
#' verifying the optimal discrimination designs.
#'
#' @details
#' For the case of two competing models with normal errors, the package searches
#' for T-optimal design introduced in Atkinson and Fedorov (1975).
#' That is, we assume one of the competing models to be the true model and
#' the fixed nominal values of the parameters.
#' The another model is called the rival model and its parameter are assumed to be in
#' a specifued parameter space.
#' The T-optimal design maximizes the minimal squared distance between
#' the true and the rival models among the parameter space of the latter.
#'
#' If the errors are non-normal, the package searches for the KL-optimal design
#' based on Lopez-Fidalgo et al. (2007).
#' The approach is similar except that the distance measure for non-normal models
#' is the Kullback-Leibler (KL) divergence.
#' This package has equipped some commonly used KL-divergence functions for convenient uses
#' and it also allows users to customize the distance measure function by themselves.
#'
#' When there are more than two competing models, this package searches for the
#' max-min optimal discrimination design in Tommasi et al. (2016).
#' To find the max-min optimal design, we need a two-stpe approach.
#' Similarly, we assume one true model, fix its nominal values and treat the remaining models as rival models.
#' The first step is to identify the T/KL-optimal designs for each pair of true model and one rival model.
#' The second step searches for the max-min optimal discrimination design by maximizing
#' the minimal discriminaiton efficiency among the efficiencies relative to each T/KL-optimal design.
#'
#' @references Our Discrimination design paper.
#' @references Atkinson, A. C. and Fedorov, V. V. (1975). The design of experiments for discriminating between two rival models. Biometrika, 62(1):57-70.
#' @references Lopez-Fidalgo, J., Tommasi, C., and Trandafir, P. C. (2007). An optimal experimental design criterion for discriminating between non-normal models. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 69(2):231-242.
#' @references Tommasi, C., Martin-Martin, R., and Lopez-Fidalgo, J. (2016). Max-min optimal discriminating designs for several statistical models. Statistics and Computing, 26(6):1163-1172.
#' @references Eberhart, R. C. and Kennedy, J. (1995). A new optimizer using particle swarm theory. In Proceedings of the sixth international symposium on micro machine and human science, pages 39-43. IEEE.
#' @references Nocedal, J. and Wright, S. (2006). Numerical Optimization. Springer.
#'
#' @docType package
#' @name DiscrimOD-package
#' @useDynLib DiscrimOD
#' @importFrom Rcpp cppFunction sourceCpp
#' @importFrom inline cxxfunction
#' @importFrom utils globalVariables
#' @importFrom stats optim lm glm
NULL
