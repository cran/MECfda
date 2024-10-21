#' @title Use UP_MEM or MP_MEM substitution to apply
#' (generalized) linear regression with one functional covariate with measurement error.
#' @description
#' The Mixed-effect model (MEM) approach is a two-stage-based method that employs functional mixed-effects models.
#' It allows us to delve into the nonlinear measurement error model,
#' where the relationship between the true and observed measurements is not constrained to be linear,
#' and the distribution assumption on the observed measurement is relaxed to encompass
#' the exponential family rather than being limited to the Gaussian distribution.
#' The MEM approach employs point-wise (UP_MEM) and multi-point-wise (MP_MEM) estimation
#' procedures to avoid potential computational complexities caused by analyses of
#' multi-level functional data and computations of potentially intractable and complex integrals.
#' @references Luan, Yuanyuan, et al. "Scalable regression calibration approaches to
#' correcting measurement error in multi-level generalized functional linear regression models
#' with heteroscedastic measurement errors." arXiv preprint arXiv:2305.12624 (2023).
#' @param data.Y Response variable, can be an atomic vector, a one-column matrix or data frame,
#' recommended form is a one-column data frame with column name.
#' @param data.W A 3-dimensional array, represents \eqn{W}, the measurement of \eqn{X}.
#' Each row represents a subject.
#' Each column represent a measurement (time) point.
#' Each layer represents an observation.
#' @param data.Z Scalar covariate(s), can be not input or NULL (when there's no scalar covariate),
#' an atomic vector (when only one scalar covariate), a matrix or data frame,
#' recommended form is a data frame with column name(s).
#' @param method The method to construct the substitution \eqn{X}.
#' Available options: 'UP_MEM', 'MP_MEM', 'average'.
#' @param t_interval A 2-element vector, represents an interval,
#' means the domain of the functional covariate. Default is c(0,1), represent interval \eqn{[0,1]}.
#' @param t_points Sequence of the measurement (time) points, default is NULL.
#' @param d The number of time points involved for MP_MEM (default and miniumn is 3).
#' @param family.W Distribution of \eqn{W} given \eqn{X},  Available options: "gaussian","poisson".
#' @param family.Y A description of the error distribution
#' and link function to be used in the model, see \code{\link[stats]{family}}.
#' @param formula.Z A formula without the response variable,
#' contains only scalar covariate(s), use the format of lme4 package if random effects exist.
#' e.g. ~ Z_1 + (1|Z_2).
#' If not assigned, include all scalar covariates and intercept term as fixed effects.
#' @param basis.type Type of function basis.
#' Can only be assigned as one type even if there is more than one functional covariates.
#' Available options: 'Fourier' or 'Bspline',
#' represent Fourier basis and b-spline basis respectively.
#' For the detailed form for Fourier and b-splines basis,
#' see \code{\link{fourier_basis_expansion}} and \code{\link{bspline_basis_expansion}}.
#' @param basis.order Indicate number of the function basis.
#' When using Fourier basis \eqn{\frac{1}{2},\sin k t, \cos k t, k = 1,\dots,K},
#' basis.order is the number \eqn{K}.
#' When using b-splines basis \eqn{\{B_{i,p}(x)\}_{i=-p}^{k}}, basis.order is the number of splines,
#' equal to \eqn{k+p+1}. (same as arguement \code{df} in \code{\link[splines]{bs}}.)
#' May set a individual number for each functional covariate.
#' When the element of this argument is less than the number of functional covariates, it will be used recursively.
#' @param bs_degree Degree of the piecewise polynomials if use b-splines basis, default is 3.
#' See \code{degree} in \code{\link[splines]{bs}}.
#' @param smooth Whether to smooth the substitution of \eqn{X}. Default is FALSE.
#' @param silent Whether not to show the state of the running of the function. Default is TRUE.
#'
#' @return Returns a \code{fcRegression} object. See \code{\link{fcRegression}}.
#' @importFrom methods hasArg
#' @importFrom magrittr `%>%`
#' @importFrom nlme lme
#' @importFrom nlme pdCompSymm
#' @importFrom glme glme
#' @importFrom mgcv gam
#' @importFrom lme4 lmer
#' @importFrom lme4 glmer
#' @import stats
#' @export
#' @examples
#' data(MECfda.data.sim.0.1)
#' res = ME.fcRegression_MEM(data.Y = MECfda.data.sim.0.1$Y,
#'                           data.W = MECfda.data.sim.0.1$W,
#'                           data.Z = MECfda.data.sim.0.1$Z,
#'                           method = 'UP_MEM',
#'                           family.W = "gaussian",
#'                           basis.type = 'Bspline')
#'
ME.fcRegression_MEM = function(data.Y, data.W, data.Z,
                               method = c('UP_MEM', 'MP_MEM', 'average'),
                               t_interval = c(0,1), t_points = NULL, d = 3,
                               family.W = c("gaussian","poisson"), family.Y = "gaussian", formula.Z,
                               basis.type = c('Fourier','Bspline'),basis.order = NULL, bs_degree = 3,
                               smooth=FALSE, silent=TRUE){
  cov.model = "us"
  n   = dim(data.W)[1]
  t   = dim(data.W)[2]
  m_w = dim(data.W)[3]
  if(is.null(basis.order)){basis.order = as.integer(ceiling(n^(1/4))+2)}
  if(is.null(t_points)){t_points = seq(t_interval[1], t_interval[2], length.out=t)}
  # if(!methods::hasArg(formula.Z)){
  #   formula.Z = NULL
  # }

  if(method=="average"){
    M_t.log = log(data.W+1)
    X.hat = apply(M_t.log, c(1,2), mean)
  }else if(method %in% c('UP_MEM', 'MP_MEM')){
    switch (method,
            'UP_MEM' = {
              fit= UP_MEM( model=family.W,smooth=smooth, data=data.W,silent = silent)
            },
            'MP_MEM' = {
              fit = MP_MEM( model=family.W, d=d,smooth=smooth, data=data.W,silent = silent)
            }
    )
    X.hat= matrix(unlist(fit), nrow= n, ncol = t)
  }else{stop("Unknown method selection")}
  X.hat = functional_variable(X=X.hat,t_0 = t_interval[1], period = t_interval[2]-t_interval[1], t_points = t_points)
  ret = fcRegression(Y = data.Y, FC = list(X.hat), Z = data.Z, formula.Z = formula.Z, family = family.Y, basis.type = basis.type, basis.order = basis.order, bs_degree = bs_degree)
  return(ret)
}


