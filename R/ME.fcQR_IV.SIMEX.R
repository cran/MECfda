#' @title Bias correction method of applying quantile linear regression to dataset
#' with one functional covariate with measurement error using instrumental variable.
#' @description
#' Perform a two-stage strategy to correct the measurement error of
#' a function-valued covariate and then fit a linear quantile regression model.
#' In the first stage, an instrumental variable is used to estimate
#' the covariance matrix associated with the measurement error.
#' In the second stage, simulation extrapolation (SIMEX) is used to correct for
#' measurement error in the function-valued covariate. \cr
#' See detailed model in the reference.
#' @references Tekwe, Carmen D., et al.
#' "Estimation of sparse functional quantile regression with measurement error: a SIMEX approach."
#' Biostatistics 23.4 (2022): 1218-1241.
#'
#' @param data.Y Response variable, can be an atomic vector, a one-column matrix or data frame,
#' recommended form is a one-column data frame with column name.
#' @param data.W A dataframe or matrix, represents \eqn{W}, the measurement of \eqn{X}.
#' Each row represents a subject. Each column represent a measurement (time) point.
#' @param data.Z Scalar covariate(s),
#' can be not input or \code{NULL} (when there's no scalar covariate),
#' an atomic vector (when only one scalar covariate), a matrix or data frame,
#' recommended form is a data frame with column name(s).
#' @param data.M A dataframe or matrix, represents \eqn{M},
#' the instrumental variable. Each row represents a subject.
#' Each column represent a measurement (time) point.
#' @param tau Quantile \eqn{\tau\in(0,1)}, default is 0.5.
#' @param t_interval A 2-element vector, represents an interval,
#' means the domain of the functional covariate. Default is c(0,1), represent interval \eqn{[0,1]}.
#' @param t_points Sequence of the measurement (time) points, default is \code{NULL}.
#' @param formula.Z A formula without the response variable,
#' contains only scalar covariate(s), no random effects. If not assigned,
#' include all scalar covariates and intercept term.
#' @param basis.type Type of funtion basis.
#' Can only be assigned as one type even if there is more than one functional covariates.
#' Available options: 'Fourier' or 'Bspline', represent Fourier basis and b-spline basis respectively.
#' For the detailed form for Fourier and b-splines basis,
#' see \code{\link{fourier_basis_expansion}} and \code{\link{bspline_basis_expansion}}.
#' @param basis.order Indicate number of the function basis.
#' When using Fourier basis \eqn{\frac{1}{2},\sin k t, \cos k t, k = 1,\dots,K},
#' basis.order is the number \eqn{K}.
#' When using b-splines basis \eqn{\{B_{i,p}(x)\}_{i=-p}^{k}},
#' basis.order is the number of splines, equal to \eqn{k+p+1}.
#' (same as arguement \code{df} in \code{\link[splines]{bs}}.)
#' May set a individual number for each functional covariate.
#' When the element of this argument is less than the number of functional covariates,
#' it will be used recursively.
#' @param bs_degree Degree of the piecewise polynomials if use b-splines basis,
#' default is 3. See \code{degree} in \code{\link[splines]{bs}}.
#'
#' @return Returns a ME.fcQR_IV.SIMEX class object.
#' It is a list that contains the following elements.
#' \item{coef.X}{A Fourier_series or bspline_series object,
#' represents the functional coefficient parameter of the functional covariate.}
#' \item{coef.Z}{The estimate of the linear coefficients of the scalar covariates.}
#' \item{coef.all}{Original estimate of linear coefficients.}
#' \item{function.basis.type}{Type of funtion basis used.}
#' \item{basis.order}{Same as in the input arguements.}
#' \item{t_interval}{A 2-element vector, represents an interval,
#' means the domain of the functional covariate.}
#' \item{t_points}{Sequence of the measurement (time) points.}
#' \item{formula}{Regression model.}
#' \item{formula.Z}{formula object contains only the scalar covariate(s).}
#' \item{zlevels}{levels of the non-continuous scalar covariate(s). }
#' @export
#' @examples
#' data(MECfda.data.sim.0.2)
#' \donttest{
#' res = ME.fcQR_IV.SIMEX(data.Y = MECfda.data.sim.0.2$Y,
#'                        data.W = MECfda.data.sim.0.2$W,
#'                        data.Z = MECfda.data.sim.0.2$Z,
#'                        data.M = MECfda.data.sim.0.2$M,
#'                        tau = 0.5,
#'                        basis.type = 'Bspline')
#' }
#' @importFrom corpcor make.positive.definite
#' @importFrom MASS mvrnorm
#' @importFrom quantreg rq
#' @importFrom Matrix forceSymmetric
#' @import stats
ME.fcQR_IV.SIMEX = function(data.Y, data.W, data.Z, data.M, tau = 0.5,
                            t_interval = c(0,1), t_points = NULL,
                            formula.Z,
                            basis.type = c('Fourier','Bspline'),basis.order = NULL, bs_degree = 3){
  data.Y = as.data.frame(data.Y)
  if(is.null(colnames(data.Y))) colnames(data.Y) = 'Y'

  if(!methods::hasArg(formula.Z)){
    formula.Z = NULL
  }
  if(!methods::hasArg(data.Z)){
    data.Z = NULL
  }
  if(!is.null(data.Z)){
    data.Z = as.data.frame(data.Z)
    if(is.null(colnames(data.Z))) colnames(data.Z) = paste('Z',1:ncol(data.Z),sep = '_')
    if(any(
      nrow(data.Y) != nrow(data.W),
      nrow(data.Y) != nrow(data.Z)
    )){
      stop("dimensionality of variables doesn't match")
    }
  }else{
    data.Z = as.data.frame(matrix(nrow = nrow(data.Y),ncol = 0))
  }

  n   = dim(data.W)[1]
  t   = dim(data.W)[2]
  if(is.null(basis.order)){basis.order = as.integer(ceiling(n^(1/4))+2)}
  if(is.null(t_points)){t_points = seq(t_interval[1], t_interval[2], length.out=t)}
  data.W = as.matrix(data.W)
  data.M = as.matrix(data.M)
  delta_t.hat = (colMeans(data.M))/(colMeans(data.W))
  M.new = data.M/matrix(rep(delta_t.hat,n), nrow=n, ncol=t, byrow = TRUE)
  data.W = functional_variable(data.W,t_0 = t_interval[1], period = t_interval[2] - t_interval[1], t_points = t_points)
  data.M = functional_variable(data.M,t_0 = t_interval[1], period = t_interval[2] - t_interval[1], t_points = t_points)
  M.new  = functional_variable(M.new ,t_0 = t_interval[1], period = t_interval[2] - t_interval[1], t_points = t_points)
  {
    switch (basis.type,
            'Fourier' = {
              W_i = fourier_basis_expansion(data.W,basis.order)
              M_i = fourier_basis_expansion(M.new, basis.order)
            },
            'Bspline' = {
              W_i = bspline_basis_expansion(data.W,basis.order,bs_degree)
              M_i = bspline_basis_expansion(M.new, basis.order,bs_degree)
            }
    )
    colnames(W_i) = paste('W_be',colnames(W_i),sep = '.')
    colnames(M_i) = paste('M_be',colnames(M_i),sep = '.')
  }
  if(is.null(formula.Z)){
    fmla = as.formula(paste0(colnames(data.Y),' ~ .'))
    scalar.covariate = colnames(data.Z)
  }else{
    scalar.covariate = all.vars(formula.Z)
    fmla = format(formula.Z)
    # fmla = stringr::str_replace(fmla,'~','')
    fmla = sub('~','',fmla)
    fmla = as.formula(paste0(colnames(data.Y),' ~ ',
                             paste0(colnames(W_i),collapse = ' + '), ' + ',
                             fmla
    ))
  }
  Sigma_xx = cov(W_i, M_i,use="complete.obs")
  Sigma_ww = var(W_i,use="complete.obs")
  Sigma_uu = Sigma_ww - Sigma_xx
  Sigma_uu = make.positive.definite(as.matrix(forceSymmetric(Sigma_uu)))
  B = 100
  lambda= seq(0.01,2.01,.05)
  gamma.simex = lapply(seq(1:B), function(b){
    ret = sapply(lambda, function(s) {
      # set.seed(b)
      U_b = mvrnorm(n, rep(0, ncol(W_i)), Sigma_uu, empirical = TRUE)
      W_lambda = W_i+ (sqrt(s)*U_b)
      data.funRegress = as.data.frame(cbind(data.Y,W_lambda,data.Z))
      model =  try(rq(formula = fmla, data = data.funRegress, tau = tau),TRUE)
      return(model[["coefficients"]])
    })
    colnames(ret) = lambda
    ret
  } )
  gamma_simex.ave = t(Reduce("+", gamma.simex)/B)
  res.coef = apply(gamma_simex.ave, 2, function(x){predict(lm(x~lambda + I(lambda^2)),newdata = data.frame(lambda = -1))})

  switch (basis.type,
          'Fourier' = {
            # col.X = stringr::str_detect(names(res.coef),'W_be') & (
            #   stringr::str_detect(names(res.coef),'sin') | stringr::str_detect(names(res.coef),'cos') | stringr::str_detect(names(res.coef),'a_0'))
            col.X = grepl('W_be',names(res.coef), fixed = TRUE) & (
              grepl('sin',names(res.coef), fixed = TRUE) | grepl('cos',names(res.coef), fixed = TRUE) | grepl('a_0',names(res.coef), fixed = TRUE))
            coef.X = Fourier_series(
              double_constant = res.coef[paste('W_be','a_0',sep = '.')],
              # sin = res.coef[stringr::str_detect(names(res.coef),'W_be') & stringr::str_detect(names(res.coef),'sin')],
              # cos = res.coef[stringr::str_detect(names(res.coef),'W_be') & stringr::str_detect(names(res.coef),'cos')],
              sin = res.coef[grepl('W_be',names(res.coef), fixed = TRUE) & grepl('sin',names(res.coef), fixed = TRUE)],
              cos = res.coef[grepl('W_be',names(res.coef), fixed = TRUE) & grepl('cos',names(res.coef), fixed = TRUE)],
              k_sin = 1:basis.order,
              k_cos = 1:basis.order,
              t_0 = t_interval[1],
              period = t_interval[2] - t_interval[1]
            )
          },
          'Bspline' = {
            # col.X = stringr::str_detect(names(res.coef),'W_be') & stringr::str_detect(names(res.coef),'bs')
            col.X = grepl('W_be',names(res.coef), fixed = TRUE) & grepl('bs',names(res.coef), fixed = TRUE)
            # coef.X = bspline_series(coef = c(res.coef[stringr::str_detect(names(res.coef),'W_be') & stringr::str_detect(names(res.coef),'bs')]),
            coef.X = bspline_series(coef = c(res.coef[grepl('W_be',names(res.coef), fixed = TRUE) & grepl('bs',names(res.coef), fixed = TRUE)]),
                                    bspline_basis = bspline_basis(
                                      Boundary.knots = c(t_interval[1],t_interval[2]),
                                      df             = basis.order,
                                      degree         = bs_degree
                                    )
            )
          }
  )
  coef.Z = res.coef[!col.X]
  data.Z.1 = data.Z[,scalar.covariate]
  zlevels = lapply(data.Z.1[,sapply(data.Z.1,is.factor),drop = FALSE],levels)
  ret = list(coef.X = coef.X, coef.Z = coef.Z, coef.all = res.coef,
             function.basis.type = basis.type,
             basis.order = basis.order,
             t_interval = t_interval,
             t_points = t_points,
             formula = fmla,
             formula.Z = formula.Z,
             zlevels = zlevels
  )
  if(basis.type == 'Bspline'){ret$bs_degree = bs_degree}
  class(ret) = 'ME.fcQR_IV.SIMEX'
  return(ret)
}
