#' @title Solve quantile regression models with functional covariate(s).
#' @description
#' Fit a quantile regression models below
#' \deqn{Q_{Y_i|X_i,Z_i}(\tau) = \sum_{l=1}^L\int_\Omega \beta_l(\tau,t) X_{li}(t) dt + (1,Z_i^T)\gamma}
#' where \eqn{Q_{Y_i}(\tau) = F_{Y_i|X_i,Z_i}^{-1}(\tau)} is the
#' \eqn{\tau}-th quantile of \eqn{Y_i} given \eqn{X_i(t)} and \eqn{Z_i},
#' \eqn{\tau\in(0,1)}.
#' \cr Model allows one or multiple functional covariate(s) as fixed effect(s),
#' and zero, one, or multiple scalar-valued covariate(s).
#'
#'
#' @param Y Response variable, can be an atomic vector, a one-column matrix or data frame,
#' recommended form is a one-column data frame with column name
#' @param FC Functional covariate(s),
#' can be a "functional_variable" object or a matrix or a data frame or a list of these object(s)
#' @param Z Scalar covariate(s), can be NULL or not input (when there's no scalar covariate),
#' an atomic vector (when only one scalar covariate), a matrix or data frame,
#' recommended form is a data frame with column name(s)
#' @param formula.Z A formula without the response variable, contains only scalar covariate(s).
#' If not assigned, include all scalar covariates and intercept term.
#' @param tau Quantile \eqn{\tau\in(0,1)}, default is 0.5. See \code{\link[quantreg]{rq}}.
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
#' @param bs_degree Degree of the piecewise polynomials if use b-splines basis, default is 3.
#' See \code{degree} in \code{\link[splines]{bs}}.
#'
#'
#' @return fcQR returns an object of class "fcQR".
#'     It is a list that contains the following elements.
#'     \item{regression_result}{Result of the regression.}
#'     \item{FC.BasisCoefficient}{A list of Fourier_series or bspline_series object(s),
#'     represents the functional linear coefficient(s) of the functional covariates.}
#'     \item{function.basis.type}{Type of funtion basis used.}
#'     \item{basis.order}{Same as in the arguemnets.}
#'     \item{data}{Original data.}
#'     \item{bs_degree}{Degree of the splines, returned only if b-splines basis is used.}
#' @export
#' @author Heyang Ji
#' @examples
#' data(MECfda.data.sim.0.0)
#' res = fcQR(FC = MECfda.data.sim.0.0$FC, Y=MECfda.data.sim.0.0$Y, Z=MECfda.data.sim.0.0$Z,
#'            basis.order = 5, basis.type = c('Bspline'))
#'
#' @importFrom methods hasArg
#' @importFrom stringr str_replace
#' @importFrom stringr str_detect
#' @importFrom quantreg rq
#' @import stats
#'
fcQR = function(Y, FC, Z, formula.Z, tau = 0.5, basis.type = c('Fourier','Bspline'), basis.order = 6L, bs_degree = 3){

  if(!methods::hasArg(Z)){
    Z = NULL
  }
  if(!methods::hasArg(formula.Z)){
    formula.Z = NULL
  }

  if(is.matrix(FC) | is.data.frame(FC)){
    FC = functional_variable(as.matrix(FC))
  }
  if(class(FC) %in% "functional_variable"){
    fc_list = list(FC)
  }else{
    fc_list = FC
  }
  rm(FC)

  if(any(basis.order != as.integer(basis.order))){
    stop("basis.order should be integer(s)")
  }
  basis.order = as.integer(basis.order)
  if(length(basis.order)!=length(fc_list)){basis.order = rep(basis.order,length(fc_list))[1:length(fc_list)]}

  Y = as.data.frame(Y)
  if(is.null(colnames(Y))) colnames(Y) = 'Y'
  if(!is.null(Z)){
    Z = as.data.frame(Z)
    if(is.null(colnames(Z))) colnames(Z) = paste('Z',1:ncol(Z),sep = '_')
    if(any(
      any(nrow(Y) != data.frame(lapply(fc_list, dim))['subject',]),
      nrow(Y) != nrow(Z)
    )){
      stop("dimensionality of variables doesn't match")
    }
  }else{
    Z = as.data.frame(matrix(nrow = nrow(Y),ncol = 0))
  }

  if(is.null(names(fc_list))){
    if(length(fc_list) == 1){
      names(fc_list) = 'Functional_Covariate'
    }else{
      names(fc_list) = paste0('FunctionalCovariate',1:length(fc_list))
    }
  }


  {
    switch (basis.type,
            'Fourier' = {
              BE = NULL
              for (i in 1:length(fc_list)) {
                X = fc_list[[i]]
                n_k = basis.order[i]
                BE_X = fourier_basis_expansion(X,n_k)
                colnames(BE_X) = paste(names(fc_list)[i],colnames(BE_X),sep = '.')
                BE = cbind(BE,BE_X)
              }
            },
            'Bspline' = {
              BE = NULL
              for (i in 1:length(fc_list)) {
                X = fc_list[[i]]
                n_k = basis.order[i]
                BE_X = bspline_basis_expansion(X,n_k,bs_degree)
                colnames(BE_X) = paste(names(fc_list)[i],colnames(BE_X),sep = '.')
                BE = cbind(BE,BE_X)
              }
            }
    )
    rm(X,i)
    data.funRegress = as.data.frame(cbind(Y,BE,Z))
    rm(BE_X)
  }


  if(is.null(formula.Z)){
    rm(BE)
    fmla = as.formula(paste0(colnames(Y),' ~ .'))
  }else{
    formula.Z = format(formula.Z)
    formula.Z = str_replace(formula.Z,'~','')
    fmla = as.formula(paste0(colnames(Y),' ~ ',
                             paste0(colnames(BE),collapse = ' + '), ' + ',
                             formula.Z
    ))
    rm(BE)
  }
  res.funRegress = rq(formula = fmla, data = data.funRegress, tau = tau)

  res.coef = coef(res.funRegress)


  {

    FP_basisCoef = list()
    i = 0
    for (fv_name in names(fc_list)) {
      i = i + 1
      n_k = basis.order[i]
      switch (basis.type,
              'Fourier' = {
                FP_basisCoef = append(FP_basisCoef,
                                      Fourier_series(
                                        double_constant = res.coef[paste(fv_name,'a_0',sep = '.')],
                                        sin = res.coef[str_detect(names(res.coef),fv_name) & str_detect(names(res.coef),'sin')],
                                        cos = res.coef[str_detect(names(res.coef),fv_name) & str_detect(names(res.coef),'cos')],
                                        k_sin = 1:n_k,
                                        k_cos = 1:n_k,
                                        t_0 = fc_list[[i]]@t_0,
                                        period = fc_list[[i]]@period
                                      ))
              },
              'Bspline' = {
                FP_basisCoef = append(FP_basisCoef,
                                      bspline_series(coef = c(res.coef[str_detect(names(res.coef),fv_name) & str_detect(names(res.coef),'bs')]),
                                                     bspline_basis = bspline_basis(
                                                       Boundary.knots = c(fc_list[[i]]@t_0,fc_list[[i]]@t_0 + fc_list[[i]]@period),


                                                       df             = n_k,
                                                       degree         = bs_degree
                                                     )
                                      ))
              }
      )

    }
    rm(i)
    names(FP_basisCoef) = names(fc_list)
  }

  data.origin = list(Y = Y, fc_list = fc_list, Z = Z)

  ret = list(regression_result = res.funRegress,
             FC.BasisCoefficient = FP_basisCoef,
             function.basis.type = basis.type,
             basis.order = basis.order,
             data = data.origin
  )
  if(basis.type == 'Bspline'){ret$bs_degree = bs_degree}
  class(ret) = c('fcQR',paste0(basis.type,'_basis'),class(res.funRegress))
  return(ret)

}
