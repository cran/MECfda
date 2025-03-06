#' @title Get MEM substitution for
#' (generalized) linear regression with one functional covariate with measurement error.
#' @description
#' The function to get the data of \eqn{\hat X_i(t)} using the mixed model based
#' measurement error bias correction method
#' proposed by Luan et al.
#' See \code{\link{ME.fcRegression_MEM}}
#' @references Luan, Yuanyuan, et al. "Scalable regression calibration approaches to
#' correcting measurement error in multi-level generalized functional linear regression models
#' with heteroscedastic measurement errors." arXiv preprint arXiv:2305.12624 (2023).
#' @param data.W A 3-dimensional array, represents \eqn{W}, the measurement of \eqn{X}.
#' Each row represents a subject.
#' Each column represent a measurement (time) point.
#' Each layer represents an observation.
#' @param method The method to construct the substitution \eqn{X}.
#' Available options: \code{'UP_MEM'}, \code{'MP_MEM'}, 'average'.
#' @param d The number of time points involved for MP_MEM (default and miniumn is 3).
#' @param family.W Distribution of \eqn{W} given \eqn{X},  Available options: \code{"gaussian"}, \code{"poisson"}.
#' @param smooth Whether to smooth the substitution of \eqn{X}. Default is \code{FALSE}.
#'
#' @returns A numeric value matrix of \eqn{\hat X_i(t)}.
#' @export
#' @examples
#' data(MECfda.data.sim.0.1)
#' X_hat = MEM_X_hat(data.W = MECfda.data.sim.0.1$W,
#'                   method = 'UP_MEM',
#'                   family.W = "gaussian")
#'
MEM_X_hat = function(data.W,
                     method = c('UP_MEM', 'MP_MEM', 'average'),
                     # t_interval = c(0,1), t_points = NULL,
                     d = 3,
                     family.W = c("gaussian","poisson"),
                     smooth=FALSE){
  # cov.model = "us"
  n   = dim(data.W)[1]
  t   = dim(data.W)[2]
  m_w = dim(data.W)[3]
  # if(is.null(t_points)){t_points = seq(t_interval[1], t_interval[2], length.out=t)}
  if(method=="average"){
    M_t.log = log(data.W+1)
    X.hat = apply(M_t.log, c(1,2), mean)
  }else if(method %in% c('UP_MEM', 'MP_MEM')){
    switch (method,
            'UP_MEM' = {
              fit= UP_MEM( model=family.W,smooth=smooth, data=data.W,silent = TRUE)
            },
            'MP_MEM' = {
              fit = MP_MEM( model=family.W, d=d,smooth=smooth, data=data.W,silent = TRUE)
            }
    )
    X.hat= matrix(unlist(fit), nrow= n, ncol = t)
  }else{stop("Unknown method selection")}
  # X.hat = functional_variable(X=X.hat,t_0 = t_interval[1], period = t_interval[2]-t_interval[1], t_points = t_points)
  return(X.hat)
}




