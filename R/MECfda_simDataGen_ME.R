#' Simulation Data Generation: Measurement Error Bias Correction of Scalar-on-function Regression
#' @description Generate data set for measurement error bias correction methods
#' for scalar-on-function regression in package MECfda
#' @param N Sample size.
#' @param J_W Number of repeated measurement (period), if applicable.
#' @param forwhich For which method of measurement error bias correction method
#' the data set is generated.
#' There are two options: 'MEM','IV','CLS','IV.SIMEX'.
#' @param t_interval A 2-element vector, represents an interval,
#' means the domain of the functional covariate.
#' Default is \code{c(0,1)}, represent interval \eqn{[0,1]}.
#' @param n_t Number of measurement time points.
#' @param seed Pseudo-random number generation seed.
#'
#' @returns return a list that possibly contains following elements.
#'    \item{Y}{An atomic vector of response variable}
#'    \item{Z}{A dataframe with a binary and a continuous scalar-valued covariate.}
#'    \item{W}{Observed values of function-valued covariate.}
#'    \item{M}{Instrumental vairable.}
#'    \item{t_interval}{Same as in the input argument.}
#'    \item{t_points}{Sequence of the measurement (time) points.}
#' @export
#' @examples
#' for (i in 1:4) {MECfda_simDataGen_ME(forwhich = c('MEM','IV','CLS','IV.SIMEX')[i])}
#' @importFrom MASS mvrnorm
MECfda_simDataGen_ME = function(N = 100, J_W = 7, forwhich = c('MEM','IV','CLS','IV.SIMEX'),
                                t_interval = c(0,1),n_t=24,seed = 0){
  set.seed(seed)
  Stationary_Gaussian_Process_Generation = function(time_points, ker_fun, mu = 0, Sigma_matrix = NULL){
    n = length(time_points)
    if(is.function(mu)){
      mu = mu(time_points)
    } else if (length(mu) == 1){
      mu = rep(mu,n)
    }
    if(is.null(Sigma_matrix)){
      ones_vector = as.matrix(rep(1,length(time_points)))
      time_points = as.matrix(time_points)
      m1 = time_points %*% t(ones_vector)
      diff_matrix = m1 - t(m1)
      Sigma_matrix = ker_fun(diff_matrix)
    }
    MASS::mvrnorm(1,mu,Sigma_matrix)
  }
  N = N
  J = J_W
  n_t = n_t
  t_0 = t_interval[1]
  period = t_interval[2] - t_interval[1]
  # t_points = t_0 + period*(1:n_t - 0.5)/n_t
  t_points = seq(0,1,length.out = n_t)
  ker = function(r, p = c(1,1)) p[1] * exp(- (r/p[2])^2 /2)
  Betafunc = function(t, met){
    switch(met,
           f1 = sin(2*pi*(t-t_0)/period),
    )
  }
  gama = c(gamma_0 = 5,
           gamma_1 = 0.2,
           gamma_2 = 0.4)

  sigma_0 = 0.02

  ### X,Y,Z
  X = exp(t(sapply(rep(1,N), function(xx)Stationary_Gaussian_Process_Generation(t_points,function(r) ker(r, c(0.8,0.2)),1))))
  Z = cbind(Z_1 = rnorm(N, 0, 1),Z_2 = rbinom(N, 1,0.6))
  Y = gama['gamma_0'] + c( period*crossprod(t(X),Betafunc(t_points,"f1"))/n_t + gama['gamma_1']*Z[,1] + gama['gamma_2']*Z[,2]+
          rnorm(N, 0,sigma_0))
  # X = functional_variable(X=X,t_0=0,period=1)
  Y = as.data.frame(Y)
  Z = as.data.frame(Z)

  ### W
  n_obs = J
  U = t(sapply(rep(1,n_obs*N), function(xx)Stationary_Gaussian_Process_Generation(t_points,function(r) ker(r, c(0.2,0.1)),0)))
  W = array(NA,dim = c(N,n_t,n_obs))
  for (obs in 1:n_obs) {
    W[,,obs] = X + U[N*(obs-1) + (1:N),]
  }
  rm(obs)

  ### M

  omega = t(sapply(rep(1,N), function(xx)Stationary_Gaussian_Process_Generation(t_points,function(r) ker(r, c(0.25,0.05)),0)))

  # delta.t = exp(Stationary_Gaussian_Process_Generation(t_points,function(r) ker(r, c(0.2,0.15)),-0.5))
  delta.t = 0.6
  M = t(t(X) * delta.t)[1,] + omega


  switch (forwhich,
          'MEM' =      { return(list(Y = Y, Z = Z, W = W,             t_interval = t_interval, t_points = t_points)) },
          'CLS' =      { return(list(Y = Y, Z = Z, W = W,             t_interval = t_interval, t_points = t_points)) },
          'IV.SIMEX' = { return(list(Y = Y, Z = Z, W = W[,,1], M = M, t_interval = t_interval, t_points = t_points)) },
          'IV' =  { Y = as.data.frame(1.6 +c( crossprod(t(X),Betafunc(t_points,"f1"))/n_t + rnorm(N, 0,sigma_0)))
                         return(list(Y = Y,        W = W[,,1], M = M, t_interval = t_interval, t_points = t_points )) }
  )
}

