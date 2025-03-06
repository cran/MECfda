#' Simulation Data Generation: Scalar-on-function Regression
#' @description Generate data set for scalar-on-function regression
#' @param N Sample size.
#' @param distribution Conditional distribution of response varaible given the covariate
#' (\eqn{Y_i|X_i(t),Z_i}).
#' There are two options: \code{'Gaussian'} and \code{'Bernoulli'}.
#' @param t_interval A 2-element vector, represents an interval,
#' means the domain of the functional covariate.
#' Default is \code{c(0,1)}, represent interval \eqn{[0,1]}.
#' @param t_points the measurement points of functional variables,
#' should be numeric vector.
#' @param n_t Number of measurement time points.
#' Overwritten if argument \code{t_points} is assigned.
#' @param seed Pseudo-random number generation seed.
#'
#' @returns return a list with following elements.
#'    \item{Y}{An atomic vector of response variable}
#'    \item{Z}{A dataframe with a binary and a continuous scalar-valued covariate.}
#'    \item{FC}{A list of two 'functional_variable' class object.}
#'    \item{t_interval}{Same as in the input argument.}
#'    \item{t_points}{Sequence of the measurement (time) points.}
#' @export
#' @examples
#' dat_sim = MECfda_simDataGen_fcReg(100,"Bernoulli")
#' res = fcRegression(FC = dat_sim$FC, Y=dat_sim$Y, Z=dat_sim$Z,
#'                    basis.order = 3, basis.type = c('Fourier'),
#'                    family = binomial(link = "logit"))
#' @importFrom MASS mvrnorm
#' @importFrom methods hasArg
MECfda_simDataGen_fcReg = function(N = 100, distribution = c('Gaussian','Bernoulli'),
                                   t_interval,t_points,n_t=100,seed = 0){
  sigmoid = function(x) 1/(1+exp(-x))
  Stationary_Gaussian_Process_Generation = function(time_points, ker_fun, mu = 0, Sigma_matrix = NULL){
    n = length(time_points)
    if(is.function(mu)){ mu = mu(time_points)
    } else if (length(mu) == 1){ mu = rep(mu,n) }
    if(is.null(Sigma_matrix)){
      ones_vector = as.matrix(rep(1,length(time_points)))
      time_points = as.matrix(time_points)
      m1 = time_points %*% t(ones_vector)
      diff_matrix = m1 - t(m1)
      Sigma_matrix = ker_fun(diff_matrix)
    }
    MASS::mvrnorm(1,mu,Sigma_matrix)
  }
  if(methods::hasArg(t_points)){
    n_t = length(t_points)
    if(!methods::hasArg(t_interval)){
      t_interval = range(t_points)
    }else{
      if(any(t_points<t_interval[1]) | any(t_points>t_interval[2])){
        stop("Value range of t_points doens't match t_interval.")
      }
    }
    t_0 = t_interval[1]
    period = t_interval[2] - t_interval[1]
  }else{
    if(!methods::hasArg(t_interval)){
      t_interval = c(0,1)
    }
    t_0 = t_interval[1]
    period = t_interval[2] - t_interval[1]
    # t_points = t_0 + period*(1:n_t - 0.5)/n_t
    t_points = seq(t_interval[1],t_interval[2],length.out = n_t)
  }

  { ## set simulation parameters
    # N = N
    # n_t = n_t
    ker = function(r, p = c(1,1)) p[1] * exp(- (r/p[2])^2 /2)
    ElogX = 1
    p_ker_X = c(0.8,0.2)
    ker_fun_X = function(r, p = p_ker_X) ker(r, p = p)
    Betafunc = function(t, met){
      switch(met,
             f1 = sin(2*pi*(t-t_0)/period),
             f2 = cos(2*pi*(t-t_0)/period),
      )
    }
    gama = c(gamma_0 = 5,
             gamma_1 = 0.2,
             gamma_2 = 0.4)

    proportion_of_zeros = 0.2
    sigma_0 = 0.02
    switch (distribution,
            'Gaussian'  = {gama['gamma_0'] = 5},
            'Bernoulli' = {gama['gamma_0'] = -0.3}
    )
  }



  set.seed(seed)
  X1 = exp(t(sapply(rep(1,N), function(xx)Stationary_Gaussian_Process_Generation(t_points,ker_fun_X,ElogX))))
  X2 = exp(t(sapply(rep(1,N), function(xx)Stationary_Gaussian_Process_Generation(t_points,ker_fun_X,ElogX))))
  # df.FC = list(as.data.frame(X1), as.data.frame(X2))
  # mt.FC = list(as.matrix(X1),     as.matrix(X2))
  # vc.FC = list(X1[1,],X2[1,])
  Z = cbind(Z_1 = rnorm(N, 0, 1),Z_2 = rbinom(N, 1,0.6))
  switch (distribution,
          'Gaussian'  = {
            Y = gama['gamma_0'] +
              c(    period*crossprod(t(X1),Betafunc(t_points,"f1"))/n_t +
                      period*crossprod(t(X2),Betafunc(t_points,"f2"))/n_t +
                      gama['gamma_1']*Z[,1] + gama['gamma_2']*Z[,2]+
                      rnorm(N, 0,sigma_0))

          },
          'Bernoulli' = {
            Y = sigmoid(gama['gamma_0'] +
                          c(period*crossprod(t(X1),Betafunc(t_points,"f1"))/n_t +
                              period*crossprod(t(X2),Betafunc(t_points,"f2"))/n_t +
                              gama['gamma_1']*Z[,1] + gama['gamma_2']*Z[,2]))>runif(N)
          }
  )

  X1 = functional_variable(X=X1,t_0=t_0,period=period,t_points = t_points)
  X2 = functional_variable(X=X2,t_0=t_0,period=period,t_points = t_points)
  ret = list(Y = Y, Z = Z, FC = list(X1, X2),t_interval = c(t_0,period),t_points = t_points)
  rm(X1, X2)
  rm(n_t,t_0,period,t_points,ker,ElogX,p_ker_X,ker_fun_X,Betafunc,gama,N,proportion_of_zeros,sigma_0)
  return(ret)

}


# dat_sim1 = MECfda_simDataGen_fcReg(100,"Gaussian",t_interval=c(0,2),n_t=101)
# res1 = fcRegression(FC = dat_sim1$FC, Y=dat_sim1$Y, Z=dat_sim1$Z,
#                     basis.order = 6, basis.type = c('Bspline'),
#                     family = gaussian(link = "identity"))
# plot(res1$FC.BasisCoefficient[[1]])
# plot(res1$FC.BasisCoefficient[[2]])
#
# dat_sim2 = MECfda_simDataGen_fcReg(1000,"Bernoulli",t_interval=c(0,2),n_t=101)
# res2 = fcRegression(FC = dat_sim2$FC, Y=dat_sim2$Y, Z=dat_sim2$Z,
#                     basis.order = 3, basis.type = c('Fourier'),
#                     family = binomial(link = "logit"))
# plot(res2$FC.BasisCoefficient[[1]])
# plot(res2$FC.BasisCoefficient[[2]])
