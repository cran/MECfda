#' @title Bias correction method of applying quantile linear regression to dataset
#' with one functional covariate with measurement error using corrected loss score method.
#' @description
#' Zhang et al. proposed a new corrected loss function for a
#' partially functional linear quantile model with functional measurement error in this manuscript.
#' They established a corrected quantile objective function of the observed measurement
#' that is an unbiased estimator of the quantile objective function
#' that would be obtained if the true measurements were available.
#' The estimators of the regression parameters are obtained
#' by optimizing the resulting corrected loss function.
#' The resulting estimator of the regression parameters is shown to be consistent.
#' @references Zhang, Mengli, et al.
#' "PARTIALLY FUNCTIONAL LINEAR QUANTILE REGRESSION WITH MEASUREMENT ERRORS."
#' Statistica Sinica 33 (2023): 2257-2280.
#' @param data.Y Response variable, can be an atomic vector, a one-column matrix or data frame,
#' recommended form is a one-column data frame with column name.
#' @param data.W A 3-dimensional array, represents \eqn{W}, the measurement of \eqn{X}.
#' Each row represents a subject.
#' Each column represent a measurement (time) point.
#' Each layer represents an observation.
#' @param data.Z Scalar covariate(s),
#' can be not input or NULL (when there's no scalar covariate),
#' an atomic vector (when only one scalar covariate),
#' a matrix or data frame, recommended form is a data frame with column name(s).
#' @param tau Quantile \eqn{\tau\in(0,1)}, default is 0.5.
#' @param t_interval A 2-element vector, represents an interval,
#' means the domain of the functional covariate. Default is c(0,1), represent interval \eqn{[0,1]}.
#' @param t_points Sequence of the measurement (time) points, default is NULL.
#' @param grid_k An atomic vector, of which each element is candidate number of basis.
#' @param grid_h A non-zero-value atomic vector, of which each element is candidate value of tunning parameter.
#' @param degree Used in computation for derivative and integral, defult is 45, large enough for most scenario.
#' @param observed_X For estimating parametric variance. Default is NULL.
#'
#' @return Returns a ME.fcQR_CLS class object. It is a list that contains the following elements.
#'    \item{estimated_beta_hat}{Estimated coefficients from corrected loss function (including functional part)}
#'    \item{estimated_beta_t}{Estimated functional curve}
#'    \item{SE_est}{Estimated parametric variance. Returned only if observed_X is not NULL.}
#'    \item{estimated_Xbasis}{The basis matrix we used}
#'    \item{res_naive}{results of naive method}
#' @export
#' @examples
#' data(MECfda.data.sim.0.1)
#' \donttest{
#' res = ME.fcQR_CLS(data.Y = MECfda.data.sim.0.1$Y,
#'                 data.W = MECfda.data.sim.0.1$W,
#'                data.Z = MECfda.data.sim.0.1$Z,
#'                tau = 0.5,
#'                grid_k = 4:7,
#'                grid_h = 1:2)
#' }
#' @importFrom MASS mvrnorm
#' @importFrom quantreg rq
#' @importFrom fda bsplineS
#' @importFrom gss gauss.quad
#' @import stats
ME.fcQR_CLS=function(data.Y, data.W, data.Z, tau = 0.5,
                     t_interval = c(0,1), t_points = NULL,
                     grid_k, grid_h, degree=45,
                     observed_X = NULL
                     ){
  bs_degree = 3
  data.Y = as.data.frame(data.Y)
  if(is.null(colnames(data.Y))) colnames(data.Y) = 'Y'
  Y = data.Y[,]



  if(!methods::hasArg(data.Z)){
    data.Z = NULL
  }
  if(!is.null(data.Z)){
    Z = as.data.frame(data.Z)
    if(is.null(colnames(Z))) colnames(Z) = paste('Z',1:ncol(Z),sep = '_')
    if(any(
      length(Y) != nrow(data.W),
      length(Y) != nrow(Z)
    )){
      stop("dimensionality of variables doesn't match")
    }
  }else{
    Z = as.data.frame(matrix(nrow = nrow(Y),ncol = 0))
  }
  Z = as.matrix(data.Z)




  observed_W = apply(data.W, 3, as.matrix, simplify = FALSE)

  repeat_num=length(observed_W)
  n=dim(observed_W[[1]])[1]
  t_length=dim(observed_W[[1]])[2]

  if(is.null(t_points)){
    observed_t=seq(t_interval[1],t_interval[2],length.out = t_length)
  }else{
    observed_t = t_points
  }

  Bs=bsplineS(observed_t, breaks = seq(t_interval[1],t_interval[2],length.out = 5), norder = bs_degree + 1)

  p_z=dim(Z)[2]


  W=matrix(0, nrow = n, ncol = t_length)
  for(i in 1:repeat_num){
    W=W+observed_W[[i]]
  }
  W=W/repeat_num

  Kernal_U=matrix(0,t_length,t_length)
  for(i in 1:repeat_num){
    Kernal_U=Kernal_U+t(observed_W[[i]]-W)%*%(observed_W[[i]]-W)
  }
  Kernal_U=Kernal_U/(n*(repeat_num-1))

  overall_mean=colSums(W)/n
  Kernal_W=(t(W)-overall_mean)%*%t(t(W)-overall_mean)/(n-1)
  Kernal_X=Kernal_W-Kernal_U/repeat_num

  delta=1/(t_length-1)
  max_k=max(grid_k)



  naive_kernal=(t(W)-colMeans(W))%*%t(t(W)-colMeans(W))/n
  naive_basis=as.matrix(delta^(-0.5)*eigen(naive_kernal)$vectors[,1:max_k])
  naive_basis=Bs%*%lm(naive_basis~Bs-1)$coefficients
  naive_score=delta*t(t(W)-colMeans(W))%*%naive_basis

  naive_predictor=cbind(Z,naive_score)

  BIC_naive=sapply(grid_k, L_BIC, Y=Y, pred=naive_predictor, tau=tau, p_z=p_z)
  K_optim_naive=grid_k[which.min(BIC_naive)]


  naive_coeff=rq(Y~naive_predictor[,1:(K_optim_naive+p_z)],tau,method = "fn")$coefficient
  naive_beta_t=as.matrix(naive_basis[,1:K_optim_naive])%*%naive_coeff[-(1:(p_z+1))]


  estimated_Xbasis=as.matrix(delta^(-0.5)*eigen(make.positive.definite(Kernal_X))$vectors[,1:K_optim_naive])
  estimated_Xbasis=Bs%*%lm(estimated_Xbasis~Bs-1)$coefficients

  estimated_Score=list()
  Score_mean=matrix(0, nrow = n, ncol = K_optim_naive)
  for(i in 1:repeat_num){
    estimated_Score[[i]]=delta*t(t(observed_W[[i]])-colMeans(observed_W[[i]]))%*%estimated_Xbasis
    Score_mean=Score_mean+estimated_Score[[i]]
  }
  Score_mean=Score_mean/repeat_num

  Sigma_hat=matrix(0, nrow = K_optim_naive, ncol = K_optim_naive)
  for(i in 1:repeat_num){
    Sigma_hat=Sigma_hat+t(estimated_Score[[i]]-Score_mean)%*%(estimated_Score[[i]]-Score_mean)
  }
  Sigma_hat=1/repeat_num*Sigma_hat/(n*(repeat_num-1))
  method_predictor=cbind(Z,Score_mean)


  num_h=length(grid_h)
  num_k=length(K_optim_naive)
  All_parameter=cbind(rep(K_optim_naive,each=num_h),rep(grid_h,times=num_k))
  All_coeff=apply(All_parameter, MARGIN=1, FUN = Multiple_grid_corrected, Y=Y,
                  estimated_Score=method_predictor, Sigma_hat=Sigma_hat, tau=tau, degree=degree, p_z=p_z)

  tunning_h=Tunning_selection_corrected(n,tau,grid_h,K_optim_naive,All_coeff,Y,
                                        method_predictor,Sigma_hat,degree=45,p_z=p_z)
  h_optim=(grid_h[which.min(tunning_h$M_1)])^2/grid_h[which.min(tunning_h$M_2)]

  initial=rq(Y~method_predictor[,1:(K_optim_naive+p_z)],tau)$coefficients
  names(initial)=NULL
  estimated_beta_hat=optim(par=initial,fn=Corrected_loss,y=Y,
                           w=as.matrix(method_predictor[,1:(K_optim_naive+p_z)]),
                           h=h_optim,Sigma_hat=as.matrix(Sigma_hat[1:K_optim_naive,1:K_optim_naive]),
                           tau=tau,degree=degree,p_z=p_z)$par
  estimated_beta_t=as.matrix(estimated_Xbasis[,1:K_optim_naive])%*%estimated_beta_hat[-(1:(p_z+1))]

  ret=list(All_coeff          = All_coeff,
           estimated_beta_hat = estimated_beta_hat,
           estimated_beta_t   = estimated_beta_t,
           tunning_h          = tunning_h,
           h_optim            = h_optim,
           estimated_Xbasis   = estimated_Xbasis)


  if(!is.null(observed_X)){
    X_score = delta*t(t(observed_X)-colMeans(observed_X))%*%estimated_Xbasis
    epsilon_est = Y-cbind(rep(1,n),Z)%*%estimated_beta_hat[1:(1+p_z)]-Score_mean%*%estimated_beta_hat[-(1:(1+p_z))]
    sigma_2_est = t(estimated_beta_hat[-(1:(1+p_z))])%*%Sigma_hat%*%estimated_beta_hat[-(1:(1+p_z))]

    Z_star_est=Z
    for(j in 1:p_z){
      Z_star_est[,j]=Z[,j]-X_score%*%solve(t(X_score)%*%X_score)%*%t(X_score)%*%Z[,j]
    }
    Z_star_est = cbind(rep(1,n), Z_star_est)

    second_diag_est = diag(sapply(epsilon_est, second_derivative, sig=sigma_2_est, h=h_optim, degree=degree, func=integrate_func2))
    D_est = t(Z_star_est)%*%second_diag_est%*%Z_star_est/n
    D_1_est = solve(D_est)


    Z_tilde_est=Z
    for(j in 1:p_z){
      Z_tilde_est[,j]=Z[,j]-Score_mean%*%solve(t(X_score)%*%X_score)%*%t(X_score)%*%Z[,j]
    }
    Z_tilde_est = cbind(rep(1,n), Z_tilde_est)
    first_diag_est = diag(tau-1/2+sapply(epsilon_est, first_derivative, sig=sigma_2_est, h=h_optim, degree=degree, func=integrate_func1))
    B_est = first_diag_est%*%Z_tilde_est
    B_est = t(B_est)%*%B_est/n

    SE_est = sqrt(diag(D_1_est%*%B_est%*%D_1_est/n))

    ret$SE_est = SE_est
  }

  res_naive=list(naive_coeff        = naive_coeff,
                 naive_beta_t       = naive_beta_t,
                 naive_basis        = naive_basis,
                 K_optim_naive      = K_optim_naive,
                 BIC_naive          = BIC_naive)
  ret$res_naive = res_naive

  class(ret) = 'ME.fcQR_CLS'
  return(ret)
}










