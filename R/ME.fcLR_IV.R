#' @title Bias correction method of applying linear regression to one functional
#' covariate with measurement error using instrumental variable.
#' @description
#' See detailed model in reference
#' @references Tekwe, Carmen D., et al.
#' "Instrumental variable approach to estimating the scalar‐on‐function regression model w
#' ith measurement error with application to energy expenditure assessment in childhood obesity."
#' Statistics in medicine 38.20 (2019): 3764-3781.
#' @param data.Y Response variable, can be an atomic vector, a one-column matrix or data frame,
#' recommended form is a one-column data frame with column name.
#' @param data.W A dataframe or matrix, represents \eqn{W}, the measurement of \eqn{X}.
#' Each row represents a subject. Each column represent a measurement (time) point.
#' @param data.M A dataframe or matrix, represents \eqn{M}, the instrumental variable.
#' Each row represents a subject. Each column represent a measurement (time) point.
#' @param t_interval A 2-element vector, represents an interval,
#' means the domain of the functional covariate.
#' Default is c(0,1), represent interval \eqn{[0,1]}.
#' @param t_points Sequence of the measurement (time) points,
#' default is NULL.
#' @param CI.bootstrap Whether to return the confidence using bootstrap method.
#' Default is FALSE.
#'
#' @return Returns a ME.fcLR_IV class object. It is a list that contains the following elements.
#'    \item{beta_tW}{Parameter estimates.}
#'    \item{CI}{Confidence interval, returnd only when CI.bootstrap is TRUE. }
#' @export
#' @examples
#' data(MECfda.data.sim.0.3)
#' res = ME.fcLR_IV(data.Y = MECfda.data.sim.0.3$Y,
#'               data.W = MECfda.data.sim.0.3$W,
#'               data.M = MECfda.data.sim.0.3$M)
#' @importFrom splines bs
#' @import stats
ME.fcLR_IV = function(data.Y, data.W, data.M,
                      t_interval = c(0,1), t_points = NULL,
                      CI.bootstrap = FALSE){

  data.Y = as.data.frame(data.Y)
  if(is.null(colnames(data.Y))) colnames(data.Y) = 'Y'
  Y = as.matrix(data.Y)


  W = as.matrix(data.W)
  M = as.matrix(data.M)


  n <- length(data.Y)
  t <- ncol(data.W)

  k = ceiling(n^{1/5})+2
  k = max(k,4)
  if(is.null(t_points)){
    a  <- seq(t_interval[1], t_interval[2],length.out=t)
    b  <- seq(t_interval[1], t_interval[2],length.out=t)
  }else{
    a = t_points
    b = t_points
  }
  bs2 <- bs(b, df = k, intercept=TRUE)


  M_i <- crossprod(t(M),bs2)/length(a)
  W_i <- crossprod(t(W),bs2)/length(a)



  M_ic <- scale(M_i,scale = TRUE)
  W_ic <- scale(W_i,scale = TRUE)
  Y_c <-  scale(Y, scale = TRUE )


  omega_mw <- crossprod(M_ic,W_ic)/n
  omega_my <- crossprod(M_ic,Y_c)/n


  c_hat <- crossprod(t(solve(crossprod(omega_mw))),crossprod(omega_mw,omega_my))
  c_hat <- c_hat

  beta_t <- crossprod(t(bs2),c_hat)


  c_hatW <- lm(Y_c ~  W_ic - 1)$coefficients
  c_predW <- tcrossprod(bs2,t(c_hatW))

  beta_tW <- crossprod(t(bs2),c_hatW)


  b3  <- seq(0, 1, length.out = t)
  bs3 <- bs(b3, df=k, intercept=TRUE)
  beta_tW <- crossprod(t(bs3),c_hatW)

  ret = list(beta_tW = beta_tW)


  if(CI.bootstrap){
    B <- 5000
    iter <- 1

    c.boot <- matrix(nrow=nrow(c_hat), ncol=B)
    beta.boot <- matrix(nrow=nrow(beta_t), ncol=B)
    beta.q1 <- matrix(nrow=nrow(beta_t), 1)
    beta.q2 <- matrix(nrow=nrow(beta_t), 1)

    repeat{

      i <- iter
      id <- sample(1:nrow(M),replace = TRUE)
      M.boot <- M[id,]
      W.boot <- W[id,]
      Y.boot <- Y[id,]


      M.boot_i <- crossprod(t(M.boot),bs2)/length(a)
      W.boot_i <- crossprod(t(W.boot),bs2)/length(a)


      M.boot_ic <-  scale(M.boot_i,scale = TRUE)
      W.boot_ic <-  scale(W.boot_i,scale = TRUE)
      Y.boot_c <-   scale(Y.boot,scale = TRUE)



      omega.boot_mw <- crossprod(M.boot_ic,W.boot_ic)/n
      omega.boot_my <- crossprod(M.boot_ic,Y.boot_c)/n


      c.boot_hat <- crossprod(t(solve(crossprod(omega.boot_mw),tol=1e-20)),crossprod(omega.boot_mw,omega.boot_my))


      beta.boot_t <- crossprod(t(bs2),c.boot_hat)

      beta.boot[,i] <- beta.boot_t


      iter <- iter+1

      i


      if(iter > B){break}


    }

    beta.q1 <- matrix(apply(beta.boot,1,quantile, probs=0.025),nrow=nrow(beta_t), 1)
    beta.q2 <- matrix(apply(beta.boot,1,quantile, probs=0.975),nrow=nrow(beta_t), 1)
    mean.beta <- matrix(apply(beta.boot,1,mean),nrow=nrow(beta_t), 1)

    CI = cbind(beta.q1,beta.q2)
    colnames(CI) = c('q1','q2')
    ret$CI = CI
  }

  return(ret)
}
