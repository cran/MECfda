L_BIC=function(K, Y, pred, tau, p_z){
  
  coef_k=rq(Y~pred[,1:(K+p_z)],tau,method = "fn")$coefficients
  names(coef_k)=NULL
  
  n=length(Y)
  eps = Y-cbind(rep(1,n),pred[,1:(K+p_z)])%*%coef_k
  rho = eps*ifelse(eps<0, tau-1, tau)
  
  re=log(mean(rho))+(K+p_z+1)*log(n)/n
  return(re)
  
}
