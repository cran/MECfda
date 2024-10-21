Multiple_grid_corrected=function(tunning,Y,estimated_Score,Sigma_hat,tau,degree=40,p_z){
  p=dim(estimated_Score)[[2]]+1
  re=rep(0,p)
  K=tunning[1]
  h=tunning[2]
  pred=as.matrix(estimated_Score[,1:(K+p_z)])
  initial=rq(Y~pred,tau,method="fn")$coefficients
  names(initial)=NULL
  sigma_hat=as.matrix(Sigma_hat[1:K,1:K])
  coeff=optim(par=initial,fn=Corrected_loss,y=Y,w=pred, h=h ,Sigma_hat=sigma_hat,tau=tau,degree=degree, p_z=p_z)$par
  names(coeff)=NULL
  re[1:(K+p_z+1)]=coeff
  return(re)
}
