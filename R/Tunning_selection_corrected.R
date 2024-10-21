Tunning_selection_corrected=function(n,tau,grid_h,grid_k,estimated_coef,Y,estimated_Score,Sigma_hat,degree=40,p_z=p_z){ 
  
  num_h=length(grid_h)
  num_k=length(grid_k)
  p=dim(estimated_Score)[[2]]+1
  W_star=list(NULL)
  length(W_star)=20
  Z_plus=matrix(0,nrow = n,ncol = p_z)
  for(j in 1:20){
    W_star[[j]]=estimated_Score+cbind(Z_plus,mvrnorm(n,rep(0,(p-1-p_z)),Sigma_hat))
  }
  
  all_parameter=cbind(rep(grid_k,each=num_h),rep(grid_h,times=num_k))
  
  
  beta_hat_star=as.data.frame(matrix(numeric(0),nrow=p))
  
  for(j in 1:20){
    beta_hat_star=cbind(beta_hat_star,apply(all_parameter, MARGIN=1, FUN = Multiple_grid_corrected, Y=Y,estimated_Score=W_star[[j]],Sigma_hat=Sigma_hat,tau=tau,degree=degree,p_z=p_z))
  }
  mean_star=t(apply(estimated_coef, MARGIN = 1, FUN=rep, times=20))
  difference_star=cbind(rep(1:20,each=num_h*num_k),cbind(rep(1:(num_h*num_k),times=20),t(beta_hat_star-mean_star)))  
  difference_star=difference_star[order(difference_star[,2]),]    
  
  M_1=sapply(1:(num_h*num_k), M_fun, da=difference_star[,-(1:2)], num_h=num_h, grid_k=grid_k,p_z=p_z)
  
  W_star2=list(NULL)
  length(W_star2)=20
  for(j in 1:20){
    W_star2[[j]]=W_star[[j]]+cbind(Z_plus,mvrnorm(n,rep(0,(p-1-p_z)),Sigma_hat))
  }
  
  beta_hat_star2=as.data.frame(matrix(numeric(0),nrow=p))
  for(j in 1:20){
    beta_hat_star2=cbind(beta_hat_star2,apply(all_parameter, MARGIN=1, FUN = Multiple_grid_corrected, Y=Y,estimated_Score=W_star2[[j]],Sigma_hat=Sigma_hat,tau=tau,degree=degree,p_z=p_z))
  }
  difference_star2=cbind(rep(1:20,each=num_h*num_k),cbind(rep(1:(num_h*num_k),times=20),t(beta_hat_star2-beta_hat_star)))
  difference_star2=difference_star2[order(difference_star2[,2]),] 
  
  M_2=sapply(1:(num_h*num_k), M_fun, da=difference_star2[,-(1:2)], num_h=num_h, grid_k=grid_k,p_z=p_z)
  
  return(list(M_1=M_1,M_2=M_2,beta_hat_star=beta_hat_star,beta_hat_star2=beta_hat_star2))
}
