Corrected_loss=function(beta,y,w,h,Sigma_hat,tau,degree=30, p_z){
  n=dim(w)[1]
  p=dim(w)[2]+1
  
  sigma_hat=as.numeric(t(beta[(p_z+2):p])%*%Sigma_hat%*%beta[(p_z+2):p])
  epsilon=y-rep(beta[1],n)-w%*%beta[2:p]
  loss=sapply(epsilon,gauss_quad,sigma2=sigma_hat,tau=tau,func=inte_part,h=h,degree=degree)/n
  return(sum(loss))
}
