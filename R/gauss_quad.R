gauss_quad=function(epsilon,sigma2,tau,func,h,degree=45){
  gauss=gauss.quad(degree,c(0,1/h))
  nodes=gauss$pt
  weight=gauss$wt
  re=epsilon*(tau-0.5)+pi^(-1)*t(sapply(nodes,func,epsilon,sigma2))%*%weight
  return(re)
}
