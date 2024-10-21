first_derivative=function(e, sig, h, degree, func){
  gauss=gauss.quad(degree,c(0,1/h))
  nodes=gauss$pt
  weight=gauss$wt
  re=pi^(-1)*t(sapply(nodes,func,e,sig))%*%weight
  return(re)
}
