inte_part=function(operator,epsilon,sigma2){
  re=(operator^(-1)*epsilon*sin(operator*epsilon)-sigma2*cos(operator*epsilon))*exp(operator^2*sigma2/2)
  return(as.numeric(re))
}
