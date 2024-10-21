integrate_func1=function(operator, e, sig){
  re=(operator^(-1)*sin(operator*e)+e*cos(operator*e)+sig*operator*sin(operator*e))*exp(operator^2*sig/2)
  return(as.numeric(re))
}
