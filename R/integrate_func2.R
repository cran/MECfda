integrate_func2=function(operator, e, sig){
  re=(2*cos(operator*e)-operator*e*sin(operator*e)+sig*operator^2*cos(operator*e))*exp(operator^2*sig/2)
  return(as.numeric(re))
}
