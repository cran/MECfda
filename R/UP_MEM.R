UP_MEM = function( model="gaussian",smooth=FALSE, data,silent = FALSE){
  # library(lme4)
  # library(dplyr)
  m_w = dim(data)[3]
  t = dim(data)[2]
  n = dim(data)[1]
  M = NULL
  i=1
  while(i <= m_w){
    M = rbind(M,data[,,i])
    i=i+1
  }
  data.M = data.frame(M = I(M),seqn= rep(1:n, m_w) )
  if(silent == FALSE) print("Step 1: Massive Univariate Mixed Models")
  if(model =="gaussian"){
    fit = data.M %>%
      dplyr::select(M) %>%
      apply( 2, function(s) {
        temp= data.frame(M=s, seqn=data.M$seqn)
        mod =  suppressMessages(lmer(M ~ (1 | seqn), data = temp))
        return(predict(mod))
      })
  }else{
    fit = data.M %>%
      dplyr::select(M) %>%
      apply( 2, function(s) {
        temp= data.frame(M=s, seqn=data.M$seqn)
        mod =  suppressMessages(glmer(M ~ (1 | seqn), data = temp, family = model))
        return(predict(mod))
      })

  }
  approx= fit[!duplicated(data.M$seqn),]
  if(smooth==TRUE){
    # library(mgcv)
    if(silent == FALSE) print("Step 2: Smoothing")
    nknots <- min(round(t/4), 35)
    argvals = seq(0,1, length.out = t)
    approx= t(apply(approx, 1, function(x) gam(x ~ s(argvals, bs = "cr", k = (nknots + 1)), method = "REML")$fitted.values))
  }
  return(approx)
}
