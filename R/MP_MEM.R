MP_MEM= function(model="gaussian", d = 3,smooth=FALSE, silent = FALSE, data){
  cov.model ="us"
  m_w = dim(data)[3]
  t = dim(data)[2]
  n = dim(data)[1]
  data.M = matrix(data, prod(dim(data)[1:3])) %>% as.data.frame()
  data.M$seqn = rep(1:n, prod(t,m_w))
  data.M$time = rep(sapply(1:t, function(s) rep(s, n)),m_w)
  data.M$day = lapply(1:m_w, function(s) rep(s, n*t)) %>% unlist()
  colnames(data.M)[1] ="M"
  if(silent == FALSE) print("Step 1: Massive Short multivariate Mixed Models")
  ind.time = 1:t
  if(d%%2 == 0){
    fit = lapply(ind.time, function(s){
      if ((s-(d-2))<1){
        ind.interval = 1:d
      }else if((s+(d-2))>t){
        ind.interval = (t-(d-1)):t
      }else {
        ind.interval = (s-d/2 +1):(s+d/2)
      }
      temp = data.M %>%
        dplyr::filter(time %in% ind.interval) %>%
        dplyr::mutate(time.cat = as.character(time))
      if(model=="gaussian"){
        if(cov.model =="us"){
          mod =  suppressMessages(lmer(M ~ (1| seqn)+(0 + time.cat | seqn), data = temp))
        }else if (cov.model =="id"){
          mod = suppressMessages(lmer(M ~ (1 | seqn)+(1|time.cat), data = temp))
        }
        else if (cov.model =="cs"){
          mod = suppressMessages(lme(M ~1, random = list(seqn= ~1, seqn =pdCompSymm(~ time.cat-1)),
                                     data = temp))
        }else {
          stop ("wrong covariance matrix of random effect")
        }
      }else{
        if(cov.model =="us"){
          mod =  suppressMessages(glmer(M ~ (1| seqn)+(0 + time.cat | seqn), data = temp, family = model))
        }else if (cov.model =="id"){
          mod = suppressMessages(glmer(M ~ (1 | seqn)+(1|time.cat), data = temp, family = model))
        }

        else if (cov.model =="cs"){
          mod = suppressMessages(glme(M ~1, random = list(seqn= ~1, seqn =pdCompSymm(~ time.cat-1)),
                                      data = temp))
        }else {
          stop ("wrong covariance matrix of random effect")
        }
      }
      pred = predict(mod, newdata = data.frame(time.cat = as.character(s), seqn = unique(temp$seqn)))
      return(pred)
    })
  }else{
    fit = lapply(ind.time, function(s){
      if ((s-(d-2))<1){
        ind.interval = 1:d
      }else if((s+(d-2))>t){
        ind.interval = (t-(d-1)):t
      }else {
        ind.interval = (s-floor(d/2)):(s+floor(d/2))
      }
      temp = data.M %>%
        dplyr::filter(time %in% ind.interval) %>%
        dplyr::mutate(time.cat = as.character(time))
      if(model=="gaussian"){
        if(cov.model =="us"){
          mod =  suppressMessages(lmer(M ~ (1| seqn)+(0 + time.cat | seqn), data = temp))
        }else if (cov.model =="id"){
          mod = suppressMessages(lmer(M ~ (1 | seqn)+(1|time.cat), data = temp))
        }
        else if (cov.model =="cs"){
          mod = suppressMessages(lme(M ~1, random = list(seqn= ~1, seqn =pdCompSymm(~ time.cat-1)),
                                     data = temp))
        }else {
          stop ("wrong covariance matrix of random effect")
        }
      }else{
        if(cov.model =="us"){
          mod =  suppressMessages(glmer(M ~ (1| seqn)+(0 + time.cat | seqn), data = temp, family = model))
        }else if (cov.model =="id"){
          mod = suppressMessages(glmer(M ~ (1 | seqn)+(1|time.cat), data = temp, family = model))
        }
        else if (cov.model =="cs"){
          mod = suppressMessages(glme(M ~1, random = list(seqn= ~1, seqn =pdCompSymm(~ time.cat-1)),
                                      data = temp))
        }else {
          stop ("wrong covariance matrix of random effect")
        }
      }
      pred = predict(mod, newdata = data.frame(time.cat = as.character(s), seqn = unique(temp$seqn)))
      return(pred)
    })
  }
  approx= matrix(unlist(fit), nrow= n, ncol = t)
  if(smooth==TRUE){
    # library(mgcv)
    if(silent == FALSE) print("Step 2: Smoothing")
    nknots <- min(round(t/4), 35)
    argvals = seq(0,1, length.out = t)
    approx= t(apply(approx, 1, function(x) gam(x ~ s(argvals, bs = "cr", k = (nknots + 1)), method = "REML")$fitted.values))
  }
  return(approx)
}
