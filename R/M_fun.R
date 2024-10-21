M_fun=function(ind,da,num_h,grid_k,p_z){
  ind1=(ind-1)*20+1
  ind2=ind1+19
  samp=da[ind1:ind2,]
  
  k_ind=(ind-1)%/%num_h+1
  p=grid_k[k_ind]+p_z+1
  samp=samp[,1:p]
  
  S=solve(var(samp))
  re=sum(diag(samp%*%S%*%t(samp)))/20
  return(as.numeric(re))
}
