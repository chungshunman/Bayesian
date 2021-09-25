
# ts_prior: Y, tau, M, p
# vbar = VB_OLS, aols = B_OLS
# a0 = A_OLS, ssig1=sigma, a02mo = VA_OLS

ts_prior<-function(Y, tau, M, p){
  yt = t(Y[(p+1):(tau+p),]); colnames(yt)=NULL
  
  m = M + p*(M**2)
  Zt=vector()
  for(i in (p+1):(tau+p)){
    ztemp = diag(M)
    for(j in 1:p){
      xlag = as.double(Y[(i-j), (1:M)])
      xtemp = matrix(0, M, M*M)
      for(jj in 1:M){
        xtemp[jj, (((jj-1)*M)+1):(jj*M)]<- xlag
      }
      ztemp = cbind(ztemp, xtemp)
    }
    Zt = rbind(Zt, ztemp)
  }
  
  
  VB_OLS = matrix(0, m, m)
  xhy = matrix(0, m, 1)
  for(i in 1:tau){
    zhat1 = Zt[(((i-1)*M)+1):(i*M), ]
    VB_OLS = VB_OLS + t(zhat1) %*% zhat1
    xhy = xhy + t(zhat1)%*%yt[, i]
  }
  ## something strange with VB_OLS
  
  VB_OLS = solve(VB_OLS)
  B_OLS = VB_OLS%*%xhy # aols is B_OLS
  
  
  sse2 = diag(0, M, M)
  for(i in 1:tau){
    zhat1 = Zt[((i-1)*M+1):(i*M), ]
    sse2 = sse2 +
      (yt[,i] - as.double(zhat1%*%B_OLS)) %*% t(yt[, i] - as.double(zhat1%*%B_OLS))
  }
  
  hbar = sse2/tau
  
  VB_OLS = matrix(0, m, m)
  
  for(i in 1:tau){
    zhat1 = Zt[((i-1)*M+1):(i*M),]
    VB_OLS = VB_OLS + t(zhat1) %*% solve(hbar) %*% zhat1
  }
  
  VB_OLS = solve(VB_OLS)
  achol = chol(hbar)
  ssig = matrix(0, M, M)
  
  for(i in 1:p){
    ssig[i,i] = achol[i,i]
    for(j in 1:p){
      achol[j,i] = achol[j,i]/ssig[i,i]
    }
  }
  
  
  achol = solve(achol)
  numa = M*(M-1)/2
  A_OLS = matrix(0, numa, 1)
  ic = 1
  
  for(i in 2:M){
    for(j in 1:(i-1)){
      A_OLS[ic, 1] = achol[i, j]
      ic = ic + 1
    }
  }
  
  ssig1 = matrix(0, M, 1)
  
  for(i in 1:p){
    ssig1[i, 1] = log(ssig[i,i]**2)
  }
  
  hbar1 = solve(tau*hbar)
  hdraw = matrix(0, M, M)
  VA_OLS = matrix(0, numa, numa)
  A_OLS = matrix(0, numa, 1)
  
  for(irep in 1:4000){
    hdraw = rWishart(1, tau, hbar1)[,,1]
    hdraw = solve(hdraw)
    achol = t(chol(hdraw))
    ssig = matrix(0, M, M)
    for(i in 1:M){
      ssig[i, i] = achol[i, i]
      for(j in 1:M){
        achol[j, i] = achol[j,i]/ssig[i,i]
      }
    }
    achol = solve(achol)
    a0draw = matrix(0, numa, 1)
    ic = 1
    for(i in 2:M){
      for(j in 1:(i-1)){
        a0draw[ic, 1] = achol[i, j]
        ic = ic+1
      }
    }
    VA_OLS = VA_OLS + a0draw %*% t(a0draw)
    A_OLS = A_OLS + a0draw
  }
  VA_OLS = VA_OLS/4000
  A_OLS = A_OLS/4000
  VA_OLS = VA_OLS - A_OLS %*% t(A_OLS)
  
  return(list(B_OLS, VB_OLS, A_OLS, VA_OLS, ssig1))
}  
