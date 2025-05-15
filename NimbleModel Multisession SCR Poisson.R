NimModel <- nimbleCode({
  #--------------------------------------------------------------
  # priors
  #--------------------------------------------------------------
  D ~ dunif(0,100) #expected D
  lam0 ~ dunif(0,10)
  sigma ~ dunif(0,100)
  #--------------------------------------------------------------
  for(g in 1:N.session){
    lambda[g] <- D*area[g] #expected N
    N[g] ~ dpois(lambda[g]) #realized N
    for(i in 1:M[g]) {
      s[g,i,1] ~ dunif(xlim[g,1],xlim[g,2])
      s[g,i,2] ~ dunif(ylim[g,1],ylim[g,2])
      lam[g,i,1:J[g]] <- GetDetectionRate(s = s[g,i,1:2], X = X[g,1:J[g],1:2], J=J[g],sigma=sigma, lam0=lam0, z=z[g,i])
      y[g,i,1:J[g]] ~ dPoissonVector(lam=lam[g,i,1:J[g]]*K1D[g,1:J[g]],z=z[g,i]) #vectorized obs mod
    }
  }
})