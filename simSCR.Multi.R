e2dist = function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

get.area = function (X, buff){
  N.session=length(X)
  area=rep(NA,N.session)
  for(g in 1:N.session){
    area[g]=diff(range(X[[g]][,1])+c(-buff[g],buff[g]))*diff(range(X[[g]][,2])+c(-buff[g],buff[g]))
  }
  return(area)
}

simSCR.Multi<-
  function(N.session=3,N=NA,lam0=NA,sigma=NA,theta=NA,K=NA,X=NA,buff=3,obstype="poisson"){
    data=vector("list",N.session)
    if(obstype=="poisson"){
      theta=rep(NA,N.session)
    }
    for(g in 1:N.session){
      data[[g]]=simSCR(N=N[g],lam0=lam0[g],sigma=sigma[g],theta=theta[g],
                       K=K[g],X=X[[g]],buff=buff[g],obstype=obstype)
    }
    return(data)
  }
