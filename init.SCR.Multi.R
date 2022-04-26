init.SCR.Multi<-
  function(data=data,M=M){
    N.session=length(data)
    if(length(M)!=N.session)stop("M and data must be of length 'N.session'")
    M.max=max(M)
    J=unlist(lapply(data,function(x){nrow(x$X)}))
    J.max=max(J)
    y=array(0,dim=c(N.session,M.max,J.max)) #maximal augmentation across sessions
    X=array(0,dim=c(N.session,J.max,2))
    xlim=ylim=matrix(0,N.session,2)
    K=n=rep(NA,N.session)
    K1D=matrix(0,N.session,J.max)
    for(g in 1:N.session){
      y[g,1:data[[g]]$n,1:J[g]]=apply(data[[g]]$y,c(1,2),sum)
      X[g,1:J[g],1:2]=data[[g]]$X
      xlim[g,]=data[[g]]$xlim
      ylim[g,]=data[[g]]$ylim
      J[g]=nrow(data[[g]]$X)
      K[g]=data[[g]]$K
      n[g]=data[[g]]$n
      K1D[g,1:J[g]]=rep(K[g],J[g])
    }
    
    s.init=array(NA,dim=c(N.session,M.max,2))
    for(g in 1:N.session){
      s.init[g,1:M[g],] <- cbind(runif(M[g],xlim[g,1],xlim[g,2]), runif(M[g],ylim[g,1],ylim[g,2])) #assign random locations
      idx=which(rowSums(y[g,,])>0) #switch for those actually caught
      for(i in idx){
        trps<- matrix(X[g,y[g,i,]>0,1:2],ncol=2,byrow=FALSE)
        if(nrow(trps)>1){
          s.init[g,i,]<- c(mean(trps[,1]),mean(trps[,2]))
        }else{
          s.init[g,i,]<- trps
        }
      }
    }
    z.init=1*(apply(y,c(1,2),sum)>0)
    
    return(list(y=y,X=X,xlim=xlim,ylim=ylim,K=K,J=J,n=n,K1D=K1D,s.init=s.init,z.init=z.init))
  }
