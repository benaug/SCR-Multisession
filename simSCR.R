e2dist = function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}
simSCR<-
  function(N=120,lam0=0.1,sigma=0.50,theta=NA,K=10,X=X,buff=3,obstype="poisson"){
    # simulate a population of activity centers
    xlim=range(X[,1])+c(-buff,buff)
    ylim=range(X[,2])+c(-buff,buff)
    s<- cbind(runif(N, xlim[1],xlim[2]), runif(N,ylim[1],ylim[2]))
    D<- e2dist(s,X)
    lamd<- lam0*exp(-D*D/(2*sigma*sigma))
    J<- nrow(X)
    
    # Capture individuals
    y=array(0,dim=c(N,J,K))
    if(obstype=="bernoulli"){
      pd=1-exp(-lamd)
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y[i,j,k]=rbinom(1,1,pd[i,j])
          }
        }
      }
    }else if(obstype=="poisson"){
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y[i,j,k]=rpois(1,lamd[i,j])
          }
        }
      } 
    }else if(obstype=="negbin"){
      for(i in 1:N){
        for(j in 1:J){
          for(k in 1:K){
            y[i,j,k]=rnbinom(1,mu=lamd[i,j],size=theta)
          }
        }
      } 
    }else{
      stop("obstype not recognized")
    }

    #discard uncaptured inds and aggregate true IDcovs for all samples, keeping track of where they came from with A matrix (used to put observed data back together)
    caught=which(apply(y,c(1),sum)>0)
    y.true=y
    y=y[caught,,]
    if(K==1){
      y=array(y,dim=c(dim(y),1))
    }
    n=length(caught)
    
    out<-list(y=y,X=X,K=K,buff=buff,obstype=obstype,s=s,n=nrow(y),K=K,
              xlim=xlim,ylim=ylim)
    return(out)
  }
