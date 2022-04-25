library(nimble)
library(coda)
source("simSCR.Multi.R")
source("simSCR.R")
source("init.SCR.Multi.R")
source("NimbleModel Multisession SCR Poisson RJ.R")
source("NimbleFunctions Multisession SCR Poisson.R")
source("sSampler Multi.R")

#make sure to run this line!
nimble:::setNimbleOption('MCMCjointlySamplePredictiveBranches', FALSE)
nimbleOptions('MCMCjointlySamplePredictiveBranches') 

#simulate some data
N.session=3
lambda=50 #expected N
lam0=rep(0.25,N.session)
sigma=rep(0.5,N.session)
K=c(5,6,7)
buff=rep(3,N.session) #state space buffer. Should be at least 3 sigma.
X=vector("list",N.session) #one trapping array per session
X[[1]]=as.matrix(expand.grid(3:11,3:11))
X[[2]]=as.matrix(expand.grid(3:10,3:10)+1)
X[[3]]=as.matrix(expand.grid(3:9,3:9)+2)

#get state space areas
area=get.area(X,buff)
area

obstype="poisson"

#Simulate some data
D=0.4
lambda=D*area #expected N
lambda

N=rpois(N.session,lambda)#realized N
N

data=simSCR.Multi(N.session=N.session,N=N,lam0=lam0,sigma=sigma,K=K,X=X,buff=buff,obstype=obstype)

#Data augmentation level for each session
M=c(125,115,105)
#initialize multisession data structures
nimbuild=init.SCR.Multi(data,M)

z.data=matrix(NA,N.session,max(M))
for(g in 1:N.session){
  z.data[g,1:nimbuild$n[g]]=1
}

#inits for nimble
N.init=rowSums(nimbuild$z.init)
Niminits <- list(z=nimbuild$z.init,s=nimbuild$s.init,lam0=1,sigma=1,N=N.init,D=0.5)

#constants for Nimble
constants<-list(N.session=N.session,M=M,J=nimbuild$J,K2D=nimbuild$K2D,xlim=nimbuild$xlim,ylim=nimbuild$ylim,area=area)

#supply data to nimble
Nimdata<-list(y=nimbuild$y,z=z.data,X=nimbuild$X)

# set parameters to monitor
parameters<-c('lambda','lam0','sigma','N','D')

nt=1 #thinning rate

# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,useConjugacy = TRUE) 

###*required* sampler replacement
z.ups=c(25,25,25) # how many z proposals per iteration per session?
J=nimbuild$J
conf$removeSampler("N")
for(g in 1:N.session){
  #nodes used for update, calcNodes + z nodes
  y.nodes <- Rmodel$expandNodeNames(paste("y[",g,",","1:",M[g],",1:",J[g],"]"))
  lam.nodes <- Rmodel$expandNodeNames(paste("lam[",g,",","1:",M[g],",1:",J[g],"]"))
  N.node <- Rmodel$expandNodeNames(paste("N[",g,"]"))
  z.nodes <- Rmodel$expandNodeNames(paste("z[",g,",","1:",M[g],"]"))
  calcNodes <- c(N.node,y.nodes,lam.nodes)
  
  conf$addSampler(target = c("N"),
                  type = 'zSampler',control = list(inds.detected=1:nimbuild$n[g],z.ups=z.ups[g],J=J[g],M=M[g],
                                                   xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],g=g,
                                                   y.nodes=y.nodes,lam.nodes=lam.nodes,N.node=N.node,z.nodes=z.nodes,
                                                   calcNodes=calcNodes),
                  silent = TRUE)
}

#replace default activity center sampler that updates x and y locations separately with a joint update
#a little more efficient. sSampler below only tunes s when z=1. Should better tune activity centers for 
#uncaptured individuals
for(g in 1:N.session){
  conf$removeSampler(paste("s[",g,",1:",M[g],", 1:2]", sep=""))
  for(i in 1:M[g]){
    # block RW option
    # conf$addSampler(target = paste("s[",g,",",i,", 1:2]", sep=""),
    #                 type = 'RW_block',control=list(adaptive=TRUE),silent = TRUE)
    conf$addSampler(target = paste("s[",g,",",i,", 1:2]", sep=""),
                    type = 'sSampler',control=list(i=i,g=g,xlim=nimbuild$xlim[g,],ylim=nimbuild$ylim[g,],scale=1),silent = TRUE)
    #scale parameter here is just the starting scale. It will be tuned.
  }
}
# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #this will run in R, used for debugging
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


# Run the model.
start.time2<-Sys.time()
Cmcmc$run(2500,reset=FALSE) #short run for demonstration. can keep running this line to get more samples
end.time<-Sys.time()
end.time-start.time  # total time for compilation, replacing samplers, and fitting
end.time-start.time2 # post-compilation run time

library(coda)
mvSamples = as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[2:nrow(mvSamples),]))
