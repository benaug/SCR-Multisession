#------------------------------------------------------------------
# Function for calculation detection rate
#------------------------------------------------------------------
GetDetectionRate <- nimbleFunction(
  run = function(s = double(1), lam0=double(0), sigma=double(0), 
                 X=double(2), J=double(0), z=double(0)){ 
    returnType(double(1))
    if(z==0) return(rep(0,J))
    if(z==1){
      d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
      ans <- lam0*exp(-d2/(2*sigma^2))
      return(ans)
    }
  }
)

dPoissonVector <- nimbleFunction(
  run = function(x = double(1), lam = double(1), z = double(0),
                 log = integer(0)) {
    returnType(double(0))
    if(z==0){
      return(0)
    }else{
      logProb <- sum(dpois(x, lambda = lam, log = TRUE))
      return(logProb)
    }
  }
)

#make dummy random vector generator to make nimble happy
rPoissonVector <- nimbleFunction(
  run = function(n = integer(0),lam = double(1),z = double(0)) {
    returnType(double(1))
    J <- nimDim(lam)[1]
    out <- numeric(J,value=0)
    return(out)
  }
)

#Required custom update for N/z
zSampler <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    g <- control$g
    J <- control$J
    M <- control$M
    inds.detected <- control$inds.detected
    z.ups <- control$z.ups
    xlim <- control$xlim
    ylim <- control$ylim
    y.nodes <- control$y.nodes
    lam.nodes <- control$lam.nodes
    N.node <- control$N.node
    z.nodes <- control$z.nodes
    calcNodes <- control$calcNodes
  },
  run = function() {
    for(up in 1:z.ups){ #how many updates per iteration?
      #propose to add/subtract 1
      updown <- rbinom(1,1,0.5) #p=0.5 is symmetric. If you change this, must account for asymmetric proposal
      reject <- FALSE #we auto reject if you select a detected individual
      if(updown==0){#subtract
        #find all z's currently on
        z.on <- which(model$z[g,1:M]==1)
        n.z.on <- length(z.on)
        pick <- rcat(1,rep(1/n.z.on,n.z.on)) #select one of these individuals
        pick <- z.on[pick]
        if(any(pick==inds.detected)){ #is this individual detected?
          reject <- TRUE #if so, we reject (could never select these inds, but then need to account for asymmetric proposal)
        }
        if(!reject){

          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick])

          #propose new N/z
          model$N[g] <<-  model$N[g] - 1
          model$z[g,pick] <<- 0

          #turn lam off
          model$calculate(lam.nodes[pick])

          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick]) #will always be 0

          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)

          if(accept) {
            mvSaved["N",1][g] <<- model[["N"]][g]
            for(j in 1:J){
              mvSaved["lam",1][g,pick,j] <<- model[["lam"]][g,pick,j]
            }
            mvSaved["z",1][g,pick] <<- model[["z"]][g,pick]
          }else{
            model[["N"]][g] <<- mvSaved["N",1][g]
            for(j in 1:J){
              model[["lam"]][g,pick,j] <<- mvSaved["lam",1][g,pick,j]
            }
            model[["z"]][g,pick] <<- mvSaved["z",1][g,pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }else{#add
        if(model$N[g] < M){ #cannot update if z maxed out. Need to raise M
          z.off <- which(model$z[g,1:M]==0)
          n.z.off <- length(z.off)
          pick <- rcat(1,rep(1/n.z.off,n.z.off)) #select one of these individuals
          pick <- z.off[pick]

          #get initial logprobs for N and y
          lp.initial.N <- model$getLogProb(N.node)
          lp.initial.y <- model$getLogProb(y.nodes[pick]) #will always be 0

          #propose new N/z
          model$N[g] <<-  model$N[g] + 1
          model$z[g,pick] <<- 1

          #turn lam on
          model$calculate(lam.nodes[pick])

          #get proposed logprobs for N and y
          lp.proposed.N <- model$calculate(N.node)
          lp.proposed.y <- model$calculate(y.nodes[pick])

          #MH step
          log_MH_ratio <- (lp.proposed.N + lp.proposed.y) - (lp.initial.N + lp.initial.y)
          accept <- decide(log_MH_ratio)
          if(accept) {
            mvSaved["N",1][g] <<- model[["N"]][g]
            for(j in 1:J){
              mvSaved["lam",1][g,pick,j] <<- model[["lam"]][g,pick,j]
            }
            mvSaved["z",1][g,pick] <<- model[["z"]][g,pick]
          }else{
            model[["N"]][g] <<- mvSaved["N",1][g]
            for(j in 1:J){
              model[["lam"]][g,pick,j] <<- mvSaved["lam",1][g,pick,j]
            }
            model[["z"]][g,pick] <<- mvSaved["z",1][g,pick]
            model$calculate(y.nodes[pick])
            model$calculate(N.node)
          }
        }
      }
    }
    #copy back to mySaved to update logProbs which was not done above
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)