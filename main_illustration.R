# main scirpt to solve the Lotka-Volterra equations for mixtures

# set working directory
setwd("/Users/Jan/Documents/Academia/UGent/Data/Karel")

# load all functions
source("LV.solve.R")

# define species parameters

        # per-capita growth rates in unstressed conditions
        mu     <- c(.1,.5)

        # carrying capacities in unstressed conditions
        K      <- c(100,500)

        # species interactions
        a      <- matrix(c(1,.5,.3,1),nrow=2)

        # initial species densities
        N0     <- c(10,10)

        # species dose response cuverse

            # species 1
          
                        # chemical 1
                        EC50.1.1    <- 50
                        beta.1.1    <- 2

                        # chemical 2
                        EC50.1.2    <- 100
                        beta.1.2    <- 5

            # species 2

                        # chemical 1
                        EC50.2.1    <- 25
                        beta.2.1    <- 1

                        # chemical 2
                        EC50.2.2    <- 150
                        beta.2.2    <- 4

# select dose response curve, mixture effects and mode of action

      drc     <- "LL"
      model   <- "CA"
      moa     <- "K"

# simulation settings

      time    <- 50
      conc    <- c(10,20)

# run simulation
pars  <- list(spec1=list(ec.50=EC50.1.1,beta=beta.1.1,ec.50=EC50.1.2,beta=beta.1.2),
          spec2=list(ec.50=EC50.2.1,beta=beta.2.1,ec.50=EC50.2.2,beta=beta.2.2))

sim   <- LV.solve(mu=mu,K=K,a=a,N0=N0,drc=drc,model=model,moa=moa,time=time,conc=conc,pars=pars) 

      #plot output
      plot(sim$control[,1:2],ylim=c(0,max(K)),type="l",xlab="Time",ylab="Density")
      points(sim$control[,c(1,3)],type="l",col="red")
      points(sim$stress[,c(1,2)],type="l",lty=2)
      points(sim$stress[,c(1,3)],type="l",lty=2,col="red")
  
