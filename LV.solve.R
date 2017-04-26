# Function to numerically solve de LV model (LV_model.R) using the 'ode' function from the 'desolve' package
#
# required input:
#                 - mu        : vector containing the per-capita growth rates in the absence of stress for all species 
#                 - K         : vector containing the carrying capacities in the absence of stress for all species 
#                 - a         : matrix containing the strength of intraspecific (diagonal elements) and intraspecific interaction (ofdiagonal elements)
#                 - time      : duration for which the model needs to be solved
#                 - N0        : initial species densities
#                 - conc      : matrix containing metal concentrations over time, each colum corresponds to a metal, note that the same order is required as for the parameters of the mixture toxicity model 
#                               a vector can be used when concentrations remain constant over time
#                 - model     : mixture toxicity model "CA"=concentration addition, "SA" = synergism/antagonism, "DR" = dose-ratio 
#                 - drc       : dose response curve "LL" = log-logistic(only log-logistic is currently implemented)
#                 - moa       : mode of action "pcg"= per-capita growth, "K"= carrying capacity
#                 - pars      : list of parameters for dose response curve (see 'inverse.drc.R' for details), each species has its own sublist
#                 - pars.mix  : list of parameters for mixture models (see 'mixture.models.R' for details)

LV.solve      <- function(...){
  
                      # load functions and packages
  
                              # load Lotka-Volterra model
                              source("LV.model.R")
          
                              # load mixture toxicity models
                              source("mixture.models.R")
                              
                              # load inverse dose-response curves
                              source("inverse.drc.R")
                          
                              # load package to solve nonlinear least squares
                              require("nleqslv")
                              
                              # load package to solve differential eqautions
                              require("deSolve")
            
                        # read input and detect missing parameters
                        input       <- list(...)
                        if(is.numeric(input$mu[1])==T)    {mu        <- input$mu}      else {warning("'mu' must be a vector!")}
                        if(is.numeric(input$a[1,1])==T)   {interact  <- a}             else {warning("'a' must be a matrix!")}
                        if(is.numeric(input$time[1])==T)  {time      <- input$time}    else {warning("'time must be a matrix!")}
                        if(is.numeric(input$conc[1])==T)  {conc      <- conc}          else {warning('conc must be a vector or matrix!')}
                        if(input$model=="CA")             {mixt.eff  <- conc.ad}       else if (input$model=="SA") {mixt.eff <- syn.ant}   else if (input$model=="DR") {mixt.eff <- dose.ratio} else {warning("invalid mixture model!")}
                        if(input$drc=="LL")               {drc       <- LL.drc}        else {warning("invalid dose response curve!")}
                        if(is.null(input$par)==F)         {par       <- input$par}     else {warning("parameter values for dose response curve are missing!")}     
                        if(is.null(input$par.mix)==F)     {par.mix   <- input$par.mix} else {if(input$model!="CA"){waring("parameters mixture model are missing!")}}
                        if(is.null(input$moa)==F)         {moa       <- input$moa}     else {warning("please speficify mode of action")}    

                        # solve LV model
                        
                              # step 1: convert concentrations to matrix if this would not be the case
                              
                              if(is.vector(conc)==T){
                                              
                                              if(is.vector(conc)==T){conc <- t(as.matrix(conc))}
                                              conc    <- t(sapply(0:time,function(i) conc))}
                        
                              # step 2: calculate mixture effects at each timestep and each species
                              n.spec   <- length(N0)
                              effects  <- sapply(1:n.spec, function(j) sapply(1:(time+1),function(i) nleqslv(.1,mixt.eff,C=conc[i,],drc=drc,par.drc=par[[j]])$x))
                              
                              mu.0     <- mu
                              K.0      <- K
                              if(moa=="mu"){  mu   <- t(sapply(1:(time+1), function(i) effects[i,]*mu.0))
                                              K    <- t(sapply(0:time,function(i) K.0))
                              }else{          K    <- t(sapply(1:(time+1), function(i) effects[i,]*K.0))
                                              mu    <- t(sapply(0:time,function(i) mu.0))}
                        
                              # step 3 solve LV model in the absence of stress
                              LV.par        <- list(mu=mu.0,K=K.0,interact=a)
                              out.control   <- ode(N0, c(0:time), LV, LV.par)
                              colnames(out.control) <- c("time",sapply(1:n.spec,function(i) paste("N",i,sep="")))
                        
                              # step 4 solve LV model for stress exposure
                              out.stress     <- array(dim=c((time+1),(n.spec+1)))
                              colnames(out.stress) <- c("time",sapply(1:n.spec,function(i) paste("N",i,sep="")))
                              out.stress[,1] <- 0:time
                              for(i in 1:(time+1)){
                                out.stress[i,2:(n.spec+1)] <- N0
                                LV.par        <- list(mu=mu[i,],K=K[i,],interact=a)
                                N0            <- ode(N0, 0:1, LV, LV.par)[2,2:(n.spec+1)]}
                                
                        # return output
                        output <- list(control=out.control,stress=out.stress)
                        return(output)}

