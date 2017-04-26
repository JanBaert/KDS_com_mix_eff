# concentration addition model
#
#   input:  x       : effect size
#           C       : vector containing metal concentrations
#           drc     : type of dose response curve required (see 'inverse.drc.R' for available types)
#           par.drc : list of parameters for dose response curves for each chemical. Note that parameter values should be grouped per metal
#
#   output: deviation from the concentration addition model

conc.ad    <- function(x,C,drc,par.drc){   
  
                                            n.chem      <- length(C)
                                    
                                            # split parameter vector per chemical
                                            par         <- split(par.drc,rep(1:ceiling(length(par.drc)/n.chem), each=n.chem))
                                            
                                            # caluclate toxic unit
                                            TU <- sapply(1:n.chem, function(i) C[i]/drc(x,par[[i]]))
                                   
                                            # calculate deviation from concentration addition model
                                            sum(TU)-1}

# synergism/antagonism model
#
#   input:  x       : effect size
#           C       : vector containing metal concentrations
#           drc     : type of dose response curve required (see 'inverse.drc.R' for available types)
#           par.drc : list of parameters for dose response curves for each chemical. Note that parameter values should be grouped per metal
#           par.mix : deviation parameter for the synergism/antagonism model 
#
#   output: deviation from the synergism/antagonism model

syn.ant    <- function(x,C,drc,par.drc,par.mix){  
  
                                            n.chem      <- length(C)
  
                                            # split parameter vector per chemical
                                            par         <- split(par.drc,rep(1:ceiling(length(par.drc)/n.chem), each=n.chem))
                                            
                                            # parameter mixture model
                                            a           <- unlist(par.mix)[1]
                                            
                                            # caluclate toxic unit
                                            TU <- sapply(1:n.chem, function(i) C[i]/drc(x,par[[i]]))
                                          
                                            # calculate deviation from synergism/antagonism model
                                            sum(TU)-exp(a*prod(TU/sum(TU)))}
                                           
  
# dose-ratio model
#
#   input:  x       : effect size
#           C       : vector containing metal concentrations
#           drc     : type of dose response curve required (see 'inverse.drc.R' for available types)
#           par.drc : list of parameters for dose response curves for each chemical. Note that parameter values should be grouped per metal
#           par.mix : deviation parameter for the synergism/antagonism model, the order should be first the same parameter as for the SA model (a), next the parameter for the dose-ratio model (b)
#
#   output: deviation from the synergism/antagonism model

dose.ratio  <- function(x,C,drc,par.drc,par.mix){  
                                            n.chem      <- length(C)
                                            
                                            # split parameter vector per chemical
                                            par         <- split(par.drc,rep(1:ceiling(length(par.drc)/n.chem), each=n.chem))
                                            
                                            # parameter mixture model
                                            a           <- unlist(par.mix)[1]
                                            b           <- unlist(par.mix)[2]
                                            
                                            # caluclate toxic unit
                                            TU <- sapply(1:n.chem, function(i) C[i]/drc(x,par[[i]]))
                                            
                                            # calculate deviation from dose-ratio model
                                            sum(TU)-exp((a+b*TU[2]/sum(TU))*prod(TU/sum(TU)))}
                                            
