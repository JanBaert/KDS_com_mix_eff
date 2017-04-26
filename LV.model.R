# lotka volterra model: function is designed to be solved numerically using de 'ode' function from the 'deSolve' package

LV   <- function(Time, State, Pars) {  
                  with(as.list(c(State, Pars)), {
                            N  <- as.vector(State)   
                            dN <- sapply(c(1:length(State)),function(x) N[x]*mu[x]*(1-(interact[x,]%*%N)/K[x]))
                            return(list(dN)) })}
