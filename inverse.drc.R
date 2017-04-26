# Log-logistic dose-response curve
#
#   input:  x     : effect size [0,1]
#           par   : list containing parameters, list elements should be names "ec.50" and "beta"
#   output: ECx   : concentration at which effect size x occurs
    
LL.drc     <- function(x,par){                              
                      ec.50       <- par$ec.50
                      beta        <- par$beta
                      if(length(par)>2){max <- par$max}else{max <- 1}
                      if(length(par)>3){min <- par$min}else{min <- 0}
                      ec.50*((max-min)/(x-min)-1)^(1/beta)}
