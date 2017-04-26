# This folder contains a series of functions allowing simulating the Lotka-Volterra equations for various mixtures of environmental stressors

# Overview of files:
#                   - LV_model.R: the Lotka-Volterra equations, implemented as a function to be solved using the ode package
#                   - Mixture.models.R: File containing the mixture effect models by Jonker et al. 2005
#                   - inverse.drc.R: file containing inverse dose-response curves (only log-logistic at this point)
#                   - LV.solve.R: main function simulating the Lotka-Volterra equations
#                   - main_illustrion.R: dummy code illustrating the functions