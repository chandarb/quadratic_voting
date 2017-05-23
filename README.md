# quadratic_voting

Numerical Simulation of Quadratic Voting 

run_simulation.R: Take in the inputed parameters and run the simulation.

simulation.R: Code for the simulation.

auxilliary.R: Auxilliary functions. Root solver, utility function and first order condition for QV, utility to vote interpolator.

initialize.R: Initialize the simulation for the inputted utility distribution and parameters.

welf_calc.R: Compute welfare after a stable point is found.

dc_funcs.R: Functions relevant to cases when there's a discontinuity in the voting function.

no_dc.R: Functions relevant to cases when there's no discontinuity in the voting function.

QV_input.json: JSON file that can be used to input the parameter settings.


INPUTS:

distribution: the utility distribution that people's valuations are drawn from

param1: the first parameter to the utility distribution.
param2: the second parameter to the utility distribution.
  - Normal: param1 is mean, param2 is standard deviation
  - Beta: Pearson Type 1 Beta distribution shape parameters
  - DPLN-Uniform Mixture: param1 is the mean of the distribution, param2 is the tail parameter for the DPLN component
  - Laplace: location and scale parameters
  - Mixture of 2 DPLN distributions: tail parameters for the left and right sides of the distribution
  
  Other parameters are set in initialize.R
  
N: number of individuals in the population

interval1: lower bound for the distribution
interval2: upper bound
  - These parameters are only used for the Uniform and Beta cases.
  
k: rate of Newton update

limiting: Boolean for whether the program should decide whether to initialize with a discontinuity
based on the limiting equation for when a discontinuity should arise. Default is true.

with_dc: If not initializing with a discontinuity based on the limiting equation, manually specify whether there should be a discontinuity or not.

tolerance: how much error should be allowed before stopping

