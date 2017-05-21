rm(list=ls())

library(PearsonDS)
library(rootSolve)
library(rjson)
library(rmutil)
library(fBasics)
library(stabledist)
library(MASS)
source("dc_funcs.R")
source("welf_calc.R")
source("initialize.R")
source("auxilliary.R")
source("no_dc.R")
source("simulation.R")
set.seed(1234567)


######### run the simulation and gather results ###########
input = fromJSON(file="QV_input.json")

# Number of individuals in the population
N = input$N
# Utility distribution
distribution = input$distribution
# first parameter input
param1 = input$param1
# second parameter input
param2 = input$param2
# range of distribution (for Uniform and Beta distributions)
range=c(input$interval1, input$interval2)
# rate of Newton update
k=input$k
# how much error should be allowed before stopping
tolerance = input$tolerance

# Boolean for whether the program should decide whether to initialize with
# a discontinuity based on the limiting equation for when a discontinuity should
# arise
limiting=input$limiting
# if not initializing the discontinuity based on the limiting equation,
# manually specify whether there should be a discontinuity or not
with_dc=input$with_dc

#### run the simulation ####
vall = main.prog(distribution, para1=param1, para2=param2, range=range, N=N, k=k,
                 limiting=limiting, with_dc=with_dc, tolerance=tolerance)

# collect welfare loss
QVI = 1 - vall$QEW
MVI = 1 - vall$MEW
LI = vall$Limiting

# create CSV with welfare results

w_mat = matrix(0, ncol = 7, nrow = iterations)
colnames(w_mat) = c("distribution", "N", "param1", "param2", "QVI",
                    "MVI", "Limiting")
w_mat[,1] = distribution
w_mat[,2] = N
w_mat[,3] = param1
w_mat[,4] = param2
w_mat[,5] = QVI
w_mat[,6] = MVI
w_mat[,7] = LI

# filename = paste0(distribution,"N",N,"a",param1,"b",param2,".csv")
# write.csv(w_mat, filename)