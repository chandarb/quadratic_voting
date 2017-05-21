rm(list=ls())

setwd("/Users/bharatchandar/Documents/QV/")
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
set.seed(1234567)

# Run the simulation.
# K is the fineness of the grid, S_size is the number of rows in sampled utility matrix,
# k controls the rate of convergence
main.prog = function(distribution, para1 = NULL, para2 = NULL, 
                     N, range = NULL, K=1000, S_size=NULL, k=.01, limiting=TRUE, 
                     with_dc=FALSE, tolerance=0.05){
  
  if(distribution != "Normal" & distribution != "Beta" & distribution != 
       "Uniform" & distribution != "Mar" & distribution != "dPLN" & 
       distribution != "uDPLN" & distribution != "DPLNu" &
        distribution != "Laplace" & distribution != "DPLN2")
    stop('The program does not support this distribution')
  
  # election sample size
  if(is.null(S_size)){
    S_size = trunc(5 * 10^6 / N)
  }

  # get distribution-specific parameters
  init_vals = initialize(distribution, para1 = para1, para2 = para2,
                        N, range = range, K=K, S_size=S_size)

  # mean of utility distribution
  mu = init_vals$mu
  # standard deviation of utility distribution
  stdv = init_vals$stdv
  # bounds for utility grid
  interval = init_vals$interval
  # sample of values from utility distribution
  sample_uN = init_vals$sample_uN
  # density for utility values
  f = init_vals$f
  # CDF for utility values
  F = init_vals$F
  # Quantile function for utility values
  Qt = init_vals$Qt
  # Sample from the distribution
  sampler = init_vals$sampler
  # Limiting inefficiency
  Limiting_i = init_vals$limiting
  
  # Define the matrix of sample utilities for N-1 voters
  sample_u = matrix(as.vector(sample_uN)[-(1:S_size)], ncol = N-1)
  
  # grid of possible utility values
  ugrid = seq(interval[1], interval[2], length = K)  
  
  # negative utilities in the utility grid
  uneg = ugrid[ugrid < 0]
  
  # number of negative values
  nneg = length(uneg)
  
  # constant used to set bounds on optimization and root solving problems
  m = 6
  
  # discontinuity in voting function
  dc0 = NA
  dc1 = NA
  
  # The value of N for which a discontinuity should likely appear
  limiting_val = 1 + (4 * stdv^2 / mu^2)
  
  # if this is greater than 2 then there should be a discontinuity
  dc_cutoff = sqrt(N-1) * mu / stdv
  if (dc_cutoff > 1 & dc_cutoff < 3){
    warning("N is close to the limiting cutoff for when a discontinuity should appear")
  }
  
  # decide whether to initialize with a discontinuity
  dc_init = FALSE
  # not using the limiting approximation for whether there's a discontinuity
  if (!limiting){
    # user wants to initialize with a discontinuity
    if (with_dc)
      dc_init = TRUE
  }
  # using the limiting approximation
  else{
    # above the limiting cutoff - initialize with discontinuity
    if (dc_cutoff > 2)
      dc_init = TRUE
  }
  
  #### Initializing with a discontinuity in the voting function ####
  if (dc_init){
    # smoothing constant used during Newton update
    A = exp(-1) + exp(-2) + exp(-3) + exp(-4)
    
    # lower bound for where discontinuity in voting function
    # should exist
    initl = interval[1]
    
    # for dpln may need to look further into the tail for the 
    # discontinuous voter
    if (distribution=="dPLN"){
      initl = initl * m
    }
      
    # initialize the discontinuity at the limiting value
    dc0 = init_dc(F, mu, N, lower=initl)
    
    # get initial voting grid. solves for linear moderate voting function
    # that leads to discontinuity at dc0.
    init = ext_init(ugrid,N,dc0,F, sample_u, m, distribution=distribution)
    
    # initial voting grid
    vgrid_0 = init$vgrid
    
    
    # initial sample voting matrix
    sample_v1 = sample_votes(ugrid, vgrid_0, sample_u, N, dc0)
    
    
    # vector of desired dc's
    u_des = c(dc0)
    # vector of true dc's
    dcs = c(dc0)
    
    # density of moderate votes
    mods = mod_dens(sample_u, sample_v1, dc0, N)
    
    # vector of desired means
    mu_des = c(mods$mu)
    # vector of true means
    mus = c(mods$mu)
    
    # numerical derivatives for discontinuous utility level
    # with respect to mean of moderate vote sums
    U_ps = c()
    # derivatives of moderate vote sums with respect to 
    # discontinuity
    Mu_ps = c()
    
    # matrix of derivatives of objective function
    J = matrix(c(0,0,0,0), ncol=2)
    
    # initial time step
    t = 2
    
    # initialize convergence criteria
    error = 1000
    mu_error = 1000
    errors = c()
    mu_errors = c()
    max_iter = 500 / k

    while(((abs(error) > tolerance) | (abs(mu_error) > tolerance) | (t < 10)) & t < max_iter){
      
      # discontinuity last period
      prev_dc = dcs[t-1]
      # mean of sum of moderate votes last period
      prev_mu = mus[t-1]
      # extremist mass
      p=F(prev_dc)
      
      mods = mod_dens(sample_u, sample_v1, prev_dc, N)
      # density with only moderates
      ge = mods$g
      # CDF with only moderates
      Ge = mods$G
      # largest absolute sum of votes in sample
      max_V = mods$max_V
      # smeared density
      gsmear = make_smear(ugrid,vgrid_0, sample_u, sample_v1, prev_dc, F, distribution, prev_dc)
      # smeared CDF
      Gsmear = CDF_smear(ugrid,vgrid_0, sample_u, sample_v1, prev_dc, F, distribution, prev_dc)
      
      # full density - mixture
      g = function(v){return((1 - ((N-1)*p))*ge(v) + (N-1)*p*gsmear(v))}
      G = function(v){return((1 - ((N-1)*p))*Ge(v) + (N-1)*p*Gsmear(v))}
            
      # new desired discontinuity
      dc1 = exact_dc(g, G, m * interval[1], 0, -m * max_V, 0)
      
      # find the roots at the desired discontinuity
      z_d = uniroot.all(foc, lower=-m*max_V,upper=0, u=dc1, g=g)
      z_d = sort(z_d)
      
      # solve for vote levels at the previous discontinuity
      z <- uniroot.all(foc,lower=-m*max_V,upper=0, u=prev_dc, g=g)
      
      # moderate utility grid
      ugmod = ugrid[ugrid > prev_dc]
      
      # extremist utility grid
      ugext = ugrid[ugrid <= prev_dc]
      
      ### If there's now only one root at last quarter's discontinuity,
      ### to solve for the mean of moderate votes given last period's
      ### discontinuity have to approximate the third root. 
      if (length(z)==1){
        warning("Previous DC only had one root")
        # the previous discontinuity now only has an extremist root
        if (z < z_d[1]){
          # utility values greater than the desired dc
          ugmod1 = ugrid[ugrid > dc1]
          # corresponding moderate vote levels
          vgmod1 = bisection(foc, lower=z_d[2], upper=m*max_V, u=ugmod1, g=g)$mid
          # part of the utility grid between the previous dc and the desired one
          ugbet = ugrid[(ugrid > prev_dc) & (ugrid <= dc1)]
          # fit a linear model over the moderate vote levels closest to the dc
          votes = vgmod1[1:6]
          utils = ugmod1[1:6]
          fit = lm(votes ~ utils, data.frame(utils, votes))
          # extrapolate over the stretch where there might not be three roots
          # to get extremist votes
          vgbet = predict(fit, data.frame(utils=ugbet))
          # estimated moderate voting grid
          vgmod = c(vgbet, vgmod1)
          
          # extremist voting grid
          vgext = bisection(foc, lower=-m*max_V, upper=z, u=ugext, g=g)$mid
        }
        # the previous discontinuity only has a moderate root
        else {

          # moderate voting grid
          vgmod = bisection(foc, lower=z, upper=m*max_V, u=ugmod, g=g)$mid
          
          # utility values less than the desired dc
          ugext1 = ugrid[ugrid <= dc1]
          # corresponding extremist vote levels
          vgext1 = bisection(foc, lower=-m*max_V, upper=z_d[1], u=ugext1, g=g)$mid
          # part of the utility grid between the previous dc and the desired one
          ugbet = ugrid[(ugrid > dc1) & (ugrid <= prev_dc)]
          
          # number of points to interpolate over
          ipoints = min(length(ugext), 7) - 1
          if (ipoints==0){
            vgbet = rep(vgext1, length(ugbet))
          }
          else{
            # fit a linear model over the extremist vote levels closest to the dc
            votes = vgext1[(length(ugext1)-ipoints):length(ugext1)]
            utils = ugext1[(length(ugext1)-ipoints):length(ugext1)]
            fit = lm(votes ~ utils, data.frame(utils, votes))
            # extrapolate over the stretch where there might not be three roots
            # to get extremist votes
            vgbet = predict(fit, data.frame(utils=ugbet))
          }
          # estimated extremist voting grid
          vgext = c(vgext1, vgbet)
        }
      }
      
      ### there are three roots at last period's discontinuity
      else{
        z = sort(z)
        # upper bound on extremist votes
        v1 = z[1]
        # lower bound on moderate votes
        v2 = z[3]
        # extremist voting grid
        vgext = bisection(foc, lower=-m*max_V, upper=v1, u=ugext, g=g)$mid
        # moderate voting grid
        vgmod = bisection(foc, lower=v2, upper=m*max_V, u=ugmod, g=g)$mid
      }
      
      # full desired voting grid
      vgrid_1 = c(vgext, vgmod)
      
      # desired sample voting matrix
      sample_vd = sample_votes(ugrid, vgrid_1, sample_u, N, prev_dc)
      # desired mean of sum of moderate votes
      mud = mod_dens(sample_u, sample_vd, prev_dc, N)$mu
      
      # to populate the vector of previous dc and mu
      # values on the first iteration, just add
      # a slight perturbation
      if (t < 3){
        # change the mass very slightly
        new_mass = p * (1 - (k / 10))
        
        # if numerically the utility values 
        # are equivalent increase the perturbation more
        while (abs(prev_dc - Qt(new_mass)) == 0){
          new_mass = (1 - (k / 10)) * new_mass
        }
        # new discontinuity
        dct = Qt(new_mass)
        # append to vectors
        dcs = c(dcs, dct)
        u_des = c(u_des, dct)
        
        # change mean very slightly
        mut = prev_mu * (1 - (k / 10^2))
        # append to vectors
        mus = c(mus, mut)
        mu_des = c(mu_des, mut)
        
      }
      
      else{
        # update vectors of desired values
        mu_des = c(mu_des, mud)
        u_des = c(u_des, dc1)
        
        ### Newton update
        
        mu_diff = mus[t-1] - mus[t-2]
        l_up = length(U_ps)
        
        # on the first few iterations just take a
        # standard derivative
        if (l_up < 3){
          U_p = (u_des[t] - u_des[t-1]) / mu_diff
        }
        
        else if (abs(mu_diff) > 10-10){
          # derivative from this period
          U_val = (u_des[t] - u_des[t-1]) / mu_diff
          # take weighted average of derivatives from previous periods
          # to increase stability. also shrinks derivative towards 0.
          U_p = (exp(-1)/A) * U_val + (exp(-2)/A) * U_ps[l_up] + 
            (exp(-3)/A)* U_ps[l_up - 1] 
        }
        
        # if the absolute value of the denominator is less than 
        # 10^-10 set it to 10^-10.
        # This makes it more numerically stable and is a conservative
        # estimate of the derivative.
        else {
          print("Mu_dem 0")
          # sign of the denominator
          sden = sign(mu_diff)
          if (sden == 0){
            sden = 1
          }
          U_val = sden * (u_des[t] - u_des[t-1]) / 10^-10
          U_p = (exp(-1)/A) * U_val + (exp(-2)/A) * U_ps[l_up] + 
            (exp(-3)/A)* U_ps[l_up - 1] 
        }
        
        # store the derivatives
        U_ps = c(U_ps, U_p)
        
        dc_diff = (dcs[t-1] - dcs[t-2])
        l_mup = length(Mu_ps)

        # same thing for derivative of the mean
        if (l_mup < 3)
          Mu_p = (mu_des[t] - mu_des[t-1]) / dc_diff

        else if (abs(dc_diff) >= 10^-10){
          Mu_val = (mu_des[t] - mu_des[t-1]) / dc_diff
          Mu_p = (exp(-1)/A) * Mu_val + (exp(-2)/A) * Mu_ps[l_mup] + 
            (exp(-3)/A)* Mu_ps[l_mup - 1] 
        }
        else{          
          print("Mu_dem 0")
          # sign of the denominator
          sden = sign(dc_diff)
          if (sden == 0){
            sden = 1
          }
          Mu_val = sden * (mu_des[t] - mu_des[t-1]) / 10^-10
          Mu_p = (exp(-1)/A) * Mu_val + (exp(-2)/A) * Mu_ps[l_mup] + 
            (exp(-3)/A)* Mu_ps[l_mup - 1] 
        }

        Mu_ps = c(Mu_ps, Mu_p)
        
        # change in desired mean. first component in the objective function
        # that we're trying to find the 0's of.
        del_mu = mus[t-1] - mu_des[t]
        # change in desired discontinuity. second component of objective
        # function. 
        del_u = dcs[t-1] - u_des[t]
        # objective function
        df = c(del_u, del_mu)

        # Jacobian of objective
        J = matrix(c(1, -Mu_p, -U_p, 1), ncol=2)
        
        # solve for update
        dp = solve(J, -k * df)
        
        # get new values for this period
        dct = dcs[t-1] + dp[1]
        mut = mus[t-1] + dp[2]
        
        dcs = c(dcs, dct)
        mus = c(mus, mut)
        
      }
      
      ### need to readjust last period's vgrid to account for the new dc
      if (dct > prev_dc){
        # part of utility grid between the previous and new dc
        ugext1 = ugrid[(ugrid > prev_dc) & (ugrid <= dct)]
        # find new upper bound for extremist votes
        z1 <- uniroot.all(foc,lower=-m*max_V,upper=0, u=dct, g=g)
        z1 = sort(z)

        vgext1 = c()
        # make the extremist portion longer
        if (length(ugext1) != 0){
          # if there's only a moderate root then repeat the values of 
          # highest extremist utility in the utility grid.
          if (length(z1)==1){
            vgext1 = rep(vgext[length(vgext)], length=length(ugext1))
          }
          # otherwise solve
          else{
            vgext1 = bisection(foc, lower=v1, upper=z1[2], u=ugext1, g=g)$mid
          }
        }
        
        # new voting grid
        vgrid_1 = c(vgext[(ugext <= prev_dc)], vgext1, vgmod[ugmod > dct])
        
        # adjust old voting grid to account for new dc - just repeat
        # last extremist vote level for stretch with no extremist vote
        vgrid_0 = c(vgrid_0[ugrid < prev_dc],
                    rep(vgrid_0[length(ugrid[ugrid <= prev_dc])],
                        length=length(ugrid[(ugrid > prev_dc) & (ugrid <= dct)])),
                    vgrid_0[ugrid > dct])
      }
      
      else if (dct <= prev_dc){
        # part of utility grid between the previous and new dc
        ugmod1 = ugrid[(ugrid > dct) & (ugrid <= prev_dc)]
        # find new lower bound for moderate votes
        z1 <- uniroot.all(foc,lower=-m*max_V,upper=0, u=dct, g=g)
        z1 = sort(z)

        vgmod1 = c()
        # make the moderate portion longer
        if (length(ugmod1) != 0){
          if (length(z1)==1){
            vgmod1 = rep(vgmod[1], length(ugmod1))
          }
          else{
            vgmod1 = bisection(foc, lower=z1[3], upper=v2, u=ugmod1, g=g)$mid
          }
        }
        
        # new voting grid
        vgrid_1 = c(vgext[ugext <= dct], vgmod1, vgmod[ugmod > prev_dc])
        # adjust old voting grid to account for new dc - just repeat
        # last moderate vote level for stretch with no moderate vote
        vgrid_0 = c(vgrid_0[ugrid <= dct],
                    rep(vgrid_0[ugrid > prev_dc][1],
                        length=length(ugmod1)),
                    vgrid_0[ugrid > prev_dc])
      }
      
      # constant to rescale the voting grid and get the desired mean
      del = (mut - prev_mu) / (mud - prev_mu)
      
      # rescale the voting grid to get the right mean
      vgrid_0 = vgrid_0 + del * (vgrid_1 - vgrid_0)
      
      # new sampled votes
      sample_v1 = sample_votes(ugrid, vgrid_0, sample_u, N, dct)
      
      t = t+1
      # error in discontinuity mass
      error = (F(dct) - F(dc1)) / (F(dct))
      # error in mean of moderate votes
      mu_error = (mut - mud) / (mut)
      errors = c(errors, error)
      mu_errors = c(mu_errors, mu_error)
      # decrease the step size every so often if convergence stops
      if (t %% (500) == 0){
        qss = quantile(errors[(length(errors)-100):length(errors)], 
                 probs=c(.25, .75))
        if ((error > qss[1]) & (error < qss[2])){
          k = k / 2
          print("Reduction in step size")
          print(paste0("k:", k))
        }
      }
      print(sprintf("t: %d dc error: %f  mu error: %f", t, error, mu_error))
      print(sprintf("mut:%f mud:%f dct:%f dc1:%f", mut, mud, dct, dc1))
    }
    
  }
  
  #### Initializing without a discontinuity in the voting function ####
  else{
  
    # initialize the voting function
    a = init_slope(ugrid, sample_u, N)
    # intiialize voting grid
    vgrid_0 = ugrid * a
    # initial voting sample
    sample_v1 = sample_u * a
    # initial actual voting function - to be used in Newton updating
    sample_v0 = sample_v1
    # initial desired voting function
    sample_vi0 = sample_v1
  
    t = 0
    
    # initialize convergence criteria
    error = 100000

    # Newton update
    c=1

    # error is the euclidean distance between the voting grids
    # on successive iterations    
    errors = c()
    
    max_iter = 500 * (1/k)
    

    while((error>tolerance & t < max_iter)){
      
      # solve for optimal vote levels
      if (N < 10000){
        V = solve_votes(sample_v1, vgrid_0, ugrid, sample_u, N)
      }
      # if N is sufficiently large use Central Limit Theorem 
      # approximation to estimate the distribution of the sum of 
      # votes. this should only be used if the first two moments exist.
      else{
        V = solve_votes_bigN(sample_v1, vgrid_0, ugrid, sample_u, N)
      }
      
      # desired sample voting matrix
      sample_vi = V$sample_vi
      
      # desired grid of votes
      vgrid_1 = V$vgrid_1
      
      if (t > 0){
        # mean l1 distance between actual votes
        d_actual = l1_dist(as.vector(sample_v0), as.vector(sample_v1))
        # mean l1 distance between desired votes
        d_desired = l1_dist(as.vector(sample_vi0), as.vector(sample_vi))
        # ratio
        c = d_actual/d_desired
      }
      
      # update last period's actual voting function
      sample_v0 = sample_v1
      # update last period's desired voting function
      sample_vi0 = sample_vi
      # limit movement in sample votes between periods
      sample_v1 = c*k*sample_vi + (1-(k*c))*sample_v1
      error = l1_dist(as.vector(sample_v1), as.vector(sample_vi))
      print(sprintf("error=%f",error))
      # update voting grid
      vgrid_0 = vgrid_1
      
      # continue to next period
      t = t+1
      print(sprintf("t = %d",t))
      errors = c(errors, error)
      
      # change rate of update if convergence stops
      if (t %% (500) == 0){
        qss = quantile(errors[(length(errors)-100):length(errors)], 
                       probs=c(.25, .75))
        print(qss)
        if ((error > qss[1]) & (error < qss[2])){
          k = k / 2
          print("Reduction in step size")
          print(paste0("k:", k))
        }
      }
  
    }
  }
  if (error > tolerance){
    warning("did not converge")
  }
  # Calculate all the votes for sample of size N
  sample_v2 = sample_votes(ugrid,vgrid_1,sample_uN,N+1, dc1)
  # sum across elections
  sample_V  = rowSums(sample_v2)
  
  if (!dc_init)
    dct = NA

  # sums of utilities
  U <- rowSums(sample_uN)
  # sums of votes
  V <- sample_V
  
  if (!is.na(dct)){
    # if the distribution is fat-tailed 
    if (distribution %in% c("DPLNu", "dPLN")){
      # create a grid of extremist values
      fineness = 10^4
      p_fine = seq(F(dct)/fineness, F(dct), length=fineness)
      fine_points = Qt(p_fine)
      # solve for the vote levels at the extremist values
      fine_ext = bisection(foc, lower=-m*max_V, upper=v1, u=fine_points, g=g)$mid
      # compute welfare by taking expectation over possible extremist vote levels
      ww = 1 - welf_fine(sampler, dct, F, ugrid, vgrid_0, 5 * S_size, N, p_fine, fine_ext)
    }
    else{
      # take numerical expectation of welfare
      ww = welf_dc(sampler, dct, F, ugrid, vgrid_0, 5 * S_size, N)
    }
    print("welfare")
    print(ww)
    QEW = ww
  }
  else{
    # different function to compute welfare if there's no discontinuity
    QEW <- welf(U, V)
  }
  # majority voting welfare
  MEW <- welf_m(U, sample_uN)
  print(sprintf("QV: %f, Majority: %f",QEW,MEW))

  # Returns utility grid(ugrid), vote function given utility(vgrid_1),
  # the name of the distribution and its parameters(distribution, para1,para2,range), the number of individuals(N),
  # the vote totals for each election (sample_V), the sample utility matrix (sample_u),
  # the welfare for QV (QEW), and the welfare for 1p1v (MEW).
  
  return(list(ugrid = ugrid, vgrid_0 = vgrid_0, distribution = distribution,para1 = para1, 
              para2 = para2, N = N, range = range, sample_V = sample_V,sample_u=sample_uN, 
              QEW=QEW, MEW=MEW, Limiting=Limiting_i, dct=dct))
}


######### run the simulation and gather results ###########
input = fromJSON(file="QV_input.json")

QVI = 0
MVI = 0
LI = 0

N = input$N
param1 = input$param1
param2 = input$param2
range=c(input$interval1, input$interval2)
k=input$k
limiting=input$limiting
with_dc=input$with_dc
distribution = input$distribution
tolerance = input$tolerance
print(N)

vall = main.prog(distribution, para1=param1, para2=param2, range=range, N=N, k=k,
          limiting=limiting, with_dc=with_dc, tolerance=tolerance)

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

# filename = paste0("data/DPLNu_seed/",distribution,"N",N,"a",param1,"b",param2,".csv")
# write.csv(w_mat, filename)
