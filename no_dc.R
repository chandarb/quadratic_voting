### solve for the optimal votes when there's no discontinuity in the voting function

solve_votes = function(sample_v1, vgrid_0, ugrid, sample_u, N){
  # vector of sums of votes
  sample_V = rowSums(sample_v1)
  
  # bounds for root solving
  upper = max(max(sample_V),max(sample_v1),max(vgrid_0),0.1)
  lower = min(min(sample_V),min(sample_v1),min(vgrid_0),-0.1)
  
  # estimate CDF of sum of votes using empirical CDF
  G = ecdf(sample_V)
  # kernel estimate for density of sum of votes
  g_d = density(sample_V)
  # function to estimate density at a point. uses linear interpolation.
  g_f = approxfun(g_d$x, g_d$y)
  g = function(x){
    d = g_f(x)
    d[is.na(d)] = 10^-12
    return(d)
  }
  
  # find optimal vote level for all values in the utility grid
  vgrid_1 = bisection(foc, lower=lower, upper=upper, u=ugrid, g=g)$mid
  
  # use mapping of utilities to votes to interpolate votes over sample utilities
  sample_vi = sample_votes(ugrid,vgrid_1,sample_u,N, NA)
  
  return(list(sample_vi=sample_vi, vgrid_1=vgrid_1))
}

### If N is sufficiently large use central limit theorem approximation instead of 
### kernel density estimator. Should only use if second moment of the distribution
### of votes exists.
solve_votes_bigN = function(sample_v1, vgrid_0, ugrid, sample_u, N){
  
  # mean sum
  mu_N = mean(rowSums(sample_v1))
  sd_N = sd(rowSums(sample_v1))
  
  upper = 5 * qnorm(1 - 10^-10, mu_N, sd_N)
  lower = 5 * qnorm(10^-10, mu_N, sd_N)
  
  # estimate CDF of sum of votes using empirical CDF
  G = function(x){return(pnorm(x, mu_N, sd_N))}
  # kernel estimate for density of sum of votes
  g = function(x){return(dnorm(x, mu_N, sd_N))}
  
  # find optimal vote level for all values in the utility grid
  vgrid_1 = bisection(foc, lower=lower, upper=upper, u=ugrid, g=g)$mid
  
  # use mapping of utilities to votes to interpolate votes over sample utilities
  sample_vi = sample_votes(ugrid,vgrid_1,sample_u,N, NA)
  
  return(list(sample_vi=sample_vi, vgrid_1=vgrid_1))
}

### This code is too slow. May revisit this in the future to add support
### for infinite variance distributions. uses generalized central limit
### theorem.
solve_votes_bigN_stable = function(sample_v1, vgrid_0, ugrid, sample_u, N){
  
  f = function(u, a, b, c, d){
    dstable(u, 2*exp(a)/(1+exp(a)), 2*exp(b)/(1+exp(b))-1, exp(c), d)
  }
  
  print("solving parameters")
  theta = fitdistr(sample_v1, f, list(a=1, b=0, c=log(mad(sample_v1)), 
                                      d=median(sample_v1)))$estimate
  print("done solving parameters")
  a = theta[1]
  b = theta[2]
  c = theta[3]
  d = theta[4]
  
  a_N = a
  b_N = b
  c_N = a * c
  d_N = a * d
  
  g = function(x){return(dstable(u, 2*exp(a_N)/(1+exp(a_N)), 
                                 2*exp(b_N)/(1+exp(b_N))-1, exp(c_N), d_N))}
  
  upper = qstable(1 - 10^-10, 2*exp(a_N)/(1+exp(a_N)), 
                  2*exp(b_N)/(1+exp(b_N))-1, exp(c_N), d_N)
  lower = qstable(10^-10, 2*exp(a_N)/(1+exp(a_N)), 
                2*exp(b_N)/(1+exp(b_N))-1, exp(c_N), d_N)
  print("solving votes")
  # find optimal vote level for all values in the utility grid
  vgrid_1 = bisection(foc, lower=lower, upper=upper, u=ugrid, g=g)$mid
  print("done solving votes")
  # use mapping of utilities to votes to interpolate votes over sample utilities
  sample_vi = sample_votes(ugrid,vgrid_1,sample_u,N, NA)
  
  return(list(sample_vi=sample_vi, vgrid_1=vgrid_1))
}

# compute average absolute distance between two grids
l1_dist = function(grid_1, grid_2){
  return(dist(rbind(grid_1, grid_2), method = "manhattan")/length(grid_1))
}

# Initialize to a linear voting function. Find the slope that differs least from the desired voting function.
init_slope = function(ugrid, sample_u, N){
  return(optimize(find_a_nodc, c(10^-3, 2), ugrid=ugrid, sample_u=sample_u, N=N)$minimum)
}

# Find the L1 distance between the current voting function and the desired one
find_a_nodc = function(a, ugrid, sample_u, N){
  # voting matrix
  sample_v1 = a * sample_u
  
  # mapping of utility values to votes
  vgrid_0 = a * ugrid
  
  # solve for the desired votes
  V = solve_votes(sample_v1, vgrid_0, ugrid, sample_u, N)
  
  sample_vi = V$sample_vi
  
  return(l1_dist(as.vector(sample_vi), as.vector(sample_v1)))
}