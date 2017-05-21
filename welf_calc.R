# functions to compute welfare


# compute welfare
# for QV pass in sums of votes
# for 1p1v call welf_m
welf = function(U, V){
  abs_U = abs(U)
  # vector of whether the optimal outcome is reached in each election
  signs = sign(U) == sign(V)
  # total utility for correct outcomes divided by total utility across all elections
  return((t(signs) %*% abs_U)[1, 1] / sum(abs_U))
}

# "flips a coin" if there's a tie in majority voting
# used in welf_m
zero = function(x){
  if (x == 0){
    return(sample(c(-1,1),1))
  }
  return(x)
}

# compute welfare of 1p1v
welf_m = function(U, sample_u){
  # whether each voter is for or against
  samp_s = sign(sample_u)
  # net votes for each election
  U_m = rowSums(samp_s)
  # decide randomly if equal votes
  U_m = sapply(U_m, zero)
  return(welf(U, U_m))
}

# compute welfare when there's a discontinuity
welf_dc = function(sampler, dc, F, ugrid, vgrid, S_size, N){
  # extremist mass
  p = F(dc)
  # sample moderates
  mods = matrix(sampler((N)*S_size, p=p, lower.tail=FALSE), nrow=S_size)
  # sample from the extremists
  exts = sampler(S_size, p=p, lower.tail=TRUE)
  # join them together to make a matrix of samples with one extremist
  sampleu = cbind(exts, mods[,1:(N-1)])
  
  # votes for moderates
  intermod = sample_votes(ugrid, vgrid, mods, N+1, dc)
  # votes for extremists
  interext = sample_votes(ugrid, vgrid, exts, 2, dc)
  # joined together
  interjoin = cbind(interext, intermod[,1:(N-1)])
  # sum of utility in each sample
  U1 = rowSums(sampleu)
  # sum of votes in each sample
  V1 = rowSums(interjoin)
  # compute welfare conditional on an extremist
  QEW1 = welf(U1, V1)
  
  # compute welfare conditional on no extremist
  U2 = rowSums(mods)
  V2 = rowSums(intermod)
  QEW2 = welf(U2, V2)
  
  # take weighted average based on probability of extremist
  return((N*p*(QEW1)) + ((1 - (N*p)) * (QEW2)))
}

# for DPLN distributions assume normality by central limit theorem
# and take a numerical integral to make it less
# sample dependent
welf_fine = function(sampler, dc, F, ugrid, vgrid, 
                     S_size, N, p_fine, fine_ext){
  # extremist mass
  p = F(dc)
  # sample moderate votes
  mods = sample_votes(ugrid, vgrid, 
                      sampler((N)*S_size, p=p, 
                              lower.tail=FALSE), N+1, dc)
  
  # mean
  mu_m = mean(mods) * N
  # standard deviation
  sd_m = sd(mods) * sqrt(N)
  
  # moderate inefficiency
  mod_comp = (1 - N*p) * pnorm(0, mu_m, sd_m)

  mu_m1 = mu_m * (N-1) / N
  sd_m1 = sd_m * sqrt(N-1) / sqrt(N)
  
  # for each point in the fine sample of extremists 
  # compute expected inefficiency
  ff = diff(c(0, p_fine))
  # inefficency at each point
  fine_int = pnorm(-fine_ext, mu_m1, sd_m1)
  # weighted average inefficiency
  ext_comp = N * sum(fine_int * ff)
  
  # sum of moderates only and extremist inefficiency
  return(mod_comp + ext_comp)
}
