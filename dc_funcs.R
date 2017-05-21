# functions used when there's a discontinuity in the voting function


# solve abs(u) / F(u) = mu * (N-1)^2
# for u to initialize the discontinuity
limiting_dc = function(u, F, c){
  return(-u - c*F(u))
}

# finds the initial discontinuity point
init_dc = function(F, mu, N, lower){
  const = mu * (N-1)^2
  return(uniroot(limiting_dc, lower=lower, upper=0, F=F, c=const)$root)
}

# estimate density and CDF of sum of moderate votes
# also returns the maximum absolute vote total, used
# to determine a lower bound on possible extremist vote levels
mod_dens = function(sample_u, sample_v, dc, N){
  # identify moderates
  mods = sample_u > dc
  # number of moderates for each election
  num_mods = rowSums(mods)
  # keep only moderate votes
  mod_votes = sample_v * mods
  # sum moderate votes
  mod_sums = rowSums(mod_votes)
  # rescale as if there are N-1 voters
  rmod_sums = (N-1) * mod_sums / num_mods
  # construct density estimate
  g_d = density(rmod_sums)
  # functional approximation
  g_f = approxfun(g_d$x, g_d$y)
  # deal with NA cases
  g = function(x){
    d = g_f(x)
    d[is.na(d)] = 10^-12
    return(d)
  }

  # numerical integral of density
  G = make_G(g_d$x, g_d$y)

  return(list(g=g, G=G, max_V=max(abs(rmod_sums)), mu=mean(rmod_sums)))
}

# since the density is a linear interpolation
# the CDF is the sum of trapezoids
make_G = function(x, y){
  # base of trapezoid
  xdiff = diff(x)
  ydiff = diff(y)
  # rectangular part of trapezoid
  riemann1 = y[1:(length(y)-1)] * xdiff
  # triangular part
  triangle = ydiff * xdiff / 2
  # area of each trapezoid
  area = riemann1 + triangle
  
  # CDF evaluated at each sampled point
  # add a point at the beginning for mass outside sample
  total_area = c(10^-12, cumsum(area))
  
  a = x[1]
  n = length(x)
  b = x[n]
  sl = (n - 1) / (b - a)
  # linear map to index in estimated density quantiles
  cst = 1 - (sl * a)
  
  # CDF at an arbitrary point
  G = function(z){
    # apply linear map to values in sample_u to get the corresponding utility index
    # in ugrid
    inds = z * sl + cst
    # for points outside of the bounds just set to lowest val or highest val
    inds[inds < 1] = 1
    inds[inds > n] = n
    # find the closest index less than
    indsf = floor(inds)
    # corresponding density estimate at that point
    lowy = y[indsf]
    # closest index greater than
    indsc = indsf+1
    indsc[indsc > n] = n
    hiy = y[indsc]
    
    # y value is weighted average
    prop = inds - indsf
    
    yz = lowy * (1 - prop) + prop * hiy
    
    # values at closest point to the left
    yf = y[indsf]
    xf = x[indsf]
    
    # triangle part of trapezoid at point
    triangle = (z - xf) * (yz - yf) / 2
    
    # add base
    trapezoid = triangle + (z - xf) * yf
    
    # add mass from values less than
    cdf_val = trapezoid + total_area[indsf]
    
    return(cdf_val)
  }
  
  return(G)
}

# Computes a "smeared" density for moderates (accounting for the possibility of extremists)
make_smear = function(ugrid,vgrid, sample_u, sample_v, dc0, F, distribution, dc){
  # bounded utility distributions
  if (distribution == "Beta" | distribution == "Uniform" | distribution == "Mar" |
        distribution == "uDPLN"){
    # there are N-2 moderates and 1 extremist
    gz = mod_dens(sample_u, sample_v, dc0, N-1)$g
    # assign the extremist the same vote level as for the minimum utility
    gsmear = function(v){
      return(gz(v - vgrid[1]))
    }
  }
  # unbounded densities. estimate an "expected" density given an extremist
  # by numerically integrating over possible extremist vote levels in the tail.
  else{
    gz = mod_dens(sample_u, sample_v, dc0, N-1)$g
    # extremist grid
    uext = ugrid[ugrid <= dc0]
    # fineness of estimation
    fine = 25
    # finely grid the extremist mass to numerically estimate the integral
    usmear = seq(min(ugrid),dc0,length=fine)
    
    ### compute the following mixture
    ### gn(v) = p1 g(v + e1) + ... + pM g(v+eM)
    
    # interpolate vote level for each utility value in the smear
    v_ext = as.vector(sample_votes(uext, vgrid[ugrid<dc0], usmear, 2, dc))
    # evaluate CDF for each utility value
    F_vals = sapply(usmear, F)
    # find mass corresponding to utility values
    u_mass = c(F_vals[1], diff(F_vals))
    # mixture proportion is share of mass
    props = u_mass / F_vals[fine]
    
    gsmear = function(v){
      tmp = 0
      for(i in 1:length(usmear)){
        tmp = tmp + gz(v - v_ext[i])*props[i]
      }
      return(tmp)
    }
  }
  return(gsmear)
}

# Computes a "smeared" CDF for moderates (accounting for the possibility of extremists)
CDF_smear = function(ugrid,vgrid, sample_u, sample_v, dc0, F, distribution, dc){
  # bounded utility distributions
  if (distribution == "Beta" | distribution == "Uniform" | distribution == "Mar" |
        distribution == "uDPLN"){
    # there are N-2 moderates and 1 extremist
    Gz = mod_dens(sample_u, sample_v, dc0, N-1)$G
    # assign the extremist the same vote level as for the minimum utility
    Gsmear = function(v){
      return(Gz(v - vgrid[1]))
    }
  }
  # unbounded densities. estimate an "expected" density given an extremist
  # by numerically integrating over possible extremist vote levels in the tail.
  else{
    Gz = mod_dens(sample_u, sample_v, dc0, N-1)$G
    # extremist grid
    uext = ugrid[ugrid <= dc0]
    # fineness of estimation
    fine = 25
    # finely grid the extremist mass to numerically estimate the integral
    usmear = seq(min(ugrid),dc0,length=fine)
    
    ### compute the following mixture
    ### Gn(v) = p1 G(v + e1) + ... + pM G(v+eM)
    
    # interpolate vote level for each utility value in the smear
    v_ext = as.vector(sample_votes(uext, vgrid[ugrid<dc0], usmear, 2, dc))
    # evaluate CDF for each utility value
    F_vals = sapply(usmear, F)
    # find mass corresponding to utility values
    u_mass = c(F_vals[1], diff(F_vals))
    # mixture proportion is share of mass
    props = u_mass / F_vals[fine]
    
    Gsmear = function(v){
      tmp = 0
      for(i in 1:length(usmear)){
        tmp = tmp + Gz(v - v_ext[i])*props[i]
      }
      return(tmp)
    }
  }
  return(Gsmear)
}

# function that computes the difference in welfare between acting as a 
# moderate and extremist for an individual with the discontinuous 
# utility level for a given choice of slope for the linear moderate
# voting function
find_a = function(at, u, ugrid, sample_u, N, p, ugext, m, distribution, F){
  # sample moderate voting matrix - ignore extremists
  sample_v = sample_votes(ugrid, at * ugrid, sample_u, N, u)
  
  mods = mod_dens(sample_u, sample_v, u, N)
  # density for sum of moderate votes
  ge = mods$g
  # CDF for sum of moderate votes
  Ge = mods$G
  # largest absolute sum of votes in sample
  max_V = mods$max_V
  # rescale - this is a mixture component for the full density
  gt = function(v){return((1-(N-1)*p)*ge(v))}
  
  # solve for the extremist vote level at the discontinuity
  z <- uniroot.all(foc,lower=-m*max_V,upper=0, u=u, g=gt)
  z = sort(z)
  # should be extremist vote level
  v = z[1]
  
  # only one root
  if (length(z)==1){
    # dc doesn't have a moderate root
    # moderates need to vote more
    if(Ge(-v) > .5){
      return(Inf)
    }
    # dc doesn't have an extremist root
    # moderates should vote less
    return(-Inf)
  }
  
  # dc has multiple roots
  
  # extremist voting grid   
  v_ext = bisection(foc, lower=-m*max_V, upper=v, u=ugext, g=gt)$mid
  # mixture component without extremist
  Gt = function(v){return((1-(N-1)*p)*Ge(v))}
  # mixture components with extremists
  gsmeart = make_smear(ugext,v_ext, sample_u, sample_v, u, F, distribution, dc=u)
  Gsmeart = CDF_smear(ugext,v_ext, sample_u, sample_v, u, F, distribution, dc=u)
  # full mixtures
  gsmt=function(v){return(gt(v) + (N-1)*p*gsmeart(v))}
  Gsm = function(v){return(Gt(v) + (N-1)*p*Gsmeart(v))}
  
  # solve for optimal moderate vote level for discontinuous voter
  # should be closer to 0 than the extremist vote
  zm = uniroot.all(foc,lower=-m*max_V,upper=0, u=u, g=gsmt)
  
  # check again if there's a moderate root
  if (length(zm) == 1)
    return(Inf)
  
  # moderate welfare at dc
  wm=util.fn(zm[3],u,Gsm)
  # extremist welfare at dc
  we=util.fn(zm[1],u,Gsm)
  
  return(we - wm)  
}

# function to find the slope of the initial linear moderate voting function
# that induces a discontinuity at the initial discontinuity point. 
# uses a root solver to find a slope where the moderate welfare equals
# the extremist welfare at the initial discontnuity point.
ext_init = function(ugrid,N,u,F, sample_u, m, distribution){

  # mass of discontinuity
  p = F(u)

  # extremist grid
  ugext = ugrid[ugrid <= u]
  
  # find the moderate slope that induces a discontinuity at u
  # need to first find good bounds for the root solver
  
  # the lower bound is too large
  # moderates are voting too much
  low = 10^-2
  while (find_a(low, u=u, ugrid=ugrid, sample_u=sample_u, 
  N=N, p=p, ugext=ugext, m=m, distribution=distribution, F=F) < 0){
    low = low / 10
  }
  # upper bound
  high = low * 1.05
  while (find_a(high, u=u, ugrid=ugrid, sample_u=sample_u, 
                N=N, p=p, ugext=ugext, m=m, distribution=distribution, F=F) > 0){
    high = high * 1.05
  }

  low = high / 1.05

  # after finding good bounds use root solver to find precise value for slope
  un = uniroot(find_a, lower=low, upper=high, u=u, ugrid=ugrid, sample_u=sample_u, 
               N=N, p=p, ugext=ugext, m=m, distribution=distribution, F=F, tol=.Machine$double.eps^0.7)
  
  at = un$root

  # sample moderate voting matrix - ignore extremists
  sample_v = sample_votes(ugrid, at * ugrid, sample_u, N, u)
  
  mods = mod_dens(sample_u, sample_v, u, N)
  # density for sum of moderate votes
  ge = mods$g
  # largest absolute sum of votes in sample
  max_V = mods$max_V
  # rescale - this is a mixture component for the full density
  gt = function(v){return((1-(N-1)*p)*ge(v))}
  
  # solve for the extremist vote level at the discontinuity
  z <- uniroot.all(foc,lower=-m*max_V,upper=0, u=u, g=gt)
  z = sort(z)
  # should be extremist vote level
  v = z[1]
  
  # extremist voting grid   
  v_ext = bisection(foc, lower=-m*max_V, upper=v, u=ugext, g=gt)$mid
  
  # return the moderate slope, extremist utility grid, extremist voting grid, and overall voting grid
  return(list(a=at,ugext=ugext,v_ext=v_ext,vgrid=c(v_ext,at*ugrid[ugrid > u])))
}

### alternative initialization method. not currently in use.
# function to get the initial voting function using the 
# limiting estimate for the initial moderate slope
limslope = function(ubar, N, mu, ugrid, sample_u, u, m,F,distribution){
  
  # extremist grid
  ugext = ugrid[ugrid <= u]
  p = F(u)
  
  # moderate slope
  at = sqrt(abs(ubar) / (mu*(N-1)))
  
  # sample moderate voting matrix - ignore extremists
  sample_v = sample_votes(ugrid, at * ugrid, sample_u, N, u)
  
  mods = mod_dens(sample_u, sample_v, u, N)
  # density with only moderates
  ge = mods$g
  # CDF with only moderates
  Ge = mods$G
  # largest absolute sum of votes in sample
  max_V = mods$max_V
  
  gt = function(v){return((1-(N-1)*p)*ge(v))}
  
  # solve for the extremist vote level at the discontinuity
  z <- uniroot.all(foc,lower=-m*max_V,upper=0, u=u, g=gt)
  z = sort(z)
  # should be extremist vote level
  v = z[1]
  
  # extremist voting grid   
  v_ext = bisection(foc, lower=-m*max_V, upper=v, u=ugext, g=gt)$mid
  
  vgrid=c(v_ext,at*ugrid[ugrid > u])
  # smeared density
  gsmear = make_smear(ugrid,vgrid, sample_u, sample_v, u, F, distribution, u)
  # smeared CDF
  Gsmear = CDF_smear(ugrid,vgrid, sample_u, sample_v, u, F, distribution, u)
  
  # full density - mixture
  g = function(v){return((1 - ((N-1)*p))*ge(v) + (N-1)*p*gsmear(v))}
  G = function(v){return((1 - ((N-1)*p))*Ge(v) + (N-1)*p*Gsmear(v))}
  
  # new desired discontinuity
  u = exact_dc(g, G, m * ubar, 0, -m * max_V, 0)
  
  mods = mod_dens(sample_u, sample_v, u, N)
  # density for sum of moderate votes
  ge = mods$g
  # largest absolute sum of votes in sample
  max_V = mods$max_V
  # rescale - this is a mixture component for the full density
  gt = function(v){return((1-(N-1)*p)*ge(v))}
  
  # solve for the extremist vote level at the discontinuity
  z <- uniroot.all(foc,lower=-m*max_V,upper=0, u=u, g=gt)
  z = sort(z)
  # should be extremist vote level
  v = z[1]
  
  # extremist voting grid   
  v_ext = bisection(foc, lower=-m*max_V, upper=v, u=ugext, g=gt)$mid
  
  return(list(a=at,ugext=ugext,v_ext=v_ext,vgrid=c(v_ext,at*ugrid[ugrid > u])))
}

# function to get the difference in welfare between acting as a moderate and acting
# as an extremist.
# The discontinuity is the utility level that causes this difference to be 0.
welf_diff = function(u, g, G, v_low, v_up){
  # find all the roots for a given utility
  z <- uniroot.all(foc,lower=v_low,upper=v_up, u=u, g=g)
  z = sort(z)
  # only one root
  if (length(z) == 1){
    # dc doesn't have a moderate root
    # being an extremist is more attractive
    if(G(-z) > .5){
      return(-100)
    }
    # dc doesn't have an extremist root
    # being a moderate is more attractive
    return(100)
  }
  # return the difference between the welfare of the extremist root and the welfare
  # of the moderate root
  return(util.fn(z[length(z)], u, G) - util.fn(z[1], u, G))
}

# finds an exact value for the discontinuity by finding the utility value for which
# the individual is indifferent between acting as a moderate and as an 
# extremist
exact_dc = function(g, G, u_low, u_up, v_low, v_up){
  return(uniroot(welf_diff, c(u_low, u_up), g=g, G=G, v_low=v_low, v_up=v_up)$root)
}
