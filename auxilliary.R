### Auxilliary functions

# Utility function given v,u, and G
util.fn = function(v,u,G){
  return(u*(1-G(-v))-v^2)
}

# first order condition for utility function
foc = function(v, u, g){
  dens = g(-v)
  dens[is.na(dens)] = 10^-12
  return(v - (u * dens / 2))
}

# solves first order condition to get optimal vote level
u.root = function(x, lower, upper, g){
  vg1 = uniroot(foc, c(lower, upper), u=x, g=g)$root
  return(vg1)
}

# bisection method to solve for roots
# allows vectorized root solving
# see http://r.789695.n4.nabble.com/vectorized-uni-root-td4648920.html
bisection <- function(f, lower, upper, ..., numiter =100, tolerance = 
                        .Machine$double.eps^0.25 ) { 
  stopifnot(length(lower) == length(upper)) 
  
  flower <- f(lower, ...); fupper <- f(upper, ...) 
  for (n in 1:numiter) { 
    mid <- (lower+upper)/2 
    fmid <- f(mid, ...) 
    if (all(abs(fmid) < tolerance) & (n > 1)) break 
    samesign <- ((fmid<0)&(flower<0))|((fmid>=0)&(flower>=0)) 
    lower <- ifelse(samesign, mid, lower ) 
    flower <- ifelse(samesign, fmid, flower ) 
    upper <- ifelse(!samesign, mid, upper ) 
    fupper <- ifelse(!samesign, fmid, fupper ) 
  } 
  return(list( mid=mid, fmid=fmid, lower=lower, upper=upper, 
               flower=flower, fupper=fupper, n=n )) 
} 

# vectorized method to interpolate votes for the sample utility matrix
# assigns same vote as nearest utility value in the utility grid
sample_votes = function(ugrid, vgrid_0, sample_u, N, dc){
  a = ugrid[1]
  n = length(ugrid)
  b = ugrid[n]
  sl = (n - 1) / (b - a)
  # linear map to index in utility grid
  cst = 1 - (sl * a)
  samp_u = as.vector(sample_u)
  # apply linear map to values in sample_u to get the corresponding utility index
  # in ugrid
  inds = samp_u * sl + cst
  # find the closest utility value in the grid for each sample utility value
  inds = round(inds)
  # have to fix cases around dc
  if (!is.na(dc)){
    # find index for discontinuity
    dc_ind = dc * sl + cst
    # closest extremist in grid
    dc_indr = floor(dc_ind)
    # some extremists were misclassified
    if (dc_ind - dc_indr > .5){
      # rounded to moderate but should be extremist
      inds[(inds==(dc_indr+1)) & (samp_u <= dc)] = dc_indr
    }
    # some moderates were misclassified
    else
      # rounded to extremist but should be moderate
      inds[(inds==dc_indr) & (samp_u > dc)] = dc_indr + 1
  }
  inds[inds < 1] = 1
  inds[inds > n] = n
  return(matrix(vgrid_0[inds], ncol=N-1))
}