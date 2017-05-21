
# auxilliary function for dPLN pdf
R = function(x){
  return((1-pnorm(x,0,1))/dnorm(x,0,1))
}

initialize = function(distribution, para1 = NULL, para2 = NULL, 
                      N, range = NULL, K=1000, S_size){
  if(distribution == "Normal"){
    # Mean of normal distribution
    mu = para1
    # Standard deviation of normal distribution
    stdv = para2
    # Interval of utility values to be considered
    interval = c(0,0)
    # Since the density out here is almost zero.
    interval[1] = mu - 5.8*stdv
    interval[2] = mu + 5.8*stdv 
    
    # Sample from the utility distribution
    sample_uN = matrix(rnorm(S_size*(N), mean = mu, sd = stdv), ncol=N)
    
    # density of utility distribution
    f = function(u){
      return(dnorm(u,mu,stdv))
    }
    # CDF of utility distribution
    F = function(u){
      return(pnorm(u,mu,stdv))
    }
    
    # quantile function
    Qt = function(p){
      return(qnorm(p, para1, para2))
    }
    
    sampler = function(n, p=NA, lower.tail=TRUE){
      if (is.na(p))
        return(rnorm(n, mu, stdv))
      if (lower.tail){
        vals = runif(n, min=0, max=p)
        return(Qt(vals))
      }
      vals = runif(n, min=p, max=1)
      return(Qt(vals))
    }
    
    limiting = 0
    if (mu != 0){
      limiting = uniroot(function(q){return(q-N*pnorm(-N* mu*q, mean=mu, sd=stdv))}, 
              interval=c(0, 1))$root
    }
  }
  if(distribution == "Beta"){
    # parameters of Beta distribution
    a = para1
    b = para2
    
    # change location and scale of the density
    interval = range
    location = interval[1]
    scale = interval[2]-location
    
    # Sample from the utility distribution
    sample_uN = matrix(rpearsonI(S_size*(N),a = a,b=b,location = location, scale = scale), ncol=N)
    
    # mean for Pearson Type 1
    mu = (a / (a + b))*(scale) + location
    # variance
    var = ((a * b)*scale^2) / ((a + b)^2 * (a + b + 1))
    # standard deviation
    stdv = sqrt(var)
    
    # density for u ~ PearsonI(a,b,l,s)
    f = function(u)
    {
      return(dpearsonI(u,a,b,location,scale))
    }
    
    # CDF for u ~ PearsonI(a,b,l,s)
    F = function(u)
    {
      return(ppearsonI(u,a,b,location,scale))
    }
    
    Qt = function(p){
      return(bisection(function(x){ return(F(x) - p)}, 
                 lower=interval[1], upper=0, tolerance=.Machine$double.eps^0.75)$mid)
    }
    
    sampler = function(n, p=NA, lower.tail=TRUE){
      if (is.na(p))
        return(rpearsonI(n,a = a,b=b,location = location, scale = scale))
      if (lower.tail){
        vals = runif(n, min=0, max=p)
        return(Qt(vals))
      }
      vals = runif(n, min=p, max=1)
      return(Qt(vals))
    }
    limiting = -interval[1] / (mu * (N-1))
  }
  
  if(distribution == "Uniform"){
    
    # sample from utility distribution
    sample_uN = matrix(runif(S_size*(N), range[1],range[2]), ncol=N)
    
    # bounds of uniform distribution
    interval = range
    
    # mean and standard deviation of uniform distribution
    mu = .5 * (range[1] + range[2])
    stdv = (range[2] - range[1]) / (sqrt(12))
    
    # Uniform density
    f = function(u){
      if(u < range[1])
        return(0)
      if(u > range[2])
        return(0)
      return(dunif(u,range[1],range[2]))
    }
    # Uniform CDF
    F = function(u){
      if(u < range[1])
        return(0)
      if(u > range[2])
        return(1)
      return(punif(u,range[1],range[2]))
    }
    
    # quantile function
    Qt = function(p){
      return(qunif(p, range[1], range[2]))
    }
    
    sampler = function(n, p=NA, lower.tail=TRUE){
      if (is.na(p))
        return(runif(n, range[1],range[2]))
      if (lower.tail){
        vals = runif(n, min=0, max=p)
        return(Qt(vals))
      }
      vals = runif(n, min=p, max=1)
      return(Qt(vals))
    }
    limiting = -interval[1] / (mu * (N-1))
  }
  
  ### Simulation of 2008 Gay Marriage vote in California assuming uniform mixture components
  if(distribution == "Mar"){
    
    # mixture proportions
    alpha1 = .616
    alpha2 = .344
    alpha3 = .033
    alpha4 = .007
    
    # against gay marriage
    min_1 = -10000
    max_1 = 0
    
    # non-LGBT in favor of gay marriage
    min_2 = 0
    max_2 = 10000
    
    # LGBT, single
    min_3 = 5000
    max_3 = 35000
    
    # LGBT, couple
    min_4 = 20000
    max_4 = 180000
    
    interval = c(1.5 * min_1, 1.5 * max_4)
    
    # mean of utility distribution
    mu = alpha1 * (.5 * (min_1 + max_1)) + alpha2 * (.5 * (min_2 + max_2)) + alpha3 * (.5 * (min_3 + max_3)) + alpha4 * (.5 * (min_4 + max_4))
    
    # second moment of uniform distribution
    sec_mom = function(a, b){
      return((1/12) * (b - a)^2 + (.5 * (a + b))^2)
    }
    
    # mean squared
    Mu2 = mu^2
    # second moment of mixture
    EX2 = alpha1 * sec_mom(min_1, max_1) + alpha2 * sec_mom(min_2, max_2) + alpha3 * sec_mom(min_3, max_3) + alpha4 * sec_mom(min_4, max_4)
    # standard deviation of mixture
    stdv = (EX2 - Mu2)^(1/2)
    
    alphas = c(alpha1, alpha2, alpha3, alpha4)
    mins = c(min_1, min_2, min_3, min_4)
    maxs = c(max_1, max_2, max_3, max_4)
    
    # draw from multinomial to get subpopulation counts
    # in sample utility matrix
    components = rmultinom(1, S_size * N, alphas)
    
    # draw sample utility values
    samp1 = runif(components[1],min_1, max_1)
    samp2 = runif(components[2],min_2, max_2)
    samp3 = runif(components[3],min_3, max_3)
    samp4 = runif(components[4],min_4, max_4)
    
    full_mat = c(samp1, samp2, samp3, samp4)
    
    # jumble and store in matrix
    sample_uN = matrix(sample(full_mat), ncol=N)
    
    # CDF restricted to people against gay marriage
    # used when adjusting movement in the discontinuity point 
    # between periods
    F = function(x){
      return(alpha1*punif(x, min_1, max_1))
    }
    
    # pdf restricted to people against gay marriage
    # used when adjusting movement in the discontinuity point 
    # between periods
    f = function(x){
      return(alpha1*dunif(x, min_1, max_1))
    }
    
    # Quantile function
    Qt = function(p){
      return(bisection(function(x){ return(F(x) - p)}, 
              lower=interval[1], upper=0, tolerance=.Machine$double.eps^0.75)$mid)
    }
    
    # function to sample from the distribution
    sampler = function(n, p=NA, lower.tail=TRUE){
      if (is.na(p)){
        return(sample(full_mat, n, replace=TRUE))
      }
      
      # if there's a discontinuity then only sample from the extremist mass  
      if (lower.tail){
        vals = runif(n, min=0, max=p)
        return(Qt(vals))
      }
      
      return(sample(full_mat[full_mat > Qt(p)], n, replace=TRUE))
    }
    
    Limiting = -min_1 / (mu * (N-1))
  }
  
  ### Simulation of 2008 Gay Marriage vote in California assuming dPLN distributed mixture components
  if (distribution == "dPLN"){
    
    # dPLN parameters fit to US Income Distribution
    a = 3 # upper tail index
    b = 1.43 # lower tail index
    s = 0.45 # ~ sd of log-income
    
    # Mixture proportions
    
    # against
    alpha1 = .52
    # non-LGBT, in favor
    alpha2 = .44
    # LGBT, single
    alpha3 = .033
    # LGBT, couple
    alpha4 = .007
    
    alphas = c(alpha1, alpha2, alpha3, alpha4)
    
    # mean willingness to pay for non-LGBT
    mean1 = 5000
    # mean for LGBT, single
    mean2 = 20000
    # mean for LGBT, couple
    mean3 = 100000
    
    m1=log(mean1)-log(((a*b)/((a-1)*(b+1))))-(s^2)/2 # parameter for voters for prop 8
    m2=log(mean1)-log(((a*b)/((a-1)*(b+1))))-(s^2)/2 # parameter for non-LGBT against prop 8
    m3=log(mean2)-log(((a*b)/((a-1)*(b+1))))-(s^2)/2 # parameter for LGBT, single
    m4=log(mean3)-log(((a*b)/((a-1)*(b+1))))-(s^2)/2 # parameter for LGBT, couple  
    
    # vector of parameters
    mvec=rep(0,4)
    mvec[1] = m1
    mvec[2] = m2
    mvec[3] = m3
    mvec[4] = m4
    
    
    
    #pdf for double pareto log-normal of subpop n
    f_dpln = function(w,n){
      m=mvec[n]
      if(n>1){
        z=w
      } else {
        z=-w
      }
    
      z[z<=0] = 1
      t1 = (a*b)/((a+b)*z)
      t2 = dnorm((log(z)-m)/s,0,1)
      r1 = R(a*s - (log(z)-m)/s)
      r2 = R(b*s + (log(z)-m)/s)
      out = t1*t2*(r1+r2)
      out[z==1] = 0
      return(out)
      
    }
    
    #CDF for double pareto log-normal
    F_dpln = function(w,n){
      vv = w
      vv[vv<=0] = 1
      
      m=mvec[n]
      t1 = pnorm((log(vv)-m)/s,0,1)
      t2 = dnorm((log(vv)-m)/s,0,1)
      x1 = a*s - (log(vv)-m)/s
      x2 = b*s + (log(vv)-m)/s
      t3 = (b*R(x1)-a*R(x2))/(a+b)
      out = (t1-t2*t3)
      out[w<=0] = 0
      return(out)
    }
    
    ##### appropriate proposal gamma distributions for each subpopulation
    
    ### parameters of gamma distribution
    ### same shape for those in favor of proposition 8 and non-LGBT opposed 
    al2 = 1.05
    be2 = (al2-1)/2602.603
    M2 = 9.7
    
    ### parameters of gamma distribution
    ### LGBT single voters
    al3 = 1.025
    be3 = (al3-1)/10410.41
    M3 = 18
    
    ### parameters of gamma distribution
    ### LGBT couples
    al4 = 1.025
    be4 = (al4-1)/52052.05
    M4 = 18
    
    # vectors of parameters
    alps = c(0,al2,al3,al4)
    bets = c(0,be2,be3,be4)
    Ms = c(0,M2,M3,M4)
    
    ### vectorized rejection sampling algorithm for dPLN
    # bign is the number of observations to draw
    # n is the subpopulation
    # note the proposal distribution is less than the dPLN density at the far end of the 
    # fat tail. We sought to make the mass of that region sufficiently small to be 
    # negligible. 
    rejection = function(bign,n){
      # for voters against gay marriage draw from the same
      # distributions as moderates in favor, flip signs.
      if(n==1){return(-rejection(bign,2))}
      # valid values drawn so far
      curr_samp = c()
      # values still needed
      vals_left = bign
      # proposal distribution parameters
      al = alps[n]
      be = bets[n]
      M = Ms[n]
      
      # keep drawing until you get full sample of values 
      while(vals_left > 0){
        u = runif(vals_left, 0, 1)
        # draw from proposal distribution
        g = rgamma(vals_left, al, be)
        # get rescaled height under curve
        M_vals = M * dgamma(g, al, be)
        # get height under curve for dPLN
        f_vals = f_dpln(g, n)
        # keep if sampled point under the dPLN curve
        accept = (u < (f_vals / M_vals))
        # set all rejected values to 0
        accept_samp = accept * g
        # add accepted values to sample
        curr_samp = c(curr_samp, accept_samp[accept_samp != 0])
        vals_left = bign - length(curr_samp)
      }
      return(curr_samp)
    }
    
    # draw from multinomial to get subpopulation counts
    # in sample utility matrix
    components = rmultinom(1, S_size * N, alphas)
    
    if ("sample_matrix_dpln.csv" %in% list.files()){
      full_samp = as.vector(as.matrix(read.csv("sample_matrix_dpln.csv")))
      full_samp = full_samp[1:(S_size * N)]
    }
    else{
      # draw sample utility values
      samp1 = rejection(components[1],1) # against prop 8
      samp2 = rejection(components[2],2) # for prop 8 non-LGBT
      samp3 = rejection(components[3],3) # LGBT single
      samp4 = rejection(components[4],4) # LGBT couple
      full_samp = c(samp1, samp2, samp3, samp4)
    }
    
    # mean for mixture
    mu = (alpha2 - alpha1) * mean1 + alpha3 * mean2 + alpha4 * mean3
    
    # variance
    sigma2 = var(full_samp)
    
    # estimated standard deviation
    stdv = sqrt(sigma2)
    
    # CDF restricted to people against gay marriage
    # used when adjusting movement in the discontinuity point 
    # between periods
    F = function(x){
      return(alpha1*(1-F_dpln(-x,2)))
    }
    
    # pdf restricted to people against gay marriage
    # used when adjusting movement in the discontinuity point 
    # between periods
    f = function(x){
      return(alpha1*(f_dpln(-x,2)))
    }
    
    # quantile for double pareto log-normal
    Qt = function(p){
      tf = function(y){return(F(y)-p)}
      rt = bisection(tf,lower=-10^8, upper=0, tolerance=.Machine$double.eps^0.75)
      return(rt$mid)
    }
    
    # bounds for utility grid
    interval = c(1.05 * min(full_samp), 1.05 * max(full_samp))
    
    # merge together and shuffle up the utility values for each 
    # subpopulation. convert to a matrix.
    sample_uN = matrix(full_samp, ncol=N)
    
    # sample from distribution
    sampler = function(n, p=NA, lower.tail=TRUE){
      if (is.na(p)){
        return(sample(full_samp, n, replace=TRUE))
      }
      
      if (lower.tail){
        vals = runif(n, min=0, max=p)
        return(Qt(vals))
      }
      
      return(sample(full_samp[full_samp > Qt(p)], n, replace=TRUE))
    }
    
    # limiting constant
    k.k = (F(-1000000)/alpha1) / 1000000^-a
    
    # constant / sqrt(N) for limiting
    limiting = k.k^(1/(1 + a)) * ((1/mu) + (1/(mu*a)))^(a/(1+a)) / sqrt(N)
  }
  
  if (distribution == "uDPLN"){
    
    # dPLN parameters fit to US Income Distribution
    a = 3 # upper tail index
    b = 1.43 # lower tail index
    s = 0.45 # ~ sd of log-income
    
    # uniform distribution parameters
    range = c(-5000,0)
    
    # Mixture proportions
    
    # uniform
    alpha1 = .5
    # dPLN
    alpha2 = .5
    
    alphas = c(alpha1, alpha2)
    
    # mean for dPLN
    mean1 = (para1 - (alpha2 * -2500)) / alpha1
    
    m1=log(mean1)-log(((a*b)/((a-1)*(b+1))))-(s^2)/2
    
    #pdf for double pareto log-normal of subpop n
    f_dpln = function(w){
      m=m1
      z=w
      z[z<=0] = 1
      t1 = (a*b)/((a+b)*z)
      t2 = dnorm((log(z)-m)/s,0,1)
      r1 = R(a*s - (log(z)-m)/s)
      r2 = R(b*s + (log(z)-m)/s)
      out = t1*t2*(r1+r2)
      out[z==1] = 0
      return(out)
    }
    
    #CDF for double pareto log-normal
    F_dpln = function(w){
      vv = w
      vv[w<=0] = 1
      
      m=m1
      t1 = pnorm((log(vv)-m)/s,0,1)
      t2 = dnorm((log(vv)-m)/s,0,1)
      x1 = a*s - (log(vv)-m)/s
      x2 = b*s + (log(vv)-m)/s
      t3 = (b*R(x1)-a*R(x2))/(a+b)
      out = t1-t2*t3
      out[w<=0] = 0
      return(out)
    }
    
    
    al1 = 1.05
    be1 = (al1-1)/2602.603

    # find multiplier of Gamma distribution that makes it twice as tall
    # as DPLN at DPLN mode
    modevals = optimize(function(x){return(-f_dpln(x))}, c(1, mean1))
    M1 = -2 * modevals$objective / dgamma(modevals$minimum, al1, be1)
    
    rejection = function(bign){
      curr_samp = c()
      # values still needed
      vals_left = bign
      # proposal distribution parameters
      al = al1
      be = be1
      M = M1
      
      # keep drawing until you get full sample of values 
      while(vals_left > 0){
        u = runif(vals_left, 0, 1)
        # draw from proposal distribution
        g = rgamma(vals_left, al, be)
        # get rescaled height under curve
        M_vals = M * dgamma(g, al, be)
        # get height under curve for dPLN
        f_vals = f_dpln(g)
        # keep if sampled point under the dPLN curve
        accept = (u < (f_vals / M_vals))
        # set all rejected values to 0
        accept_samp = accept * g
        # add accepted values to sample
        curr_samp = c(curr_samp, accept_samp[accept_samp != 0])
        vals_left = bign - length(curr_samp)
      }
      return(curr_samp)
    }
    
    
    # draw from multinomial to get subpopulation counts
    # in sample utility matrix
    components = rmultinom(1, S_size * N, alphas)
    
    samp1 = rejection(components[1])
    
    # estimated second moment of utility distribution for non-LGBT
    secmom1 = var(samp1) + mean1^2
    
    # sample from utility distribution
    samp2 = runif(components[2], range[1],range[2])
    
    # mean and standard deviation of uniform distribution
    mean2 = .5 * (range[1] + range[2])
    secmom2 = mean2^2 + (range[2] - range[1])^2 / 12 
    
    # mean for mixture
    mu = alpha1 * mean1 + alpha2 * mean2
    
    # estimated second moment for mixture
    secmom = alpha1 * secmom1 + alpha2 * secmom2
    
    stdv = sqrt(secmom - mu^2)
    
    # Uniform density
    f = function(u){
      if(u < range[1])
        return(0)
      if(u > range[2])
        return(0)
      return(alpha1*dunif(u,range[1],range[2]))
    }
    # Uniform CDF
    F = function(u){
      if(u < range[1])
        return(0)
      if(u > range[2])
        return(1)
      return(alpha1*punif(u,range[1],range[2]))
    }
    
    # quantile function
    Qt = function(p){
      return(qunif(p / alpha1, range[1], range[2]))
    }
    
    min_1 = min(samp1)
    max_1 = max(samp1)
    
    # bounds for utility grid
    interval = c(range[1], 1.05 * max(max_1, range[2]))
    
    full_samp = c(samp1, samp2)
    
    # merge together and shuffle up the utility values for each 
    # subpopulation. convert to a matrix.
    sample_uN = matrix(sample(full_samp), ncol=N)
    
    sampler = function(n, p=NA, lower.tail=TRUE){
      if (is.na(p)){
        return(sample(full_samp, n, replace=TRUE))
      }
      
      if (lower.tail){
        vals = runif(n, min=0, max=p)
        return(Qt(vals))
      }
      
      return(sample(full_samp[full_samp > Qt(p)], n, replace=TRUE))
    }
    limiting = -interval[1] / (mu * (N-1))
  }
  
  if (distribution == "DPLNu"){
    
    a = para2
    #a = 3 # upper tail index
    b = 1.43 # lower tail index
    s = 0.45 # ~ sd of log-income
    
    # uniform distribution parameters
    range = c(0,10000)
    
    # Mixture proportions
    
    # uniform
    alpha1 = .5
    # dPLN
    alpha2 = .5
    
    alphas = c(alpha1, alpha2)
    
    # mean for dPLN
    mean1 = (para1 - (alpha2 * 5000)) / alpha1
    mmean1 = -mean1
    
    m1=log(mmean1)-log(((a*b)/((a-1)*(b+1))))-(s^2)/2
    
    #pdf for double pareto log-normal of subpop n
    f_dpln = function(w){
      m=m1
      z=-w
      z[z<=0] = 1
      t1 = (a*b)/((a+b)*z)
      t2 = dnorm((log(z)-m)/s,0,1)
      r1 = R(a*s - (log(z)-m)/s)
      r2 = R(b*s + (log(z)-m)/s)
      out = t1*t2*(r1+r2)
      out[z==1] = 0
      return(out)
    }
    
    #CDF for double pareto log-normal
    F_dpln = function(w){
      vv = -w
      vv[vv<=0] = 1
      
      m=m1
      t1 = pnorm((log(vv)-m)/s,0,1)
      t2 = dnorm((log(vv)-m)/s,0,1)
      x1 = a*s - (log(vv)-m)/s
      x2 = b*s + (log(vv)-m)/s
      t3 = (b*R(x1)-a*R(x2))/(a+b)
      out = 1 - (t1-t2*t3)
      out[vv<=0] = 0
      return(out)
    }
    
    modevals = optimize(function(x){return(-f_dpln(x))}, c(mean1, -1))
    
    al1 = 1.05
    be1 = (al1-1)/2602.603
    
    # find multiplier of Gamma distribution that makes it twice as tall
    # as DPLN at DPLN mode
    M1 = -2 * modevals$objective / dgamma(-modevals$minimum, al1, be1)
    
    rejection = function(bign){
      # for voters against gay marriage draw from the same
      # distributions as moderates in favor, flip signs.
      # valid values drawn so far
      curr_samp = c()
      # values still needed
      vals_left = bign
      # proposal distribution parameters
      al = al1
      be = be1
      M = M1
      
      # keep drawing until you get full sample of values 
      while(vals_left > 0){
        u = runif(vals_left, 0, 1)
        # draw from proposal distribution
        g = rgamma(vals_left, al, be)
        # get rescaled height under curve
        M_vals = M * dgamma(g, al, be)
        # get height under curve for dPLN
        f_vals = f_dpln(-g)
        # keep if sampled point under the dPLN curve
        accept = (u < (f_vals / M_vals))
        # set all rejected values to 0
        accept_samp = accept * g
        # add accepted values to sample
        curr_samp = c(curr_samp, accept_samp[accept_samp != 0])
        vals_left = bign - length(curr_samp)
      }
      return(-curr_samp)
    }
    
    
    # draw from multinomial to get subpopulation counts
    # in sample utility matrix
    components = rmultinom(1, S_size * N, alphas)
    
    samp1 = rejection(components[1])
    
    # estimated second moment of utility distribution for non-LGBT
    secmom1 = var(samp1) + mean1^2
    
    # sample from utility distribution
    samp2 = runif(components[2], range[1],range[2])
    
    # mean and standard deviation of uniform distribution
    mean2 = .5 * (range[1] + range[2])
    secmom2 = mean2^2 + (range[2] - range[1])^2 / 12 
    
    # mean for mixture
    mu = alpha1 * mean1 + alpha2 * mean2
    
    # estimated second moment for mixture
    secmom = alpha1 * secmom1 + alpha2 * secmom2
    
    stdv = sqrt(secmom - mu^2)
    
    # CDF restricted to people against gay marriage
    # used when adjusting movement in the discontinuity point 
    # between periods
    F = function(x){
      return(alpha1*F_dpln(x))
    }
    
    # pdf restricted to people against gay marriage
    # used when adjusting movement in the discontinuity point 
    # between periods
    f = function(x){
      return(alpha1*(f_dpln(x)))
    }
    
    # quantile for double pareto log-normal
    Qt = function(p){
      tf = function(y){return(F(y)-p)}
      rt = bisection(tf,lower=-10^8, upper=0, tolerance=.Machine$double.eps^0.75)
      return(rt$mid)
    }
    
    min_1 = min(samp1)
    max_1 = max(samp1)
    
    # bounds for utility grid
    interval = c(1.05*min_1, 1.05*range[2])
    
    full_samp = c(samp1, samp2)
    
    # merge together and shuffle up the utility values for each 
    # subpopulation. convert to a matrix.
    sample_uN = matrix(sample(full_samp), ncol=N)
    
    
    sampler = function(n, p=NA, lower.tail=TRUE){
      if (is.na(p)){
        return(sample(full_samp, n, replace=TRUE))
      }
      
      if (lower.tail){
        vals = runif(n, min=0, max=p)
        return(Qt(vals))
      }
      
      return(sample(full_samp[full_samp > Qt(p)], n, replace=TRUE))
    }
    
    k.k = (F(-1000000)/alpha1) / 1000000^-a

    # constant / sqrt(N) for limiting
    limiting = k.k^(1/(1 + a)) * ((1/mu) + (1/(mu*a)))^(a/(1+a)) / sqrt(N)
    
  }
  
  if (distribution == "Laplace"){
    
    m = para1
    s = para2
    # Mean of normal distribution
    mu = para1
    # Standard deviation of normal distribution
    stdv = sqrt(2) * para2
    # Interval of utility values to be considered
    interval = c(0,0)
    # Since the density 3 standard deviations away is almost zero.
    interval[1] = mu - 5.8*stdv
    interval[2] = mu + 5.8*stdv 
    
    # Sample from the utility distribution
    sample_uN = matrix(rlaplace(S_size*(N), m=m, s=s), ncol=N)
    
    # density of utility distribution
    f = function(u){
      return(dlaplace(u,m,s))
    }
    # CDF of utility distribution
    F = function(u){
      return(plaplace(u,m,s))
    }
    
    # quantile function
    Qt = function(p){
      return(qlaplace(p, m, s))
    }
    
    sampler = function(n, p=NA, lower.tail=TRUE){
      if (is.na(p))
        return(rlaplace(n, m, s))
      if (lower.tail){
        vals = runif(n, min=0, max=p)
        return(Qt(vals))
      }
      vals = runif(n, min=p, max=1)
      return(Qt(vals))
    }
    
    # don't know what the limiting inefficiency is
    limiting = 0
  }
  
  if (distribution == "DPLN2"){
    
    a1 = para1
    a2 = para2
    as = c(a1,a2)
    b = 1.43 # lower tail index
    s = 0.45 # ~ sd of log-income
    
    # Mixture proportions
    
    # uniform
    alpha1 = .5
    # dPLN
    alpha2 = .5
    
    alphas = c(alpha1, alpha2)
    
    # mean for dPLN
    mean1 = -3000
    mmean1 = -mean1
    mean2 = 5000
    
    m1=log(mmean1)-log(((a1*b)/((a1-1)*(b+1))))-(s^2)/2
    m2=log(mean2)-log(((a2*b)/((a2-1)*(b+1))))-(s^2)/2
    
    mvec = c(m1, m2)
    
    #pdf for double pareto log-normal of subpop n
    f_dpln = function(w, n){
      m=mvec[n]
      a=as[n]
      if (n > 1){
        z = w
      } 
      else {
        z = -w
      }
      
      z[z<=0] = 1
      t1 = (a*b)/((a+b)*z)
      t2 = dnorm((log(z)-m)/s,0,1)
      r1 = R(a*s - (log(z)-m)/s)
      r2 = R(b*s + (log(z)-m)/s)
      out = t1*t2*(r1+r2)
      out[z==1] = 0
      return(out)
    }
    
    #CDF for double pareto log-normal
    F_dpln = function(w, n){
      if (n == 1){
        vv = -w
      }
      else{
        vv = w
      }
      
      vv[vv<=0] = 1
      a = as[n]
      m = mvec[n]
      
      t1 = pnorm((log(vv)-m)/s,0,1)
      t2 = dnorm((log(vv)-m)/s,0,1)
      x1 = a*s - (log(vv)-m)/s
      x2 = b*s + (log(vv)-m)/s
      t3 = (b*R(x1)-a*R(x2))/(a+b)
      out = (t1-t2*t3)
      out[vv==1] = 0
      if (n == 1)
        return(1 - out)
      return(out)
    }
    
    ## find mode of density function. used to find appropriate 
    ## multiplier for proposal gamma distribution.
    
    modevals1 = optimize(function(x){return(-f_dpln(x, 1))}, c(mean1, -1))
    
    al1 = 1.05
    be1 = (al1-1)/2602.603
    
    # M1 * gamma / dpln = 2 at mode
    M1 = -2 * modevals1$objective / dgamma(-modevals1$minimum, al1, be1)
    
    
    modevals2 = optimize(function(x){return(-f_dpln(x, 2))}, c(1, mean2))
    
    al1 = 1.05
    be1 = (al1-1)/2602.603
    
    M2 = -2 * modevals2$objective / dgamma(modevals2$minimum, al1, be1)
    
    Ms = c(M1, M2)
    
    rejection = function(bign, n){
      # valid values drawn so far
      curr_samp = c()
      # values still needed
      vals_left = bign
      # proposal distribution parameters
      al = al1
      be = be1
      M = Ms[n]
      
      # keep drawing until you get full sample of values 
      while(vals_left > 0){
        u = runif(vals_left, 0, 1)
        # draw from proposal distribution
        g = rgamma(vals_left, al, be)
        # get rescaled height under curve
        M_vals = M * dgamma(g, al, be)
        
        if(n == 1){
          g=-g
        }
        
        # get height under curve for dPLN
        f_vals = f_dpln(g, n)
        # keep if sampled point under the dPLN curve
        accept = (u < (f_vals / M_vals))
        # set all rejected values to 0
        accept_samp = accept * g
        accept_samp[is.na(accept_samp)] = 0
        # add accepted values to sample
        curr_samp = c(curr_samp, accept_samp[accept_samp != 0])
        vals_left = bign - length(curr_samp)
      }
      return(curr_samp)
    }
    
    
    # draw from multinomial to get subpopulation counts
    # in sample utility matrix
    components = rmultinom(1, S_size * N, alphas)
    
    samp1 = rejection(components[1], 1)
    samp2 = rejection(components[2], 2)
    full_samp = c(samp1, samp2)
    
    # mean for mixture
    mu = alpha1 * mean1 + alpha2 * mean2
    
    var_dpln = function(a, b, s, v){
      p1 = (a * b * exp(2*v + s^2)) / ((a - 1)^2 * (b + 1)^2)
      p2 = (((a - 1)^2 * (b + 1)^2) / ((a - 2) * (b + 2))) * exp(t^2)
      p3 = a * b
      return(p1 * (p2 - p3))
    }
    
    sigma2 = var(full_samp)
    stdv = sqrt(sigma2)
    
    # CDF restricted to people against gay marriage
    # used when adjusting movement in the discontinuity point 
    # between periods
    F = function(x){
      return(alpha1*F_dpln(x, 1))
    }
    
    # pdf restricted to people against gay marriage
    # used when adjusting movement in the discontinuity point 
    # between periods
    f = function(x){
      return(alpha1*(f_dpln(x, 1)))
    }
    
    # quantile for double pareto log-normal
    Qt = function(p){
      tf = function(y){return(F(y)-p)}
      rt = bisection(tf,lower=-10^8, upper=0, tolerance=.Machine$double.eps^0.75)
      return(rt$mid)
    }
    
    min_1 = min(samp1)
    max_1 = max(samp2)
    
    # bounds for utility grid
    interval = c(1.05*min_1, 1.05*max_1)
    
    # merge together and shuffle up the utility values for each 
    # subpopulation. convert to a matrix.
    sample_uN = matrix(sample(full_samp), ncol=N)
    
    
    sampler = function(n, p=NA, lower.tail=TRUE){
      if (is.na(p)){
        return(sample(full_samp, n, replace=TRUE))
      }
      
      if (lower.tail){
        vals = runif(n, min=0, max=p)
        return(Qt(vals))
      }
      
      return(sample(full_samp[full_samp > Qt(p)], n, replace=TRUE))
    }
    
    k.k = (F(-1000000)/alpha1) / 1000000^-a1
    
    # constant / sqrt(N) for limiting
    limiting = k.k^(1/(1 + a1)) * ((1/mu) + (1/(mu*a1)))^(a1/(1+a1)) / sqrt(N)
    
  }
  
  
  
  return(list(mu=mu,stdv=stdv,interval=interval,sample_uN=sample_uN,f=f,F=F, Qt=Qt,
              sampler=sampler, limiting=limiting))
}


