
library(bayesm)
library(parallel)
library(abind)

rmnorm <- function (n, mu, root)  #modified from 'mvtnorm' library function rmvnorm, returns n multivariate normal rvs
{
  z = matrix(rnorm(n * length(mu)), nrow = n) %*% root
  z = sweep(z, 2, mu, "+")
  t(z)
}

# Proposed algorithm
hierBayesParallel <- function(Data, Prior, Mcmc, Options)
{
  Data = partition_data(Data, Options$S1)
  if(is.null(Mcmc$drawcomp)) {Mcmc$drawcomp=FALSE}
  if (Options$Sample != 1)
  {
    DataSample = vector(mode="list", length=Options$S1)
    for (i in 1:Options$S1)
    {
      Ns = Data[[i]]$Ns
      DataSample[[i]] = list(p = Data[[i]]$p, 
                             lgtdata = sample(Data[[i]]$lgtdata, size = max(1, round(Options$Sample * Ns))), 
                             Ns = Ns, N = Data[[i]]$N)
    }
    out_S = mclapply(DataSample, FUN=rhierMnlRwMixtureParallel, Prior = Prior, Mcmc = Mcmc,
                     mc.cores=Options$cpus1, mc.set.seed = FALSE)
    rm(DataSample)
  } else {
    out_S = mclapply(Data, FUN=rhierMnlRwMixtureParallel, Prior = Prior, Mcmc = Mcmc,
                     mc.cores=Options$cpus1, mc.set.seed = FALSE)
  }
  if (Options$S1 != Options$S2) {
    Data = partition_data(Data, Options$S2) 
  }
  Draws = combineDraws(out_S) 
  out2_S = mclapply(Data, FUN=IndepMetropParallel, Draws, Mcmc, Options,
                    mc.cores=Options$cpus2, mc.set.seed = FALSE)
  R = Mcmc$R2
  betadraw = do.call("abind", list(lapply(out2_S, function(x) x$betadraw[,,1:R]), along=1))
  attr(betadraw,"dimnames")<-NULL
  if (Mcmc$drawcomp == TRUE)
  {
    return(list(betadraw=betadraw, stats=stats, 
                compdraw=lapply(out_S, function(x) { return(x$compdraw) })))
  } else {
    return(list(betadraw=betadraw))
  }
}

# Partition data into S shards 
partition_data <- function(Data, S)
{
  if (is.list(Data[[1]]) == TRUE) 
  {
    S_in = length(Data)
  } else {
    S_in = 1
  }
  if (S_in == S)
  {
    return(Data)
  }
  if (is.list(Data[[1]]) == TRUE) 
  {
    S_in = length(Data)
    lgtdata = NULL
    for (i in 1:S_in)
    {
      lgtdata = c(lgtdata, Data[[i]]$lgtdata)
    }
    p = Data[[1]]$p
    N = Data[[1]]$N
    Nunits = Data[[1]]$Nunits
    Data = list(p = p, lgtdata=lgtdata, Ns=N, N=N, Nunits=Nunits, Nsunits=Nunits)
  }
  if (S == 1)
  {
    return(Data)
  }
  p = Data$p
  lgtdata = Data$lgtdata
  N = Data$N
  Nunits = Data$Nunits
  Data_S = NULL
  size = rep(floor(N/S), S)
  if (N != sum(size))
  {
    for (i in 1:(N-sum(size)))
    {
      size[i] = size[i] + 1
    }
  }
  index1 = 1
  Ncum = 0
  for (i in 1:S)
  {
    index2 = index1 + size[i] - 1
    if (Ncum < Nunits)
    {
      Nsunits = min(Nunits-Ncum, size[i])
      Ncum = Ncum + Nsunits
    } else {
      Nsunits = 0
    }
    Data_S[[i]] = list(p=p, lgtdata=lgtdata[index1:index2], Ns=size[i], N=N, Nunits=Nunits, Nsunits=Nsunits)
    index1 = index2 + 1
  }
  return(Data_S)
}

# Combine draws from Stage one of proposed algorithm
combineDraws <- function(draws_S)
{
  R = draws_S[[1]]$R
  S1 = length(draws_S)
  draws = do.call("rbind", lapply(draws_S, function(x) x$betadraw))
  ndraws = dim(draws)[1]
  draws = draws[sample(1:ndraws, R, replace = FALSE),]
  return(draws)
}

# Simulate data function
simulate <- function(N, p, k, Ti_min, Ti_max, mu, sigma)
{
  beta = matrix(0, nrow = N, ncol = k)
  for (i in 1:N)
  {
    beta[i,] = rmnorm(1, mu, chol(sigma))
  }
  Ti = round(runif(N, Ti_min, Ti_max))
  lgtdata = vector(mode = "list", length = N)
  for (i in 1:N)
  {
    e = matrix(-log(-log(runif(Ti[i] * p))), nrow=Ti[i], ncol=p)  #matrix of type 1 extreme value errors
    x = matrix(rlnorm(Ti[i]*p*(k-p+1)), nrow=Ti[i], ncol=p*(k-p+1))   #prices for p alternatives
    u = matrix(nrow=Ti[i], ncol=p)
    for (j in 1:p) {
      u[,j] = x[,seq.int(from=j,by=p,length.out=k-p+1)] * beta[i,p:k] + e[,j]   #utilities for p alternatives
    }
    y = max.col(u)  #vector of alternatives chosen
    X = createX(p=p, na=k-p+1, nd=NULL, Xa=x, Xd=NULL, base=1)
    lgtdata[[i]] = list(X = X, y = y)
  }
  return(list(p = p, lgtdata = lgtdata, beta = beta))
}

# Simulate data
simulate_data <- function(N, p, k, Ti_min, Ti_max, mu, sigma)
{
  lgtdata = vector(mode = "list", length = N)
  out = simulate(N, p, k, Ti_min, Ti_max, mu, sigma)
  return(list(p = p, lgtdata = out$lgtdata, beta = out$beta))
}


# Simulate data
RNGkind("L'Ecuyer-CMRG")
set.seed(1)
N = 10000 #number of units (panelists)
Nunits = 1000  #number of units to return
p = 4     #number of choice alternatives
k = 4     #number of parameters to estimate
Ti_min = 5   #min number of observations per unit
Ti_max = 5   #min number of observations per unit
mu = c(c(1, 2, 3), rnorm(k-p+1))
sigma = diag(k)

out_simulate = simulate_data(N, p, k, Ti_min, Ti_max, mu, sigma)
beta_true = out_simulate$beta 
Data = list(p=out_simulate$p, lgtdata=out_simulate$lgtdata, Ns=N, N=N, Nsunits=Nunits, Nunits=Nunits)


# Run proposed algorithm
R1 = 100    # Stage one iterations
Bpercent = 0.2    # Percent burn-in
B1 = floor(Bpercent * R1)    #Stage one burn-in iterations
R2 = R1    # Stage two iterations
B2 = B1    #Stage two burn-in
M = 4     #max number of available processors
S1 = M    # Stage one shards
S2 = M    # Stage two shards
cpus1 = S1
cpus2 = S2
K1 = 1    # Keep every K1th draw in stage one
K2 = 1    # Keep every K2th draw in stage two
Sample = 1   # Stage one subsampling rate
Data = partition_data(Data, 1)
Prior = list(ncomp = 1)
Mcmc = list(R = R1, R2 = R2, B = B1, B2 = B2, keep = K1, keep2 = K2, 
            s = NULL, w = NULL)
Options = list(cpus1=S1, cpus2=S2, S1=S1, S2=S2, M=M, Sample=Sample, Nunits=Nunits,
               DiscardBetas = FALSE, vectorize=TRUE, P=1)
out <- hierBayesParallel(Data = Data, Prior = Prior, Mcmc = Mcmc, Options = Options)
str(out)


