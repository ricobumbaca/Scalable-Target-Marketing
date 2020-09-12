rhierMnlRwMixtureParallel=function(Data,Prior,Mcmc){
#
# revision history:
#   changed 12/17/04 by rossi to fix bug in drawdelta when there is zero/one unit
#   in a mixture component
#   added loglike output, changed to reflect new argument order in llmnl, mnlHess 9/05
#   changed weighting scheme to (1-w)logl_i + w*Lbar (normalized) 12/05
#   3/07 added classes
#   9/08 changed Dirichlet a check
#   W. Taylor 4/15 - added nprint option to MCMC argument  
#   Rico Bumbaca 10/10/2019
#
# purpose: run hierarchical mnl logit model with mixture of normals 
#   using RW and cov(RW inc) = (hess_i + Vbeta^-1)^-1
#   uses normal approximation to pooled likelihood
#
# Arguments:
#   Data contains a list of (p,lgtdata, and possibly Z)
#      p is number of choice alternatives
#      lgtdata is a list of lists (one list per unit)
#          lgtdata[[i]]=list(y,X)
#             y is a vector indicating alternative chosen
#               integers 1:p indicate alternative
#             X is a length(y)*p x nvar matrix of values of
#               X vars including intercepts
#             Z is an length(lgtdata) x nz matrix of values of variables
#               note: Z should NOT contain an intercept
#   Prior contains a list of (deltabar,Ad,mubar,Amu,nu,V,ncomp) 
#      ncomp is the number of components in normal mixture
#           if elements of Prior (other than ncomp) do not exist, defaults are used
#   Mcmc contains a list of (s,c,R,keep,nprint)
#
# Output:  as list containing
#   Deltadraw R/keep  x nz*nvar matrix of draws of Delta, first row is initial value
#   betadraw is nlgt x nvar x R/keep array of draws of betas
#   probdraw is R/keep x ncomp matrix of draws of probs of mixture components
#   compdraw is a list of list of lists (length R/keep)
#      compdraw[[rep]] is the repth draw of components for mixtures
#   loglike  log-likelikelhood at each kept draw
#
# Priors:
#    beta_i = D %*% z[i,] + u_i
#       u_i ~ N(mu_ind[i],Sigma_ind[i])
#       ind[i] ~multinomial(p)
#       p ~ dirichlet (a)
#       D is a k x nz array
#          delta= vec(D) ~ N(deltabar,A_d^-1)
#    mu_j ~ N(mubar,A_mu^-1(x)Sigma_j)
#    Sigma_j ~ IW(nu,V^-1)
#    ncomp is number of components
#
# MCMC parameters
#   s is the scaling parameter for the RW inc covariance matrix; s^2 Var is inc cov
#      matrix
#   w is parameter for weighting function in fractional likelihood
#      w is the weight on the normalized pooled likelihood 
#   R is number of draws
#   keep is thinning parameter, keep every keepth draw
#   nprint - print estimated time remaining on every nprint'th draw
#
#  check arguments
#
if(missing(Data)) {pandterm("Requires Data argument -- list of p,lgtdata, and (possibly) Z")}
  if(is.null(Data$p)) {pandterm("Requires Data element p (# chce alternatives)") }
  p=Data$p
  if(is.null(Data$lgtdata)) {pandterm("Requires Data element lgtdata (list of data for each unit)")}
  nlgt=length(Data$lgtdata)
  nvar=ncol(Data$lgtdata[[1]]$X)
  drawdelta=TRUE
if(is.null(Data$Z)) { cat("Z not specified",fill=TRUE); #fsh() ; 
  drawdelta=FALSE}
  else {if (nrow(Data$Z) != nlgt) {pandterm(paste("Nrow(Z) ",nrow(Z),"ne number logits ",nlgt))}
      else {Z=Data$Z}}
  if(drawdelta) {
     nz=ncol(Z)
     colmeans=apply(Z,2,mean)
     if(sum(colmeans) > .00001) 
       {pandterm(paste("Z does not appear to be de-meaned: colmeans= ",colmeans))}
  }
#
# check lgtdata for validity
#
  #ypooled=NULL
  #Xpooled=NULL
  ny = 0
  nX = 0
  for (i in 1:nlgt)
  {
    ny = ny + length(Data$lgtdata[[i]]$y)
    nX = nX + nrow(Data$lgtdata[[i]]$X)
  }
  indy = 1
  indX = 1
  ypooled = vector(length=ny)
  Xpooled = matrix(0, nrow=nX, ncol=nvar)

if(!is.null(Data$lgtdata[[1]]$X)) {oldncol=ncol(Data$lgtdata[[1]]$X)}
for (i in 1:nlgt) 
{
    if(is.null(Data$lgtdata[[i]]$y)) {pandterm(paste("Requires element y of lgtdata[[",i,"]]"))}
    if(is.null(Data$lgtdata[[i]]$X)) {pandterm(paste("Requires element X of lgtdata[[",i,"]]"))}
  #    ypooled=c(ypooled,lgtdata[[i]]$y)
  ypooled[indy:(indy+length(Data$lgtdata[[i]]$y)-1)] = Data$lgtdata[[i]]$y
  indy = indy + length(Data$lgtdata[[i]]$y)
    nrowX=nrow(Data$lgtdata[[i]]$X)
    if((nrowX/p) !=length(Data$lgtdata[[i]]$y)) {pandterm(paste("nrow(X) ne p*length(yi); exception at unit",i))}
    newncol=ncol(Data$lgtdata[[i]]$X)
    if(newncol != oldncol) {pandterm(paste("All X elements must have same # of cols; exception at unit",i))}
    #    Xpooled=rbind(Xpooled,lgtdata[[i]]$X)
    Xpooled[indX:(indX+nrow(Data$lgtdata[[i]]$X)-1),] = Data$lgtdata[[i]]$X
    indX = indX + nrow(Data$lgtdata[[i]]$X)
    oldncol=newncol
}

nvar=ncol(Xpooled)
levely=as.numeric(levels(as.factor(ypooled)))
if(length(levely) != p) {pandterm(paste("y takes on ",length(levely)," values -- must be = p"))}
bady=FALSE
for (i in 1:p )
{
    if(levely[i] != i) bady=TRUE
}
cat("Table of Y values pooled over all units",fill=TRUE)
print(table(ypooled))
if (bady) 
  {pandterm("Invalid Y")}
#
# check on prior
#
if(missing(Prior)) 
{pandterm("Requires Prior list argument (at least ncomp)")} 
if(is.null(Prior$ncomp)) {pandterm("Requires Prior element ncomp (num of mixture components)")} else {ncomp=Prior$ncomp}
if(is.null(Prior$mubar)) {mubar=matrix(rep(0,nvar),nrow=1)} else { mubar=matrix(Prior$mubar,nrow=1)}
  if(ncol(mubar) != nvar) {pandterm(paste("mubar must have ncomp cols, ncol(mubar)= ",ncol(mubar)))}
if(is.null(Prior$Amu)) {Amu=matrix(BayesmConstant.A,ncol=1)} else {Amu=matrix(Prior$Amu,ncol=1)}
  if(ncol(Amu) != 1 | nrow(Amu) != 1) {pandterm("Am must be a 1 x 1 array")}
if(is.null(Prior$nu)) {nu=nvar+BayesmConstant.nuInc}  else {nu=Prior$nu}
  if(nu < 1) {pandterm("invalid nu value")}
if(is.null(Prior$V)) {V=nu*diag(nvar)} else {V=Prior$V}
  if(sum(dim(V)==c(nvar,nvar)) !=2) pandterm("Invalid V in prior")
if(is.null(Prior$Ad) & drawdelta) {Ad=BayesmConstant.A*diag(nvar*nz)} else {Ad=Prior$Ad}
if(drawdelta) {if(ncol(Ad) != nvar*nz | nrow(Ad) != nvar*nz) {pandterm("Ad must be nvar*nz x nvar*nz")}}
if(is.null(Prior$deltabar)& drawdelta) {deltabar=rep(0,nz*nvar)} else {deltabar=Prior$deltabar}
  if(drawdelta) {if(length(deltabar) != nz*nvar) {pandterm("deltabar must be of length nvar*nz")}}
if(is.null(Prior$a)) { a=rep(BayesmConstant.a,ncomp)} else {a=Prior$a}
if(length(a) != ncomp) {pandterm("Requires dim(a)= ncomp (no of components)")}
bada=FALSE
   for(i in 1:ncomp) { if(a[i] < 0) bada=TRUE}
  if(bada) pandterm("invalid values in a vector")
#
# check on Mcmc
#
if(missing(Mcmc)) 
  {pandterm("Requires Mcmc list argument")}
else 
   { 
    if(is.null(Mcmc$s)) {s=BayesmConstant.RRScaling/sqrt(nvar)} else {s=Mcmc$s}
    if(is.null(Mcmc$w)) {w=BayesmConstant.w}  else {w=Mcmc$w}
    if(is.null(Mcmc$keep)) {keep=BayesmConstant.keep} else {keep=Mcmc$keep}
    if(is.null(Mcmc$R)) {pandterm("Requires R argument in Mcmc list")} else {
      R=Mcmc$R
      }
    if(is.null(Mcmc$nprint)) {nprint=BayesmConstant.nprint} else {nprint=Mcmc$nprint}
      if(nprint<0) {pandterm('nprint must be an integer greater than or equal to 0')}
     if(is.null(Mcmc$drawcomp)) {Mcmc$drawcomp=FALSE}
   }

#
# print out problem
#
cat(" ",fill=TRUE)
cat("Starting MCMC Inference for Hierarchical Logit:",fill=TRUE)
cat("   Normal Mixture with",ncomp,"components for first stage prior",fill=TRUE)
cat(paste("  ",p," alternatives; ",nvar," variables in X"),fill=TRUE)
cat(paste("   for ",nlgt," cross-sectional units"),fill=TRUE)
cat(" ",fill=TRUE)
cat("Prior Parms: ",fill=TRUE)
cat("nu =",nu,fill=TRUE)
cat("V ",fill=TRUE)
print(V)
cat("mubar ",fill=TRUE)
print(mubar)
cat("Amu ", fill=TRUE)
print(Amu)
cat("a ",fill=TRUE)
print(a)
if(drawdelta) 
{
   cat("deltabar",fill=TRUE)
   print(deltabar)
   cat("Ad",fill=TRUE)
   print(Ad)
}
cat(" ",fill=TRUE)
cat("MCMC Parms: ",fill=TRUE)
cat(paste("s=",round(s,3)," w= ",w," R= ",R," keep= ",keep," nprint= ",nprint),fill=TRUE)
cat("",fill=TRUE)

oldbetas = matrix(double(nlgt * nvar), ncol = nvar)

#--------------------------------------------------------------------------------------------------
#
#  create functions needed
#
llmnlFract=
function(beta,y,X,betapooled,rootH,w,wgt){
z=as.vector(rootH%*%(beta-betapooled))
return((1-w)*llmnl(beta,y,X)+w*wgt*(-.5*(z%*%z)))
}
#-------------------------------------------------------------------------------------------------------
#
# intialize compute quantities for Metropolis
#
cat("initializing Metropolis candidate densities for ",nlgt," units ...",fill=TRUE)
#fsh()
#
#  now go thru and computed fraction likelihood estimates and hessians
#
#       Lbar=log(pooled likelihood^(n_i/N))
#
#       fraction loglike = (1-w)*loglike_i + w*Lbar
#
betainit=c(rep(0,nvar))
#
#  compute pooled optimum
#
out=optim(betainit,llmnl,method="BFGS",control=list( fnscale=-1,trace=0,reltol=1e-6), 
     X=Xpooled,y=ypooled)
betapooled=out$par
H=mnlHess(betapooled,ypooled,Xpooled)
rm(ypooled)
rm(Xpooled)
invisible(gc())
rootH=chol(H)
hess = vector(mode="list", length=nlgt)
for (i in 1:nlgt) 
{
   wgt=length(Data$lgtdata[[i]]$y)/ny
   out=optim(betapooled,llmnlFract,method="BFGS",control=list( fnscale=-1,trace=0,reltol=1e-4), 
   X=Data$lgtdata[[i]]$X,y=Data$lgtdata[[i]]$y,betapooled=betapooled,rootH=rootH,w=w,wgt=wgt)
   if(out$convergence == 0) 
     { hess[[i]]=mnlHess(out$par,Data$lgtdata[[i]]$y,Data$lgtdata[[i]]$X)
     oldbetas[i,] = out$par}
   else
     { hess[[i]]=diag(nvar)
     oldbetas[i,]=c(rep(0,nvar))}
   if(i%%50 ==0) cat("  completed unit #",i,fill=TRUE)
#   fsh()
}
#
#  initialize values
#
# set initial values for the indicators
#     ind is of length(nlgt) and indicates which mixture component this obs
#     belongs to.
#
ind=NULL
ninc=floor(nlgt/ncomp)
for (i in 1:(ncomp-1)) {ind=c(ind,rep(i,ninc))}
if(ncomp != 1) {ind = c(ind,rep(ncomp,nlgt-length(ind)))} else {ind=rep(1,nlgt)}
#
# initialize probs
#
oldprob=rep(1/ncomp,ncomp)
#
#initialize delta
#
if (drawdelta){
  olddelta = rep(0,nz*nvar)
} else { #send placeholders to the _loop function if there is no Z matrix
  olddelta = 0
  Z = matrix(0)
  deltabar = 0
  Ad = matrix(0)
}

###################################################################
# Wayne Taylor
# 09/22/2014
###################################################################
# compdraw =  rhierMnlRwMixtureParallel_rcpp_loop(lgtdata, Z,
#                                      deltabar, Ad, mubar, Amu,
#                                      nu, V, s,
#                                      R, keep, nprint, drawdelta,
#                                      as.matrix(olddelta), a, oldprob, oldbetas, ind,
#                                      Mcmc$B)
draws =  rhierMnlRwMixtureParallel_rcpp_loop(Data$lgtdata, Z,
                                                deltabar, Ad, mubar, Amu,
                                                nu, V, s,
                                                R, keep, nprint, drawdelta,
                                                as.matrix(olddelta), a, oldprob, oldbetas, hess, ind,
                                                Mcmc$B, Mcmc$R2, Data$Ns, Data$N, Mcmc$drawcomp)
####################################################################

# if(drawdelta){
#   attributes(draws$Deltadraw)$class=c("bayesm.mat","mcmc")
#   attributes(draws$Deltadraw)$mcpar=c(1,R,keep)}
# attributes(draws$betadraw)$class=c("bayesm.hcoef")
# attributes(draws$nmix)$class="bayesm.nmix"

# Ndraws = ceiling(Mcmc$R2 * (Data$Ns/Data$N))
# betadraw = matrix(0, nrow=Ndraws, ncol=nvar)
# r = sample(1:length(compdraw), size = Ndraws, replace = TRUE)
# for (i in 1:Ndraws)
# {
#   root = backsolve(compdraw[[r[i]]][[1]]$rooti, diag(nvar))
#   betadraw[i,] = compdraw[[r[i]]][[1]]$mu + t(root) %*% matrix(rnorm(nvar), ncol=1)
# } 
# return(list(betadraw = betadraw, R = Mcmc$R2))

return(draws)
}