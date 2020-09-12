#include "bayesm.h"
  
//[[Rcpp::export]]
List rhierMnlRwMixtureParallel_rcpp_loop(List const& lgtdata, mat const& Z,
                                  vec const& deltabar, mat const& Ad, mat const& mubar, mat const& Amu,
                                  int const& nu, mat const& V, double s,
                                  int R, int keep, int nprint, bool drawdelta,
                                  mat olddelta,  vec const& a, vec oldprob, mat oldbetas, List const& hess, vec ind,
                                  int B, int R2, int Ns, int N, bool drawcomp){

// Wayne Taylor 10/01/2014
// Rico Bumbaca 10/10/2019

  int nlgt = lgtdata.size();
  int nvar = V.n_cols;
  int nz = Z.n_cols;
 
  mat rootpi, betabar, ucholinv, incroot;
  int mkeep;
  mnlMetropOnceOut metropout_struct;
  List lgtdatai, nmix;
  
  // convert List to std::vector of struct
  std::vector<moments> lgtdata_vector;
  moments lgtdatai_struct;
  for (int lgt = 0; lgt<nlgt; lgt++){
    lgtdatai = lgtdata[lgt];
    
    lgtdatai_struct.y = as<vec>(lgtdatai["y"]);
    lgtdatai_struct.X = as<mat>(lgtdatai["X"]);
    lgtdatai_struct.hess = as<mat>(hess[lgt]);
    lgtdata_vector.push_back(lgtdatai_struct);    
  }
   
  // allocate space for draws
  vec oldll = zeros<vec>(nlgt);
//  cube betadraw(nlgt, nvar, floor(R/keep));
//  mat probdraw(floor(R/keep), oldprob.size());
  mat probdraw(floor((R-B)/keep), oldprob.size()); 
//  vec loglike(floor(R/keep));
//  mat Deltadraw(1,1); if(drawdelta) Deltadraw.zeros(floor(R/keep), nz*nvar);//enlarge Deltadraw only if the space is required
  List compdraw(floor((R-B)/keep));
  
  if (nprint>0) startMcmcTimer();
    
  for (int rep = 0; rep<R; rep++){
    
    //first draw comps,ind,p | {beta_i}, delta
    // ind,p need initialization comps is drawn first in sub-Gibbs
    List mgout;
    if(drawdelta) {
      olddelta.reshape(nvar,nz);
      mgout = rmixGibbs (oldbetas-Z*trans(olddelta),mubar,Amu,nu,V,a,oldprob,ind);
    } else {
      mgout = rmixGibbs(oldbetas,mubar,Amu,nu,V,a,oldprob,ind);
    }
    
    List oldcomp = mgout["comps"];
    oldprob = as<vec>(mgout["p"]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
    ind = as<vec>(mgout["z"]);
    
    //now draw delta | {beta_i}, ind, comps
    if(drawdelta) olddelta = drawDelta(Z,oldbetas,ind,oldcomp,deltabar,Ad);
    
    //loop over all LGT equations drawing beta_i | ind[i],z[i,],mu[ind[i]],rooti[ind[i]]
      for(int lgt = 0; lgt<nlgt; lgt++){
        List oldcomplgt = oldcomp[ind[lgt]-1];
        rootpi = as<mat>(oldcomplgt[1]);
        
        //note: beta_i = Delta*z_i + u_i  Delta is nvar x nz
        if(drawdelta){
          olddelta.reshape(nvar,nz);
          betabar = as<vec>(oldcomplgt[0])+olddelta*vectorise(Z(lgt,span::all));
        } else {
          betabar = as<vec>(oldcomplgt[0]);
        }
        
        if (rep == 0) oldll[lgt] = llmnl(vectorise(oldbetas(lgt,span::all)),lgtdata_vector[lgt].y,lgtdata_vector[lgt].X);
        
        //compute inc.root
        ucholinv = solve(trimatu(chol(lgtdata_vector[lgt].hess+rootpi*trans(rootpi))), eye(nvar,nvar)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
        incroot = chol(ucholinv*trans(ucholinv));
                
        metropout_struct = mnlMetropOnce(lgtdata_vector[lgt].y,lgtdata_vector[lgt].X,vectorise(oldbetas(lgt,span::all)),
                                         oldll[lgt],s,incroot,betabar,rootpi);
         
         oldbetas(lgt,span::all) = trans(metropout_struct.betadraw);
         oldll[lgt] = metropout_struct.oldll;  
      }
      
    //print time to completion and draw # every nprint'th draw
    if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);
    
    if((rep+1-B)>0 & (rep+1-B)%keep==0){
      mkeep = (rep+1-B)/keep;
//      betadraw.slice(mkeep-1) = oldbetas;
      probdraw(mkeep-1, span::all) = trans(oldprob);
//      loglike[mkeep-1] = sum(oldll);
//      if(drawdelta) Deltadraw(mkeep-1, span::all) = trans(vectorise(olddelta));
      compdraw[mkeep-1] = oldcomp;
    }
  }
  
  if (nprint>0) endMcmcTimer();

  //draw from posterior predictive density
  int ndraws = ceil((float)R2 * (((float)Ns)/((float)N)));
  mat betadraw(ndraws, nvar);
  int ncompdraw = floor((R-B)/keep);
  ivec r = randi(ndraws, distr_param(0,ncompdraw-1));
  
  List compdrawi, compdrawi0;
  mat root;
  vec mu;
  int indp;

  for(int i = 0; i<ndraws; i++){
    compdrawi = compdraw[r[i]];
    indp = rmultinomF(trans(probdraw(r[i], span::all)));
    compdrawi0 = compdrawi[indp-1]; 
    mat rooti = compdrawi0["rooti"];
    vec mu = compdrawi0["mu"];
    root = solve(trimatu(rooti), eye(nvar,nvar));
    betadraw.row(i) = trans(mu + trans(root) * as<vec>(rnorm(nvar)));
  }
  
  if (drawcomp) {
    return(List::create(
        Named("betadraw") = betadraw,
        Named("compdraw") = compdraw,
        Named("probdraw") = probdraw,
        Named("R") = R2));
  } else {
    return(List::create(
        Named("betadraw") = betadraw,
        Named("R") = R2));
  }
}
