#include "bayesm.h"

vec llmnlParallel(mat const& draws, vec const& y, mat const& X);

//[[Rcpp::export]]
List IndepMetropParallel_rcpp_loop(List const& data, mat const& draws, 
                                   List const& mcmc, List const& options) {
  
// Rico Bumbaca 10/10/2019
  
  int R, k, N, keep, mkeep, P;
  keep = mcmc["keep2"];
  R = draws.n_rows;
  k = draws.n_cols;
  P = options["P"];
  
  int nunits;
  nunits = data["Nsunits"];
  
  List lgtdata = data["lgtdata"];
  List lgtdatai;
  
  N = lgtdata.length();
  //  cube betadraw(N, k, floor(R/keep));
  cube betadraw(nunits, k, floor(R/keep));
  
  rowvec betainit(k);
  betainit = mean(draws, 0);
  rowvec betap(k);
  
  if (options["vectorize"]){
    int n = N/P;   // N must be a multiple of P
    vec ratio(n);    
    uvec ind(n);
    vec unif(n);
    mat llratio(n,2);
    llratio.col(1) = zeros<vec>(n);
    
    vec logc(n), logp(n);
    mat currentdraw(n, k);
    mat loglike(n, R);
    for (int p=0; p<P; p++) {
      mat betac = repmat(betainit, n, 1);
      for (int i=0; i<n; i++) {
        lgtdatai = lgtdata[p * n + i];
        vec y = lgtdatai["y"];
        mat X = lgtdatai["X"];
        logc[i] = llmnl(trans(betainit), y, X);
        loglike.row(i) = trans(llmnlParallel(trans(draws), y, X));
      }
      for (int rep=0; rep<R; rep++) {
        betap = draws.row(rep);
        logp = loglike.col(rep);
//        ratio = exp(logp - logc);
        llratio.col(0) = logp - logc;
        ratio = exp(min(llratio, 1));
        unif = runif(n);
        ind = find((unif < ratio) == 1);
        currentdraw = betac;
        if (ind.n_elem != 0) {
          for (int i=0; i<ind.n_elem; i++) {
            currentdraw.row(ind[i]) = betap;
            betac.row(ind[i]) = betap;
            //        logc(ind[i]) = logp(ind[i]);
          }
          logc(ind) = logp(ind);
        }
        if((rep+1)>0 & (rep+1)%keep==0){
          mkeep = (rep+1)/keep;
          //        betadraw.subcube(span(p*n,p*n+n-1), span(), span(mkeep-1,mkeep-1)) = currentdraw;
          if (nunits>0) betadraw.slice(mkeep-1) = currentdraw.rows(0,nunits-1);
        }
      }
    }
  } else {
    double ratio, unif;
    mat currentdraw(N, k);
    vec logc(N), logp(N);
    mat betac = repmat(betainit, N, 1);
    for (int i=0; i<N; i++){
      lgtdatai = lgtdata[i];
      vec y = lgtdatai["y"];
      mat X = lgtdatai["X"];
      logc(i) = llmnl(trans(betainit), y, X);
    }
    for (int rep=0; rep<R; rep++) {
      betap = draws.row(rep);
      for (int i=0; i<N; i++) {
        lgtdatai = lgtdata[i];
        vec y = lgtdatai["y"];
        mat X = lgtdatai["X"];
        logp(i) = llmnl(trans(betap), y, X);
        if (logp(i) - logc(i) > 0) 
          ratio = 1;
        else 
          ratio = exp(logp(i) - logc(i));
        unif = runif(1)[0];
        if (unif < ratio) {
          currentdraw.row(i) = betap;
          betac.row(i) = betap;
          logc(i) = logp(i);
        } else {
          currentdraw.row(i) = betac.row(i);
        }
      }
      if((rep+1)>0 & (rep+1)%keep==0){
        mkeep = (rep+1)/keep;
        //        betadraw.subcube(span(), span(), span(mkeep-1,mkeep-1)) = currentdraw;
        if (nunits>0) betadraw.slice(mkeep-1) = currentdraw.rows(0,nunits-1);
      }
    }
  }
  if (nunits>0){
    return(List::create(
        Named("betadraw") = betadraw));
  } else {
    return(List::create(
        Named("betadraw") = R_NilValue)); //sets the value to NULL in R
  }
}


