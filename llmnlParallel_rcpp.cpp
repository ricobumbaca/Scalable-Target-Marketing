#include "bayesm.h"


//[[Rcpp::export]]
vec llmnlParallel(mat const& beta, vec const& y, mat const& X){
  
// Rico Bumbaca 10/10/19

  // Evaluates log-likelihood for the multinomial logit model for m betas
  int m = beta.n_cols; 
  int n = y.size();    
  int j = X.n_rows/n;  
  mat Xbeta = X * beta;
  int l = n*m;
  Xbeta.reshape(j, l);
  Xbeta = Xbeta.t();
  
  vec xby(l);
  for(int i = 0; i<l; i++) {
    xby[i] = Xbeta.at(i, y[i%n]-1);
  }
  
  // int l = n*m;
  // vec xby(l);
  // myFunctorClass function1(l, n, y);
  //    uvec indices = linspace<uvec>(0, l-1, l);
  //    indices.transform(function1);
  //    xby=Xbeta(indices);
  
  //   int l = n*m;
  //   vec xby(l);
  //   uvec indices(l);
  //   for (int i = 0; i<l; i++) {
  //     indices[i] = (y[i%n]-1)*l+i;
  //   }
  //   xby = Xbeta(indices);
  
  //  imat yi(n,j);
  //  yi.zeros();
  //  for (int i = 0; i<n; i++) {
  //    yi.at(i,y[i%n]-1) = 1;
  //  }
  //  imat yim = conv_to<imat>::from(repmat(conv_to<mat>::from(yi), m, 1));
  //  vec xby(n*m);
  //  xby = sum(Xbeta % yim, 1);
  
  //  imat yi(j,n);
  //  yi.zeros();
  //  for (int i = 0; i<n; i++) {
  //    yi.at(y[i%n]-1,i) = 1;
  //  }
  //  vec xby(n*m);
  //  for(int i = 0; i<(n*m);i++) {
  //    xby[i] = dot(Xbeta.row(i), yi.col(i%n));
  //  }
  
  //  uvec ym = conv_to<uvec>::from(repmat(conv_to<mat>::from(y), m, 1));
  //  vec xby(n*m);
  //  for(int i = 0; i<(n*m);i++) {
  //    xby[i] = Xbeta.at(i,ym[i]-1);
  //  }
  
  //  uvec ym = conv_to<uvec>::from(repmat(conv_to<mat>::from(y), m, 1));
  //  uvec seq = linspace<uvec>(0, (n*m)-1, n*m);
  //  mat XbetaSub = Xbeta.submat( seq, ym-1 );
  //  vec xby = XbetaSub.diag();
  
  Xbeta = exp(Xbeta);
  
  uvec iota(j);
  iota.ones();
  vec denom = log(Xbeta * iota);
  
  //  vec denom = log(sum(Xbeta,1));
  
  mat xby_denom = xby - denom;
  xby_denom.reshape(n, m);
  return(conv_to<vec>::from(sum(xby_denom)));
}
