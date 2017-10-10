#include <Rcpp.h>
#include "knock_interface.h"

using namespace Rcpp;
using namespace knockoffs;

//' Wrapper for DMC knockoffs
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector knockoffDMC_wrapper(SEXP X_, SEXP pInit_, SEXP Q_, SEXP n_, SEXP p_, SEXP K_, SEXP seed_) {
  int n = as<int>(n_);
  int p = as<int>(p_);
  int K = as<int>(K_);
  ivector2 X = numToIntVec2(X_, n);
  vector pInit = numToVec(pInit_);
  vector3 Q = numToVec3(Q_, p-1, K);

  int seed = as<int>(seed_);
  KnockoffDMC knockoffs(pInit, Q, seed);
  ivector2 Xk = knockoffs.sample(X);

  // Convert Xk back to R matrix
  NumericVector Xk_ = NumericVector( Dimension(n,p));
  for(int i=0; i<n; i++) {
    for(int j=0; j<p; j++) {
      Xk_[j*n+i] = Xk[i][j];
    }
  }
  return(Xk_);
}

//' Wrapper for HMM knockoffs
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector knockoffHMM_wrapper(SEXP X_, SEXP pInit_, SEXP Q_, SEXP pEmit_, SEXP n_, SEXP p_, SEXP K_, SEXP M_, SEXP seed_) {
  int n = as<int>(n_);
  int p = as<int>(p_);
  int K = as<int>(K_);
  int M = as<int>(M_);
  ivector2 X = numToIntVec2(X_, n);
  vector pInit = numToVec(pInit_);
  vector3 Q = numToVec3(Q_, p-1, K);
  vector3 pEmit = numToVec3(pEmit_, p, M);

  int seed = as<int>(seed_);
  KnockoffHMM knockoffs(pInit, Q, pEmit, seed);
  ivector2 Xk = knockoffs.sample(X);

  // Convert Xk back to R matrix
  NumericVector Xk_ = NumericVector( Dimension(n,p));
  for(int i=0; i<n; i++) {
    for(int j=0; j<p; j++) {
      Xk_[j*n+i] = Xk[i][j];
    }
  }
  return(Xk_);
}
