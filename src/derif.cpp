// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
//' Column gradient for the loss function based on Pearson correlation distance (faster Rcpp version)
//'
//' @description A specific implementation of column gradient to be used in the main coordinate gradient descent algorithm scBatchCpp
//' @param c The original p*n batch effect data with n subjects and p RNA-seq measurements.
//' @param d The n*n distance matrix obtained by QuantNorm.
//' @param w An initial n*n weight matrix to conduct linear transformation. Default to be identity matrix if not specified.
//' @param core Pre-calculated component from main algorithm in order to save computation time.
//' @param idx Random selected columns input from the main algorithm.
//' @author Teng Fei. Email: tfei@emory.edu
//' @export
//' @import RcppArmadillo
// [[Rcpp::export]]
Rcpp::List derif(arma::mat c, arma::mat w, arma::mat d, arma::mat core, arma::uvec idx){
  int n = c.n_cols;
  arma::mat nc = c * w;
  arma::mat a = 1-cor(nc);
  arma::mat amd = a-d;
  double f; f = accu(pow(amd,2));

  arma::mat W = w.t() * core * w;
  arma::uvec idxw = arma::find(W < 0);
  W(idxw) = W(idxw)- W(idxw);
  W = sqrt((W+W.t())/2);
  arma::mat M = core * w;
  arma::mat invW = 1/W;
  arma::mat invWdiag = invW.diag();
  arma::vec ln(n); ln.fill(1);

  arma::mat a1 = M/(ln*W.diag().t());
  arma::vec b1 = vectorise(M/(ln*pow(W.diag(),3).t()));
  arma::mat b20 = pow(W,2)/(ln*pow(W.diag(),2).t());
  arma::vec b2 = vectorise(b20.t());

  arma::mat df(n,n); df.fill(0);
  int lidx = idx.n_elem;

  for(int k = 0; k < lidx; ++k) {
    int i = idx(k)+1;
    arma::mat a11 = invWdiag(i-1)*a1;
    arma::mat D = a11 - (b1.subvec(n*(i-1),n*(i)-1)*b2.subvec(n*(i-1),n*(i)-1).t());
    df.col(i-1) = D.t() * amd.col(i-1);
    arma::mat damd = D.t() * amd;
    df(i-1,i-1) = sum(damd.diag());
  }

  df = -df;

  Rcpp::List ret;
  ret["f"] = f;
  ret["df"] = df;
  return ret;
}
