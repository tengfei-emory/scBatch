//#include <Rcpp.h>

#include <RcppArmadillo.h>
// [[Rcpp::depends( RcppArmadillo)]]
using namespace arma;

//' Correct RNA-seq count matrix subject to batch effects by sample distance matrix correction (faster Rcpp version)
//'
//' @description  Use gradient descent algorithm to correct count matrix by linear transformation based on known sample matrix correction.
//' @param c The original p*n batch effect data with n subjects and p RNA-seq measurements.
//' @param d The n*n distance matrix obtained by QuantNorm.
//' @param w An initial n*n weight matrix to conduct linear transformation. Default to be identity matrix if not specified.
//' @param m Number of groups to be divided for coordinate gradient descent. 1 < m <= n. Default to be 0.1n if not specified.
//' @param max Maximum number of the iteration if the tolerance is not reached.
//' @param step Step size of the gradient descent algorithm.
//' @param tol Stopping criteria of the algorithm. The algorithm stops if the step size is smaller than tol.
//' @param derif Function to compute gradient of the loss function.
//' @return Returns the corrected count matrix.
//' @author Teng Fei. Email: tfei@emory.edu
//' @references Fei et al (2018), Mitigating the adverse impact of batch effects in sample pattern detection, Bioinformatics, <https://doi.org/10.1093/bioinformatics/bty117>.
//' @export
// [[Rcpp::export]]
arma::mat scBatchCpp(arma::mat c, arma::mat w, arma::mat d, int m, double max, double step, double tol, Rcpp::Function derif){
  //NumericMatrix d, NumericMatrix core, NumericVector idx) {
  arma::wall_clock timer;

  int p = c.n_rows; int n = c.n_cols;
  arma::mat core = c.t()*(eye(p,p)-ones(p,p)/p)*(eye(p,p)-ones(p,p)/p)*c;
  //core = t(count.mat)%*%t((diag(p) - matrix(1,p,p)/p))%*%(diag(p) - matrix(1,p,p)/p)%*%count.mat
  core = (core+core.t())/2;

  for(int i = 0; i < max; ++i) {
    arma::vec group = randi<vec>(n,arma::distr_param(0,m-1));
    for(int k = 0; k < m; ++k){
      timer.tic();
      arma::uvec idx = arma::find(group == k);
      Rcpp::List fdf = derif(c,w,d,core,idx);
      double f = fdf["f"];
      arma::mat df = fdf["df"];
      for(int j = 0; j < 5; ++j){
        arma::mat u = w - step*df;
        arma::vec ln(n); ln.fill(1);
        u = u/(ln*(arma::max(arma::abs(u),0)));
        arma::mat nc = c*u;
        arma::mat A = 1-cor(nc);
        double fnew = accu(pow(A-d,2));

        if (fnew >= f){
          step = 0.5*step;
        }else{
          step = 1.5*step;
          w = u;
          double n = timer.toc();

          cout << k << " time elapsed: " << n << " L: " << fnew << " step size: " << step << endl;
          //cout << fnew << endl;
          break;
        }
      }
    }
    if (step < tol){
      break;
    }
  }

  //Rcpp::List ret;
  //ret["f"] = f;
  //ret["df"] = df;
  //return ret;

  arma::mat nc = c*w;

  return nc;
}
