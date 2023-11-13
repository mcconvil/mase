// [[Rcpp::depends(RcppEigen)]]

#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Eigen;
using namespace Rcpp;

// [[Rcpp::export]]
SEXP get_coefs(Eigen::Map<Eigen::MatrixXd> X_samp,
               Eigen::Map<Eigen::VectorXd> Y,
               Eigen::Map<Eigen::MatrixXd> W){
  
  Eigen::VectorXd out = (X_samp.transpose() * W * X_samp).inverse() * X_samp.transpose() * W * Y;
  
  return Rcpp::wrap(out);
  
}

