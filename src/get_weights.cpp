// [[Rcpp::depends(RcppEigen)]]

#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Eigen;
using namespace Rcpp;

// [[Rcpp::export]]
SEXP get_weights(Eigen::Map<Eigen::MatrixXd> X_pop,
                 Eigen::Map<Eigen::MatrixXd> X_samp,
                 Eigen::Map<Eigen::MatrixXd> W,
                 Eigen::Map<Eigen::MatrixXd> one_mat){
  
  Eigen::MatrixXd p = one_mat + ((X_pop - X_samp.transpose() * W.diagonal()).transpose() * (X_samp.transpose() * W * X_samp).inverse() * X_samp.transpose());
  Eigen::MatrixXd out = p * W;
  
  return Rcpp::wrap(out);
}
