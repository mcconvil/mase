// [[Rcpp::depends(RcppEigen)]]

#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Eigen;
using namespace Rcpp;

// [[Rcpp::export]]
SEXP get_weights_greg(Eigen::Map<Eigen::MatrixXd> X_pop,
                      Eigen::Map<Eigen::MatrixXd> X_samp,
                      Eigen::Map<Eigen::MatrixXd> W,
                      Eigen::Map<Eigen::MatrixXd> one_mat){
  
  Eigen::MatrixXd p = one_mat + ((X_pop - X_samp.transpose() * W.diagonal()).transpose() * (X_samp.transpose() * W * X_samp).inverse() * X_samp.transpose());
  Eigen::MatrixXd out = p * W;
  
  return Rcpp::wrap(out);
}


// [[Rcpp::export]]
SEXP get_weights_modGreg(Eigen::Map<Eigen::MatrixXd> X_pop_dom,
                         Eigen::Map<Eigen::MatrixXd> X_samp_dom,
                         Eigen::Map<Eigen::MatrixXd> W_dom,
                         Eigen::Map<Eigen::MatrixXd> const1,
                         Eigen::Map<Eigen::MatrixXd> const2,
                         Eigen::Map<Eigen::MatrixXd> weighted_indic_mat){
  
  Eigen::MatrixXd out = weighted_indic_mat + (((X_pop_dom - X_samp_dom.transpose() * W_dom.diagonal()).transpose() * const1) * const2);
  
  return Rcpp::wrap(out);
}

