// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// get_weights_greg
SEXP get_weights_greg(Eigen::Map<Eigen::MatrixXd> X_pop, Eigen::Map<Eigen::MatrixXd> X_samp, Eigen::Map<Eigen::MatrixXd> W, Eigen::Map<Eigen::MatrixXd> one_mat);
RcppExport SEXP _mase_get_weights_greg(SEXP X_popSEXP, SEXP X_sampSEXP, SEXP WSEXP, SEXP one_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X_pop(X_popSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X_samp(X_sampSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type W(WSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type one_mat(one_matSEXP);
    rcpp_result_gen = Rcpp::wrap(get_weights_greg(X_pop, X_samp, W, one_mat));
    return rcpp_result_gen;
END_RCPP
}
// get_weights_modGreg
SEXP get_weights_modGreg(Eigen::Map<Eigen::MatrixXd> X_pop_dom, Eigen::Map<Eigen::MatrixXd> X_samp_dom, Eigen::Map<Eigen::MatrixXd> W_dom, Eigen::Map<Eigen::MatrixXd> const1, Eigen::Map<Eigen::MatrixXd> const2, Eigen::Map<Eigen::MatrixXd> weighted_indic_mat);
RcppExport SEXP _mase_get_weights_modGreg(SEXP X_pop_domSEXP, SEXP X_samp_domSEXP, SEXP W_domSEXP, SEXP const1SEXP, SEXP const2SEXP, SEXP weighted_indic_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X_pop_dom(X_pop_domSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X_samp_dom(X_samp_domSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type W_dom(W_domSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type const1(const1SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type const2(const2SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type weighted_indic_mat(weighted_indic_matSEXP);
    rcpp_result_gen = Rcpp::wrap(get_weights_modGreg(X_pop_dom, X_samp_dom, W_dom, const1, const2, weighted_indic_mat));
    return rcpp_result_gen;
END_RCPP
}
// const_comp1
SEXP const_comp1(Eigen::Map<Eigen::MatrixXd> X_samp, Eigen::Map<Eigen::MatrixXd> W);
RcppExport SEXP _mase_const_comp1(SEXP X_sampSEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X_samp(X_sampSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(const_comp1(X_samp, W));
    return rcpp_result_gen;
END_RCPP
}
// get_coefs
SEXP get_coefs(Eigen::Map<Eigen::MatrixXd> X_samp, Eigen::Map<Eigen::VectorXd> Y, Eigen::Map<Eigen::MatrixXd> W);
RcppExport SEXP _mase_get_coefs(SEXP X_sampSEXP, SEXP YSEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X_samp(X_sampSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(get_coefs(X_samp, Y, W));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mase_get_weights_greg", (DL_FUNC) &_mase_get_weights_greg, 4},
    {"_mase_get_weights_modGreg", (DL_FUNC) &_mase_get_weights_modGreg, 6},
    {"_mase_const_comp1", (DL_FUNC) &_mase_const_comp1, 2},
    {"_mase_get_coefs", (DL_FUNC) &_mase_get_coefs, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_mase(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
