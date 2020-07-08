// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// cumcifun
Rcpp::List cumcifun(Rcpp::List& param, bool T_F, double lower, double upper);
RcppExport SEXP _art_cumcifun(SEXP paramSEXP, SEXP T_FSEXP, SEXP lowerSEXP, SEXP upperSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type param(paramSEXP);
    Rcpp::traits::input_parameter< bool >::type T_F(T_FSEXP);
    Rcpp::traits::input_parameter< double >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< double >::type upper(upperSEXP);
    rcpp_result_gen = Rcpp::wrap(cumcifun(param, T_F, lower, upper));
    return rcpp_result_gen;
END_RCPP
}
// cumcicurve
Eigen::VectorXd cumcicurve(Rcpp::List& param, bool T_F, double lower, double upper, int npar);
RcppExport SEXP _art_cumcicurve(SEXP paramSEXP, SEXP T_FSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP nparSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type param(paramSEXP);
    Rcpp::traits::input_parameter< bool >::type T_F(T_FSEXP);
    Rcpp::traits::input_parameter< double >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< double >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< int >::type npar(nparSEXP);
    rcpp_result_gen = Rcpp::wrap(cumcicurve(param, T_F, lower, upper, npar));
    return rcpp_result_gen;
END_RCPP
}
// eval_incirate
Eigen::VectorXd eval_incirate(Rcpp::List& param, bool T_F, Eigen::VectorXd t);
RcppExport SEXP _art_eval_incirate(SEXP paramSEXP, SEXP T_FSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type param(paramSEXP);
    Rcpp::traits::input_parameter< bool >::type T_F(T_FSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(eval_incirate(param, T_F, t));
    return rcpp_result_gen;
END_RCPP
}
// eval_accuracy
Eigen::VectorXd eval_accuracy(Rcpp::List& param, Eigen::VectorXd t);
RcppExport SEXP _art_eval_accuracy(SEXP paramSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(eval_accuracy(param, t));
    return rcpp_result_gen;
END_RCPP
}
// rcppeigen_hello_world
Eigen::MatrixXd rcppeigen_hello_world();
RcppExport SEXP _art_rcppeigen_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcppeigen_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcppeigen_outerproduct
Eigen::MatrixXd rcppeigen_outerproduct(const Eigen::VectorXd& x);
RcppExport SEXP _art_rcppeigen_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcppeigen_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcppeigen_innerproduct
double rcppeigen_innerproduct(const Eigen::VectorXd& x);
RcppExport SEXP _art_rcppeigen_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcppeigen_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcppeigen_bothproducts
Rcpp::List rcppeigen_bothproducts(const Eigen::VectorXd& x);
RcppExport SEXP _art_rcppeigen_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcppeigen_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_art_cumcifun", (DL_FUNC) &_art_cumcifun, 4},
    {"_art_cumcicurve", (DL_FUNC) &_art_cumcicurve, 5},
    {"_art_eval_incirate", (DL_FUNC) &_art_eval_incirate, 3},
    {"_art_eval_accuracy", (DL_FUNC) &_art_eval_accuracy, 2},
    {"_art_rcppeigen_hello_world", (DL_FUNC) &_art_rcppeigen_hello_world, 0},
    {"_art_rcppeigen_outerproduct", (DL_FUNC) &_art_rcppeigen_outerproduct, 1},
    {"_art_rcppeigen_innerproduct", (DL_FUNC) &_art_rcppeigen_innerproduct, 1},
    {"_art_rcppeigen_bothproducts", (DL_FUNC) &_art_rcppeigen_bothproducts, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_art(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}