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
// eval_log_like
double eval_log_like(Rcpp::List& param, bool T_F, double t);
RcppExport SEXP _art_eval_log_like(SEXP paramSEXP, SEXP T_FSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List& >::type param(paramSEXP);
    Rcpp::traits::input_parameter< bool >::type T_F(T_FSEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(eval_log_like(param, T_F, t));
    return rcpp_result_gen;
END_RCPP
}
// get_loglike
NumericVector get_loglike(NumericMatrix lambda_, NumericMatrix theta_, NumericMatrix z_, NumericMatrix w_, NumericMatrix gamma_, List& param_);
RcppExport SEXP _art_get_loglike(SEXP lambda_SEXP, SEXP theta_SEXP, SEXP z_SEXP, SEXP w_SEXP, SEXP gamma_SEXP, SEXP param_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type lambda_(lambda_SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta_(theta_SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type z_(z_SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type w_(w_SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gamma_(gamma_SEXP);
    Rcpp::traits::input_parameter< List& >::type param_(param_SEXP);
    rcpp_result_gen = Rcpp::wrap(get_loglike(lambda_, theta_, z_, w_, gamma_, param_));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_gen_surv_pp
Eigen::MatrixXd rcpp_gen_surv_pp(NumericMatrix lambda_, NumericMatrix theta_, NumericMatrix z_, NumericMatrix w_, NumericMatrix gamma_, List& param_, int item);
RcppExport SEXP _art_rcpp_gen_surv_pp(SEXP lambda_SEXP, SEXP theta_SEXP, SEXP z_SEXP, SEXP w_SEXP, SEXP gamma_SEXP, SEXP param_SEXP, SEXP itemSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type lambda_(lambda_SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta_(theta_SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type z_(z_SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type w_(w_SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gamma_(gamma_SEXP);
    Rcpp::traits::input_parameter< List& >::type param_(param_SEXP);
    Rcpp::traits::input_parameter< int >::type item(itemSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_gen_surv_pp(lambda_, theta_, z_, w_, gamma_, param_, item));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_art_cumcifun", (DL_FUNC) &_art_cumcifun, 4},
    {"_art_cumcicurve", (DL_FUNC) &_art_cumcicurve, 5},
    {"_art_eval_incirate", (DL_FUNC) &_art_eval_incirate, 3},
    {"_art_eval_accuracy", (DL_FUNC) &_art_eval_accuracy, 2},
    {"_art_eval_log_like", (DL_FUNC) &_art_eval_log_like, 3},
    {"_art_get_loglike", (DL_FUNC) &_art_get_loglike, 6},
    {"_art_rcpp_gen_surv_pp", (DL_FUNC) &_art_rcpp_gen_surv_pp, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_art(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
