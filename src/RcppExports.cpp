// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/dynmix.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// F
vec F(vec x, int N, mat Rx);
RcppExport SEXP _dynmix_F(SEXP xSEXP, SEXP NSEXP, SEXP RxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< mat >::type Rx(RxSEXP);
    rcpp_result_gen = Rcpp::wrap(F(x, N, Rx));
    return rcpp_result_gen;
END_RCPP
}
// J
mat J(vec x, int m, mat Rxj);
RcppExport SEXP _dynmix_J(SEXP xSEXP, SEXP mSEXP, SEXP RxjSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< mat >::type Rxj(RxjSEXP);
    rcpp_result_gen = Rcpp::wrap(J(x, m, Rxj));
    return rcpp_result_gen;
END_RCPP
}
// newt
vec newt(vec x0, int Nn, mat Rxn);
RcppExport SEXP _dynmix_newt(SEXP x0SEXP, SEXP NnSEXP, SEXP RxnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< int >::type Nn(NnSEXP);
    Rcpp::traits::input_parameter< mat >::type Rxn(RxnSEXP);
    rcpp_result_gen = Rcpp::wrap(newt(x0, Nn, Rxn));
    return rcpp_result_gen;
END_RCPP
}
// mKIapprox
mat mKIapprox(mat w, mat vold);
RcppExport SEXP _dynmix_mKIapprox(SEXP wSEXP, SEXP voldSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type w(wSEXP);
    Rcpp::traits::input_parameter< mat >::type vold(voldSEXP);
    rcpp_result_gen = Rcpp::wrap(mKIapprox(w, vold));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dynmix_F", (DL_FUNC) &_dynmix_F, 3},
    {"_dynmix_J", (DL_FUNC) &_dynmix_J, 3},
    {"_dynmix_newt", (DL_FUNC) &_dynmix_newt, 3},
    {"_dynmix_mKIapprox", (DL_FUNC) &_dynmix_mKIapprox, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_dynmix(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}