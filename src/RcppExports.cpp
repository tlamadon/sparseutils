// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// sparseDiagCross
NumericVector sparseDiagCross(S4 A, NumericMatrix M, S4 B, int n);
RcppExport SEXP _sparseutils_sparseDiagCross(SEXP ASEXP, SEXP MSEXP, SEXP BSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type M(MSEXP);
    Rcpp::traits::input_parameter< S4 >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(sparseDiagCross(A, M, B, n));
    return rcpp_result_gen;
END_RCPP
}
// sparseDiagDot
NumericVector sparseDiagDot(S4 A, NumericVector D, S4 B, int n);
RcppExport SEXP _sparseutils_sparseDiagDot(SEXP ASEXP, SEXP DSEXP, SEXP BSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type D(DSEXP);
    Rcpp::traits::input_parameter< S4 >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(sparseDiagDot(A, D, B, n));
    return rcpp_result_gen;
END_RCPP
}
// sparseSandwichTrace
NumericVector sparseSandwichTrace(S4 L, NumericMatrix M, S4 Q, S4 R, int n, int n0);
RcppExport SEXP _sparseutils_sparseSandwichTrace(SEXP LSEXP, SEXP MSEXP, SEXP QSEXP, SEXP RSEXP, SEXP nSEXP, SEXP n0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type L(LSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type M(MSEXP);
    Rcpp::traits::input_parameter< S4 >::type Q(QSEXP);
    Rcpp::traits::input_parameter< S4 >::type R(RSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type n0(n0SEXP);
    rcpp_result_gen = Rcpp::wrap(sparseSandwichTrace(L, M, Q, R, n, n0));
    return rcpp_result_gen;
END_RCPP
}
// sparseSandwichDiag
NumericVector sparseSandwichDiag(S4 L, NumericMatrix M, S4 Q, S4 R, int n, int n0);
RcppExport SEXP _sparseutils_sparseSandwichDiag(SEXP LSEXP, SEXP MSEXP, SEXP QSEXP, SEXP RSEXP, SEXP nSEXP, SEXP n0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type L(LSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type M(MSEXP);
    Rcpp::traits::input_parameter< S4 >::type Q(QSEXP);
    Rcpp::traits::input_parameter< S4 >::type R(RSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type n0(n0SEXP);
    rcpp_result_gen = Rcpp::wrap(sparseSandwichDiag(L, M, Q, R, n, n0));
    return rcpp_result_gen;
END_RCPP
}
// denseTraceProd
double denseTraceProd(NumericMatrix A, NumericMatrix B);
RcppExport SEXP _sparseutils_denseTraceProd(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(denseTraceProd(A, B));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sparseutils_sparseDiagCross", (DL_FUNC) &_sparseutils_sparseDiagCross, 4},
    {"_sparseutils_sparseDiagDot", (DL_FUNC) &_sparseutils_sparseDiagDot, 4},
    {"_sparseutils_sparseSandwichTrace", (DL_FUNC) &_sparseutils_sparseSandwichTrace, 6},
    {"_sparseutils_sparseSandwichDiag", (DL_FUNC) &_sparseutils_sparseSandwichDiag, 6},
    {"_sparseutils_denseTraceProd", (DL_FUNC) &_sparseutils_denseTraceProd, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_sparseutils(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
