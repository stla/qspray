// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/qspray.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// lexLeadingArma
unsigned int lexLeadingArma(const arma::umat& M);
RcppExport SEXP _qspray_lexLeadingArma(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::umat& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(lexLeadingArma(M));
    return rcpp_result_gen;
END_RCPP
}
// qsprayDivisionRcpp
Rcpp::List qsprayDivisionRcpp(Rcpp::List Powers1, Rcpp::StringVector coeffs1, Rcpp::List Powers2, Rcpp::StringVector coeffs2);
RcppExport SEXP _qspray_qsprayDivisionRcpp(SEXP Powers1SEXP, SEXP coeffs1SEXP, SEXP Powers2SEXP, SEXP coeffs2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type Powers1(Powers1SEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type coeffs1(coeffs1SEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type Powers2(Powers2SEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type coeffs2(coeffs2SEXP);
    rcpp_result_gen = Rcpp::wrap(qsprayDivisionRcpp(Powers1, coeffs1, Powers2, coeffs2));
    return rcpp_result_gen;
END_RCPP
}
// BBdivisionRcpp
Rcpp::List BBdivisionRcpp(Rcpp::List Powers, Rcpp::StringVector coeffs, Rcpp::List gs, Rcpp::List LTgs, int d);
RcppExport SEXP _qspray_BBdivisionRcpp(SEXP PowersSEXP, SEXP coeffsSEXP, SEXP gsSEXP, SEXP LTgsSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type Powers(PowersSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type coeffs(coeffsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type gs(gsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type LTgs(LTgsSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(BBdivisionRcpp(Powers, coeffs, gs, LTgs, d));
    return rcpp_result_gen;
END_RCPP
}
// evalQxspray
Rcpp::StringVector evalQxspray(const Rcpp::List Powers, const Rcpp::StringVector coeffs, const Rcpp::StringVector v_re, const Rcpp::StringVector v_im);
RcppExport SEXP _qspray_evalQxspray(SEXP PowersSEXP, SEXP coeffsSEXP, SEXP v_reSEXP, SEXP v_imSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type Powers(PowersSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector >::type coeffs(coeffsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector >::type v_re(v_reSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector >::type v_im(v_imSEXP);
    rcpp_result_gen = Rcpp::wrap(evalQxspray(Powers, coeffs, v_re, v_im));
    return rcpp_result_gen;
END_RCPP
}
// qspray_maker
Rcpp::List qspray_maker(const Rcpp::List& Powers, const Rcpp::StringVector& coeffs);
RcppExport SEXP _qspray_qspray_maker(SEXP PowersSEXP, SEXP coeffsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Powers(PowersSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type coeffs(coeffsSEXP);
    rcpp_result_gen = Rcpp::wrap(qspray_maker(Powers, coeffs));
    return rcpp_result_gen;
END_RCPP
}
// qspray_deriv
Rcpp::List qspray_deriv(const Rcpp::List& Powers, const Rcpp::StringVector& coeffs, const Rcpp::IntegerVector& n);
RcppExport SEXP _qspray_qspray_deriv(SEXP PowersSEXP, SEXP coeffsSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Powers(PowersSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type coeffs(coeffsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(qspray_deriv(Powers, coeffs, n));
    return rcpp_result_gen;
END_RCPP
}
// qspray_add
Rcpp::List qspray_add(const Rcpp::List& Powers1, const Rcpp::StringVector& coeffs1, const Rcpp::List& Powers2, const Rcpp::StringVector& coeffs2);
RcppExport SEXP _qspray_qspray_add(SEXP Powers1SEXP, SEXP coeffs1SEXP, SEXP Powers2SEXP, SEXP coeffs2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Powers1(Powers1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type coeffs1(coeffs1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Powers2(Powers2SEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type coeffs2(coeffs2SEXP);
    rcpp_result_gen = Rcpp::wrap(qspray_add(Powers1, coeffs1, Powers2, coeffs2));
    return rcpp_result_gen;
END_RCPP
}
// qspray_subtract
Rcpp::List qspray_subtract(const Rcpp::List& Powers1, const Rcpp::StringVector& coeffs1, const Rcpp::List& Powers2, const Rcpp::StringVector& coeffs2);
RcppExport SEXP _qspray_qspray_subtract(SEXP Powers1SEXP, SEXP coeffs1SEXP, SEXP Powers2SEXP, SEXP coeffs2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Powers1(Powers1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type coeffs1(coeffs1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Powers2(Powers2SEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type coeffs2(coeffs2SEXP);
    rcpp_result_gen = Rcpp::wrap(qspray_subtract(Powers1, coeffs1, Powers2, coeffs2));
    return rcpp_result_gen;
END_RCPP
}
// qspray_mult
Rcpp::List qspray_mult(const Rcpp::List& Powers1, const Rcpp::StringVector& coeffs1, const Rcpp::List& Powers2, const Rcpp::StringVector& coeffs2);
RcppExport SEXP _qspray_qspray_mult(SEXP Powers1SEXP, SEXP coeffs1SEXP, SEXP Powers2SEXP, SEXP coeffs2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Powers1(Powers1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type coeffs1(coeffs1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Powers2(Powers2SEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type coeffs2(coeffs2SEXP);
    rcpp_result_gen = Rcpp::wrap(qspray_mult(Powers1, coeffs1, Powers2, coeffs2));
    return rcpp_result_gen;
END_RCPP
}
// qspray_equality
bool qspray_equality(const Rcpp::List& Powers1, const Rcpp::StringVector& coeffs1, const Rcpp::List& Powers2, const Rcpp::StringVector& coeffs2);
RcppExport SEXP _qspray_qspray_equality(SEXP Powers1SEXP, SEXP coeffs1SEXP, SEXP Powers2SEXP, SEXP coeffs2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Powers1(Powers1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type coeffs1(coeffs1SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Powers2(Powers2SEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type coeffs2(coeffs2SEXP);
    rcpp_result_gen = Rcpp::wrap(qspray_equality(Powers1, coeffs1, Powers2, coeffs2));
    return rcpp_result_gen;
END_RCPP
}
// qspray_power
Rcpp::List qspray_power(const Rcpp::List& Powers, const Rcpp::StringVector& coeffs, unsigned int n);
RcppExport SEXP _qspray_qspray_power(SEXP PowersSEXP, SEXP coeffsSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Powers(PowersSEXP);
    Rcpp::traits::input_parameter< const Rcpp::StringVector& >::type coeffs(coeffsSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(qspray_power(Powers, coeffs, n));
    return rcpp_result_gen;
END_RCPP
}
// lexLeadingIndexCPP
int lexLeadingIndexCPP(const Rcpp::List& Powers);
RcppExport SEXP _qspray_lexLeadingIndexCPP(SEXP PowersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Powers(PowersSEXP);
    rcpp_result_gen = Rcpp::wrap(lexLeadingIndexCPP(Powers));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_qspray_lexLeadingArma", (DL_FUNC) &_qspray_lexLeadingArma, 1},
    {"_qspray_qsprayDivisionRcpp", (DL_FUNC) &_qspray_qsprayDivisionRcpp, 4},
    {"_qspray_BBdivisionRcpp", (DL_FUNC) &_qspray_BBdivisionRcpp, 5},
    {"_qspray_evalQxspray", (DL_FUNC) &_qspray_evalQxspray, 4},
    {"_qspray_qspray_maker", (DL_FUNC) &_qspray_qspray_maker, 2},
    {"_qspray_qspray_deriv", (DL_FUNC) &_qspray_qspray_deriv, 3},
    {"_qspray_qspray_add", (DL_FUNC) &_qspray_qspray_add, 4},
    {"_qspray_qspray_subtract", (DL_FUNC) &_qspray_qspray_subtract, 4},
    {"_qspray_qspray_mult", (DL_FUNC) &_qspray_qspray_mult, 4},
    {"_qspray_qspray_equality", (DL_FUNC) &_qspray_qspray_equality, 4},
    {"_qspray_qspray_power", (DL_FUNC) &_qspray_qspray_power, 3},
    {"_qspray_lexLeadingIndexCPP", (DL_FUNC) &_qspray_lexLeadingIndexCPP, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_qspray(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
