// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// setDebug
void setDebug(bool debug_lvl);
RcppExport SEXP STASNet_setDebug(SEXP debug_lvlSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< bool >::type debug_lvl(debug_lvlSEXP);
    setDebug(debug_lvl);
    return R_NilValue;
END_RCPP
}
