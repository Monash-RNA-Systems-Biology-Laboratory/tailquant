// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// quality_clip
double quality_clip(std::string seq, std::string qual, char c, double penalty);
RcppExport SEXP _tailquant_quality_clip(SEXP seqSEXP, SEXP qualSEXP, SEXP cSEXP, SEXP penaltySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< std::string >::type qual(qualSEXP);
    Rcpp::traits::input_parameter< char >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type penalty(penaltySEXP);
    rcpp_result_gen = Rcpp::wrap(quality_clip(seq, qual, c, penalty));
    return rcpp_result_gen;
END_RCPP
}
// scan
NumericVector scan(std::string s, char c, int n, double penalty);
RcppExport SEXP _tailquant_scan(SEXP sSEXP, SEXP cSEXP, SEXP nSEXP, SEXP penaltySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type s(sSEXP);
    Rcpp::traits::input_parameter< char >::type c(cSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type penalty(penaltySEXP);
    rcpp_result_gen = Rcpp::wrap(scan(s, c, n, penalty));
    return rcpp_result_gen;
END_RCPP
}
// scan_suffix
NumericVector scan_suffix(std::string s, char c, std::string suffix, int n, double penalty, double suffix_penalty);
RcppExport SEXP _tailquant_scan_suffix(SEXP sSEXP, SEXP cSEXP, SEXP suffixSEXP, SEXP nSEXP, SEXP penaltySEXP, SEXP suffix_penaltySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type s(sSEXP);
    Rcpp::traits::input_parameter< char >::type c(cSEXP);
    Rcpp::traits::input_parameter< std::string >::type suffix(suffixSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type penalty(penaltySEXP);
    Rcpp::traits::input_parameter< double >::type suffix_penalty(suffix_penaltySEXP);
    rcpp_result_gen = Rcpp::wrap(scan_suffix(s, c, suffix, n, penalty, suffix_penalty));
    return rcpp_result_gen;
END_RCPP
}
// scan_from
NumericVector scan_from(std::string s, char c, int start, int n, double penalty);
RcppExport SEXP _tailquant_scan_from(SEXP sSEXP, SEXP cSEXP, SEXP startSEXP, SEXP nSEXP, SEXP penaltySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type s(sSEXP);
    Rcpp::traits::input_parameter< char >::type c(cSEXP);
    Rcpp::traits::input_parameter< int >::type start(startSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type penalty(penaltySEXP);
    rcpp_result_gen = Rcpp::wrap(scan_from(s, c, start, n, penalty));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tailquant_quality_clip", (DL_FUNC) &_tailquant_quality_clip, 4},
    {"_tailquant_scan", (DL_FUNC) &_tailquant_scan, 4},
    {"_tailquant_scan_suffix", (DL_FUNC) &_tailquant_scan_suffix, 6},
    {"_tailquant_scan_from", (DL_FUNC) &_tailquant_scan_from, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_tailquant(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
