// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// compareCodonPositions
bool compareCodonPositions(std::string a, std::string b, const char mutNameDelimiter);
RcppExport SEXP _mutscan_compareCodonPositions(SEXP aSEXP, SEXP bSEXP, SEXP mutNameDelimiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type a(aSEXP);
    Rcpp::traits::input_parameter< std::string >::type b(bSEXP);
    Rcpp::traits::input_parameter< const char >::type mutNameDelimiter(mutNameDelimiterSEXP);
    rcpp_result_gen = Rcpp::wrap(compareCodonPositions(a, b, mutNameDelimiter));
    return rcpp_result_gen;
END_RCPP
}
// digestFastqsCpp
List digestFastqsCpp(std::string fastqForward, std::string fastqReverse, bool mergeForwardReverse, bool revComplForward, bool revComplReverse, int skipForward, int skipReverse, int umiLengthForward, int umiLengthReverse, unsigned int constantLengthForward, unsigned int constantLengthReverse, unsigned int variableLengthForward, unsigned int variableLengthReverse, std::string adapterForward, std::string adapterReverse, Rcpp::StringVector wildTypeForward, Rcpp::StringVector wildTypeReverse, std::string constantForward, std::string constantReverse, double avePhredMinForward, double avePhredMinReverse, int variableNMaxForward, int variableNMaxReverse, int umiNMax, unsigned int nbrMutatedCodonsMaxForward, unsigned int nbrMutatedCodonsMaxReverse, CharacterVector forbiddenMutatedCodonsForward, CharacterVector forbiddenMutatedCodonsReverse, double mutatedPhredMinForward, double mutatedPhredMinReverse, std::string mutNameDelimiter, bool verbose);
RcppExport SEXP _mutscan_digestFastqsCpp(SEXP fastqForwardSEXP, SEXP fastqReverseSEXP, SEXP mergeForwardReverseSEXP, SEXP revComplForwardSEXP, SEXP revComplReverseSEXP, SEXP skipForwardSEXP, SEXP skipReverseSEXP, SEXP umiLengthForwardSEXP, SEXP umiLengthReverseSEXP, SEXP constantLengthForwardSEXP, SEXP constantLengthReverseSEXP, SEXP variableLengthForwardSEXP, SEXP variableLengthReverseSEXP, SEXP adapterForwardSEXP, SEXP adapterReverseSEXP, SEXP wildTypeForwardSEXP, SEXP wildTypeReverseSEXP, SEXP constantForwardSEXP, SEXP constantReverseSEXP, SEXP avePhredMinForwardSEXP, SEXP avePhredMinReverseSEXP, SEXP variableNMaxForwardSEXP, SEXP variableNMaxReverseSEXP, SEXP umiNMaxSEXP, SEXP nbrMutatedCodonsMaxForwardSEXP, SEXP nbrMutatedCodonsMaxReverseSEXP, SEXP forbiddenMutatedCodonsForwardSEXP, SEXP forbiddenMutatedCodonsReverseSEXP, SEXP mutatedPhredMinForwardSEXP, SEXP mutatedPhredMinReverseSEXP, SEXP mutNameDelimiterSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fastqForward(fastqForwardSEXP);
    Rcpp::traits::input_parameter< std::string >::type fastqReverse(fastqReverseSEXP);
    Rcpp::traits::input_parameter< bool >::type mergeForwardReverse(mergeForwardReverseSEXP);
    Rcpp::traits::input_parameter< bool >::type revComplForward(revComplForwardSEXP);
    Rcpp::traits::input_parameter< bool >::type revComplReverse(revComplReverseSEXP);
    Rcpp::traits::input_parameter< int >::type skipForward(skipForwardSEXP);
    Rcpp::traits::input_parameter< int >::type skipReverse(skipReverseSEXP);
    Rcpp::traits::input_parameter< int >::type umiLengthForward(umiLengthForwardSEXP);
    Rcpp::traits::input_parameter< int >::type umiLengthReverse(umiLengthReverseSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type constantLengthForward(constantLengthForwardSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type constantLengthReverse(constantLengthReverseSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type variableLengthForward(variableLengthForwardSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type variableLengthReverse(variableLengthReverseSEXP);
    Rcpp::traits::input_parameter< std::string >::type adapterForward(adapterForwardSEXP);
    Rcpp::traits::input_parameter< std::string >::type adapterReverse(adapterReverseSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type wildTypeForward(wildTypeForwardSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type wildTypeReverse(wildTypeReverseSEXP);
    Rcpp::traits::input_parameter< std::string >::type constantForward(constantForwardSEXP);
    Rcpp::traits::input_parameter< std::string >::type constantReverse(constantReverseSEXP);
    Rcpp::traits::input_parameter< double >::type avePhredMinForward(avePhredMinForwardSEXP);
    Rcpp::traits::input_parameter< double >::type avePhredMinReverse(avePhredMinReverseSEXP);
    Rcpp::traits::input_parameter< int >::type variableNMaxForward(variableNMaxForwardSEXP);
    Rcpp::traits::input_parameter< int >::type variableNMaxReverse(variableNMaxReverseSEXP);
    Rcpp::traits::input_parameter< int >::type umiNMax(umiNMaxSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nbrMutatedCodonsMaxForward(nbrMutatedCodonsMaxForwardSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nbrMutatedCodonsMaxReverse(nbrMutatedCodonsMaxReverseSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type forbiddenMutatedCodonsForward(forbiddenMutatedCodonsForwardSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type forbiddenMutatedCodonsReverse(forbiddenMutatedCodonsReverseSEXP);
    Rcpp::traits::input_parameter< double >::type mutatedPhredMinForward(mutatedPhredMinForwardSEXP);
    Rcpp::traits::input_parameter< double >::type mutatedPhredMinReverse(mutatedPhredMinReverseSEXP);
    Rcpp::traits::input_parameter< std::string >::type mutNameDelimiter(mutNameDelimiterSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(digestFastqsCpp(fastqForward, fastqReverse, mergeForwardReverse, revComplForward, revComplReverse, skipForward, skipReverse, umiLengthForward, umiLengthReverse, constantLengthForward, constantLengthReverse, variableLengthForward, variableLengthReverse, adapterForward, adapterReverse, wildTypeForward, wildTypeReverse, constantForward, constantReverse, avePhredMinForward, avePhredMinReverse, variableNMaxForward, variableNMaxReverse, umiNMax, nbrMutatedCodonsMaxForward, nbrMutatedCodonsMaxReverse, forbiddenMutatedCodonsForward, forbiddenMutatedCodonsReverse, mutatedPhredMinForward, mutatedPhredMinReverse, mutNameDelimiter, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mutscan_compareCodonPositions", (DL_FUNC) &_mutscan_compareCodonPositions, 3},
    {"_mutscan_digestFastqsCpp", (DL_FUNC) &_mutscan_digestFastqsCpp, 32},
    {NULL, NULL, 0}
};

RcppExport void R_init_mutscan(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
