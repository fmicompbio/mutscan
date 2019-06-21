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
// test_mergeReadPairPartial
List test_mergeReadPairPartial(std::string seqF, std::vector<int> qualF, std::string seqR, std::vector<int> qualR, size_t minOverlap, size_t maxOverlap, double maxFracMismatchOverlap, bool greedy);
RcppExport SEXP _mutscan_test_mergeReadPairPartial(SEXP seqFSEXP, SEXP qualFSEXP, SEXP seqRSEXP, SEXP qualRSEXP, SEXP minOverlapSEXP, SEXP maxOverlapSEXP, SEXP maxFracMismatchOverlapSEXP, SEXP greedySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type seqF(seqFSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type qualF(qualFSEXP);
    Rcpp::traits::input_parameter< std::string >::type seqR(seqRSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type qualR(qualRSEXP);
    Rcpp::traits::input_parameter< size_t >::type minOverlap(minOverlapSEXP);
    Rcpp::traits::input_parameter< size_t >::type maxOverlap(maxOverlapSEXP);
    Rcpp::traits::input_parameter< double >::type maxFracMismatchOverlap(maxFracMismatchOverlapSEXP);
    Rcpp::traits::input_parameter< bool >::type greedy(greedySEXP);
    rcpp_result_gen = Rcpp::wrap(test_mergeReadPairPartial(seqF, qualF, seqR, qualR, minOverlap, maxOverlap, maxFracMismatchOverlap, greedy));
    return rcpp_result_gen;
END_RCPP
}
// findClosestRefSeq
int findClosestRefSeq(std::string varSeq, Rcpp::StringVector wtSeq);
RcppExport SEXP _mutscan_findClosestRefSeq(SEXP varSeqSEXP, SEXP wtSeqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type varSeq(varSeqSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type wtSeq(wtSeqSEXP);
    rcpp_result_gen = Rcpp::wrap(findClosestRefSeq(varSeq, wtSeq));
    return rcpp_result_gen;
END_RCPP
}
// digestFastqsCpp
List digestFastqsCpp(std::string fastqForward, std::string fastqReverse, bool mergeForwardReverse, size_t minOverlap, size_t maxOverlap, double maxFracMismatchOverlap, bool greedyOverlap, bool revComplForward, bool revComplReverse, int skipForward, int skipReverse, int umiLengthForward, int umiLengthReverse, int constantLengthForward, int constantLengthReverse, int variableLengthForward, int variableLengthReverse, std::string adapterForward, std::string adapterReverse, std::string primerForward, std::string primerReverse, Rcpp::StringVector wildTypeForward, Rcpp::StringVector wildTypeReverse, std::string constantForward, std::string constantReverse, double avePhredMinForward, double avePhredMinReverse, int variableNMaxForward, int variableNMaxReverse, int umiNMax, unsigned int nbrMutatedCodonsMaxForward, unsigned int nbrMutatedCodonsMaxReverse, CharacterVector forbiddenMutatedCodonsForward, CharacterVector forbiddenMutatedCodonsReverse, double mutatedPhredMinForward, double mutatedPhredMinReverse, std::string mutNameDelimiter, int maxNReads, bool verbose);
RcppExport SEXP _mutscan_digestFastqsCpp(SEXP fastqForwardSEXP, SEXP fastqReverseSEXP, SEXP mergeForwardReverseSEXP, SEXP minOverlapSEXP, SEXP maxOverlapSEXP, SEXP maxFracMismatchOverlapSEXP, SEXP greedyOverlapSEXP, SEXP revComplForwardSEXP, SEXP revComplReverseSEXP, SEXP skipForwardSEXP, SEXP skipReverseSEXP, SEXP umiLengthForwardSEXP, SEXP umiLengthReverseSEXP, SEXP constantLengthForwardSEXP, SEXP constantLengthReverseSEXP, SEXP variableLengthForwardSEXP, SEXP variableLengthReverseSEXP, SEXP adapterForwardSEXP, SEXP adapterReverseSEXP, SEXP primerForwardSEXP, SEXP primerReverseSEXP, SEXP wildTypeForwardSEXP, SEXP wildTypeReverseSEXP, SEXP constantForwardSEXP, SEXP constantReverseSEXP, SEXP avePhredMinForwardSEXP, SEXP avePhredMinReverseSEXP, SEXP variableNMaxForwardSEXP, SEXP variableNMaxReverseSEXP, SEXP umiNMaxSEXP, SEXP nbrMutatedCodonsMaxForwardSEXP, SEXP nbrMutatedCodonsMaxReverseSEXP, SEXP forbiddenMutatedCodonsForwardSEXP, SEXP forbiddenMutatedCodonsReverseSEXP, SEXP mutatedPhredMinForwardSEXP, SEXP mutatedPhredMinReverseSEXP, SEXP mutNameDelimiterSEXP, SEXP maxNReadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type fastqForward(fastqForwardSEXP);
    Rcpp::traits::input_parameter< std::string >::type fastqReverse(fastqReverseSEXP);
    Rcpp::traits::input_parameter< bool >::type mergeForwardReverse(mergeForwardReverseSEXP);
    Rcpp::traits::input_parameter< size_t >::type minOverlap(minOverlapSEXP);
    Rcpp::traits::input_parameter< size_t >::type maxOverlap(maxOverlapSEXP);
    Rcpp::traits::input_parameter< double >::type maxFracMismatchOverlap(maxFracMismatchOverlapSEXP);
    Rcpp::traits::input_parameter< bool >::type greedyOverlap(greedyOverlapSEXP);
    Rcpp::traits::input_parameter< bool >::type revComplForward(revComplForwardSEXP);
    Rcpp::traits::input_parameter< bool >::type revComplReverse(revComplReverseSEXP);
    Rcpp::traits::input_parameter< int >::type skipForward(skipForwardSEXP);
    Rcpp::traits::input_parameter< int >::type skipReverse(skipReverseSEXP);
    Rcpp::traits::input_parameter< int >::type umiLengthForward(umiLengthForwardSEXP);
    Rcpp::traits::input_parameter< int >::type umiLengthReverse(umiLengthReverseSEXP);
    Rcpp::traits::input_parameter< int >::type constantLengthForward(constantLengthForwardSEXP);
    Rcpp::traits::input_parameter< int >::type constantLengthReverse(constantLengthReverseSEXP);
    Rcpp::traits::input_parameter< int >::type variableLengthForward(variableLengthForwardSEXP);
    Rcpp::traits::input_parameter< int >::type variableLengthReverse(variableLengthReverseSEXP);
    Rcpp::traits::input_parameter< std::string >::type adapterForward(adapterForwardSEXP);
    Rcpp::traits::input_parameter< std::string >::type adapterReverse(adapterReverseSEXP);
    Rcpp::traits::input_parameter< std::string >::type primerForward(primerForwardSEXP);
    Rcpp::traits::input_parameter< std::string >::type primerReverse(primerReverseSEXP);
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
    Rcpp::traits::input_parameter< int >::type maxNReads(maxNReadsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(digestFastqsCpp(fastqForward, fastqReverse, mergeForwardReverse, minOverlap, maxOverlap, maxFracMismatchOverlap, greedyOverlap, revComplForward, revComplReverse, skipForward, skipReverse, umiLengthForward, umiLengthReverse, constantLengthForward, constantLengthReverse, variableLengthForward, variableLengthReverse, adapterForward, adapterReverse, primerForward, primerReverse, wildTypeForward, wildTypeReverse, constantForward, constantReverse, avePhredMinForward, avePhredMinReverse, variableNMaxForward, variableNMaxReverse, umiNMax, nbrMutatedCodonsMaxForward, nbrMutatedCodonsMaxReverse, forbiddenMutatedCodonsForward, forbiddenMutatedCodonsReverse, mutatedPhredMinForward, mutatedPhredMinReverse, mutNameDelimiter, maxNReads, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mutscan_compareCodonPositions", (DL_FUNC) &_mutscan_compareCodonPositions, 3},
    {"_mutscan_test_mergeReadPairPartial", (DL_FUNC) &_mutscan_test_mergeReadPairPartial, 8},
    {"_mutscan_findClosestRefSeq", (DL_FUNC) &_mutscan_findClosestRefSeq, 2},
    {"_mutscan_digestFastqsCpp", (DL_FUNC) &_mutscan_digestFastqsCpp, 39},
    {NULL, NULL, 0}
};

RcppExport void R_init_mutscan(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
