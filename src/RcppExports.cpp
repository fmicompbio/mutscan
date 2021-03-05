// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// levenshtein_distance
int levenshtein_distance(const std::string& str1, const std::string& str2);
RcppExport SEXP _mutscan_levenshtein_distance(SEXP str1SEXP, SEXP str2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type str1(str1SEXP);
    Rcpp::traits::input_parameter< const std::string& >::type str2(str2SEXP);
    rcpp_result_gen = Rcpp::wrap(levenshtein_distance(str1, str2));
    return rcpp_result_gen;
END_RCPP
}
// hamming_distance
int hamming_distance(const std::string& str1, const std::string& str2);
RcppExport SEXP _mutscan_hamming_distance(SEXP str1SEXP, SEXP str2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type str1(str1SEXP);
    Rcpp::traits::input_parameter< const std::string& >::type str2(str2SEXP);
    rcpp_result_gen = Rcpp::wrap(hamming_distance(str1, str2));
    return rcpp_result_gen;
END_RCPP
}
// hamming_shift_distance
int hamming_shift_distance(const std::string& str1, const std::string& str2, int max_abs_shift);
RcppExport SEXP _mutscan_hamming_shift_distance(SEXP str1SEXP, SEXP str2SEXP, SEXP max_abs_shiftSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type str1(str1SEXP);
    Rcpp::traits::input_parameter< const std::string& >::type str2(str2SEXP);
    Rcpp::traits::input_parameter< int >::type max_abs_shift(max_abs_shiftSEXP);
    rcpp_result_gen = Rcpp::wrap(hamming_shift_distance(str1, str2, max_abs_shift));
    return rcpp_result_gen;
END_RCPP
}
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
// test_decomposeRead
List test_decomposeRead(const std::string sseq, const std::string squal, const std::string elements, const std::vector<int> elementLengths, const std::vector<std::string> primerSeqs, std::string umiSeq, std::string varSeq, std::string varQual, std::string constSeq, std::string constQual, int nNoPrimer, int nReadWrongLength);
RcppExport SEXP _mutscan_test_decomposeRead(SEXP sseqSEXP, SEXP squalSEXP, SEXP elementsSEXP, SEXP elementLengthsSEXP, SEXP primerSeqsSEXP, SEXP umiSeqSEXP, SEXP varSeqSEXP, SEXP varQualSEXP, SEXP constSeqSEXP, SEXP constQualSEXP, SEXP nNoPrimerSEXP, SEXP nReadWrongLengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type sseq(sseqSEXP);
    Rcpp::traits::input_parameter< const std::string >::type squal(squalSEXP);
    Rcpp::traits::input_parameter< const std::string >::type elements(elementsSEXP);
    Rcpp::traits::input_parameter< const std::vector<int> >::type elementLengths(elementLengthsSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string> >::type primerSeqs(primerSeqsSEXP);
    Rcpp::traits::input_parameter< std::string >::type umiSeq(umiSeqSEXP);
    Rcpp::traits::input_parameter< std::string >::type varSeq(varSeqSEXP);
    Rcpp::traits::input_parameter< std::string >::type varQual(varQualSEXP);
    Rcpp::traits::input_parameter< std::string >::type constSeq(constSeqSEXP);
    Rcpp::traits::input_parameter< std::string >::type constQual(constQualSEXP);
    Rcpp::traits::input_parameter< int >::type nNoPrimer(nNoPrimerSEXP);
    Rcpp::traits::input_parameter< int >::type nReadWrongLength(nReadWrongLengthSEXP);
    rcpp_result_gen = Rcpp::wrap(test_decomposeRead(sseq, squal, elements, elementLengths, primerSeqs, umiSeq, varSeq, varQual, constSeq, constQual, nNoPrimer, nReadWrongLength));
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
int findClosestRefSeq(std::string& varSeq, std::vector<std::string>& wtSeq, size_t upperBoundMismatch, int& sim);
RcppExport SEXP _mutscan_findClosestRefSeq(SEXP varSeqSEXP, SEXP wtSeqSEXP, SEXP upperBoundMismatchSEXP, SEXP simSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type varSeq(varSeqSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string>& >::type wtSeq(wtSeqSEXP);
    Rcpp::traits::input_parameter< size_t >::type upperBoundMismatch(upperBoundMismatchSEXP);
    Rcpp::traits::input_parameter< int& >::type sim(simSEXP);
    rcpp_result_gen = Rcpp::wrap(findClosestRefSeq(varSeq, wtSeq, upperBoundMismatch, sim));
    return rcpp_result_gen;
END_RCPP
}
// findClosestRefSeqEarlyStop
int findClosestRefSeqEarlyStop(std::string& varSeq, std::vector<std::string>& wtSeq, size_t upperBoundMismatch, int& sim);
RcppExport SEXP _mutscan_findClosestRefSeqEarlyStop(SEXP varSeqSEXP, SEXP wtSeqSEXP, SEXP upperBoundMismatchSEXP, SEXP simSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string& >::type varSeq(varSeqSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string>& >::type wtSeq(wtSeqSEXP);
    Rcpp::traits::input_parameter< size_t >::type upperBoundMismatch(upperBoundMismatchSEXP);
    Rcpp::traits::input_parameter< int& >::type sim(simSEXP);
    rcpp_result_gen = Rcpp::wrap(findClosestRefSeqEarlyStop(varSeq, wtSeq, upperBoundMismatch, sim));
    return rcpp_result_gen;
END_RCPP
}
// digestFastqsCpp
List digestFastqsCpp(std::vector<std::string> fastqForwardVect, std::vector<std::string> fastqReverseVect, bool mergeForwardReverse, size_t minOverlap, size_t maxOverlap, double maxFracMismatchOverlap, bool greedyOverlap, bool revComplForward, bool revComplReverse, std::string elementsForward, std::vector<int> elementLengthsForward, std::string elementsReverse, std::vector<int> elementLengthsReverse, std::string adapterForward, std::string adapterReverse, std::vector<std::string> primerForward, std::vector<std::string> primerReverse, std::vector<std::string> wildTypeForward, std::vector<std::string> wildTypeForwardNames, std::vector<std::string> wildTypeReverse, std::vector<std::string> wildTypeReverseNames, std::vector<std::string> constantForward, std::vector<std::string> constantReverse, double avePhredMinForward, double avePhredMinReverse, int variableNMaxForward, int variableNMaxReverse, int umiNMax, int nbrMutatedCodonsMaxForward, int nbrMutatedCodonsMaxReverse, int nbrMutatedBasesMaxForward, int nbrMutatedBasesMaxReverse, CharacterVector forbiddenMutatedCodonsForward, CharacterVector forbiddenMutatedCodonsReverse, bool useTreeWTmatch, double mutatedPhredMinForward, double mutatedPhredMinReverse, std::string mutNameDelimiter, int constantMaxDistForward, int constantMaxDistReverse, double variableCollapseMaxDist, int variableCollapseMinReads, double variableCollapseMinRatio, double umiCollapseMaxDist, std::string filteredReadsFastqForward, std::string filteredReadsFastqReverse, int maxNReads, bool verbose, int nThreads, int chunkSize);
RcppExport SEXP _mutscan_digestFastqsCpp(SEXP fastqForwardVectSEXP, SEXP fastqReverseVectSEXP, SEXP mergeForwardReverseSEXP, SEXP minOverlapSEXP, SEXP maxOverlapSEXP, SEXP maxFracMismatchOverlapSEXP, SEXP greedyOverlapSEXP, SEXP revComplForwardSEXP, SEXP revComplReverseSEXP, SEXP elementsForwardSEXP, SEXP elementLengthsForwardSEXP, SEXP elementsReverseSEXP, SEXP elementLengthsReverseSEXP, SEXP adapterForwardSEXP, SEXP adapterReverseSEXP, SEXP primerForwardSEXP, SEXP primerReverseSEXP, SEXP wildTypeForwardSEXP, SEXP wildTypeForwardNamesSEXP, SEXP wildTypeReverseSEXP, SEXP wildTypeReverseNamesSEXP, SEXP constantForwardSEXP, SEXP constantReverseSEXP, SEXP avePhredMinForwardSEXP, SEXP avePhredMinReverseSEXP, SEXP variableNMaxForwardSEXP, SEXP variableNMaxReverseSEXP, SEXP umiNMaxSEXP, SEXP nbrMutatedCodonsMaxForwardSEXP, SEXP nbrMutatedCodonsMaxReverseSEXP, SEXP nbrMutatedBasesMaxForwardSEXP, SEXP nbrMutatedBasesMaxReverseSEXP, SEXP forbiddenMutatedCodonsForwardSEXP, SEXP forbiddenMutatedCodonsReverseSEXP, SEXP useTreeWTmatchSEXP, SEXP mutatedPhredMinForwardSEXP, SEXP mutatedPhredMinReverseSEXP, SEXP mutNameDelimiterSEXP, SEXP constantMaxDistForwardSEXP, SEXP constantMaxDistReverseSEXP, SEXP variableCollapseMaxDistSEXP, SEXP variableCollapseMinReadsSEXP, SEXP variableCollapseMinRatioSEXP, SEXP umiCollapseMaxDistSEXP, SEXP filteredReadsFastqForwardSEXP, SEXP filteredReadsFastqReverseSEXP, SEXP maxNReadsSEXP, SEXP verboseSEXP, SEXP nThreadsSEXP, SEXP chunkSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type fastqForwardVect(fastqForwardVectSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type fastqReverseVect(fastqReverseVectSEXP);
    Rcpp::traits::input_parameter< bool >::type mergeForwardReverse(mergeForwardReverseSEXP);
    Rcpp::traits::input_parameter< size_t >::type minOverlap(minOverlapSEXP);
    Rcpp::traits::input_parameter< size_t >::type maxOverlap(maxOverlapSEXP);
    Rcpp::traits::input_parameter< double >::type maxFracMismatchOverlap(maxFracMismatchOverlapSEXP);
    Rcpp::traits::input_parameter< bool >::type greedyOverlap(greedyOverlapSEXP);
    Rcpp::traits::input_parameter< bool >::type revComplForward(revComplForwardSEXP);
    Rcpp::traits::input_parameter< bool >::type revComplReverse(revComplReverseSEXP);
    Rcpp::traits::input_parameter< std::string >::type elementsForward(elementsForwardSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type elementLengthsForward(elementLengthsForwardSEXP);
    Rcpp::traits::input_parameter< std::string >::type elementsReverse(elementsReverseSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type elementLengthsReverse(elementLengthsReverseSEXP);
    Rcpp::traits::input_parameter< std::string >::type adapterForward(adapterForwardSEXP);
    Rcpp::traits::input_parameter< std::string >::type adapterReverse(adapterReverseSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type primerForward(primerForwardSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type primerReverse(primerReverseSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type wildTypeForward(wildTypeForwardSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type wildTypeForwardNames(wildTypeForwardNamesSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type wildTypeReverse(wildTypeReverseSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type wildTypeReverseNames(wildTypeReverseNamesSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type constantForward(constantForwardSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type constantReverse(constantReverseSEXP);
    Rcpp::traits::input_parameter< double >::type avePhredMinForward(avePhredMinForwardSEXP);
    Rcpp::traits::input_parameter< double >::type avePhredMinReverse(avePhredMinReverseSEXP);
    Rcpp::traits::input_parameter< int >::type variableNMaxForward(variableNMaxForwardSEXP);
    Rcpp::traits::input_parameter< int >::type variableNMaxReverse(variableNMaxReverseSEXP);
    Rcpp::traits::input_parameter< int >::type umiNMax(umiNMaxSEXP);
    Rcpp::traits::input_parameter< int >::type nbrMutatedCodonsMaxForward(nbrMutatedCodonsMaxForwardSEXP);
    Rcpp::traits::input_parameter< int >::type nbrMutatedCodonsMaxReverse(nbrMutatedCodonsMaxReverseSEXP);
    Rcpp::traits::input_parameter< int >::type nbrMutatedBasesMaxForward(nbrMutatedBasesMaxForwardSEXP);
    Rcpp::traits::input_parameter< int >::type nbrMutatedBasesMaxReverse(nbrMutatedBasesMaxReverseSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type forbiddenMutatedCodonsForward(forbiddenMutatedCodonsForwardSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type forbiddenMutatedCodonsReverse(forbiddenMutatedCodonsReverseSEXP);
    Rcpp::traits::input_parameter< bool >::type useTreeWTmatch(useTreeWTmatchSEXP);
    Rcpp::traits::input_parameter< double >::type mutatedPhredMinForward(mutatedPhredMinForwardSEXP);
    Rcpp::traits::input_parameter< double >::type mutatedPhredMinReverse(mutatedPhredMinReverseSEXP);
    Rcpp::traits::input_parameter< std::string >::type mutNameDelimiter(mutNameDelimiterSEXP);
    Rcpp::traits::input_parameter< int >::type constantMaxDistForward(constantMaxDistForwardSEXP);
    Rcpp::traits::input_parameter< int >::type constantMaxDistReverse(constantMaxDistReverseSEXP);
    Rcpp::traits::input_parameter< double >::type variableCollapseMaxDist(variableCollapseMaxDistSEXP);
    Rcpp::traits::input_parameter< int >::type variableCollapseMinReads(variableCollapseMinReadsSEXP);
    Rcpp::traits::input_parameter< double >::type variableCollapseMinRatio(variableCollapseMinRatioSEXP);
    Rcpp::traits::input_parameter< double >::type umiCollapseMaxDist(umiCollapseMaxDistSEXP);
    Rcpp::traits::input_parameter< std::string >::type filteredReadsFastqForward(filteredReadsFastqForwardSEXP);
    Rcpp::traits::input_parameter< std::string >::type filteredReadsFastqReverse(filteredReadsFastqReverseSEXP);
    Rcpp::traits::input_parameter< int >::type maxNReads(maxNReadsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type nThreads(nThreadsSEXP);
    Rcpp::traits::input_parameter< int >::type chunkSize(chunkSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(digestFastqsCpp(fastqForwardVect, fastqReverseVect, mergeForwardReverse, minOverlap, maxOverlap, maxFracMismatchOverlap, greedyOverlap, revComplForward, revComplReverse, elementsForward, elementLengthsForward, elementsReverse, elementLengthsReverse, adapterForward, adapterReverse, primerForward, primerReverse, wildTypeForward, wildTypeForwardNames, wildTypeReverse, wildTypeReverseNames, constantForward, constantReverse, avePhredMinForward, avePhredMinReverse, variableNMaxForward, variableNMaxReverse, umiNMax, nbrMutatedCodonsMaxForward, nbrMutatedCodonsMaxReverse, nbrMutatedBasesMaxForward, nbrMutatedBasesMaxReverse, forbiddenMutatedCodonsForward, forbiddenMutatedCodonsReverse, useTreeWTmatch, mutatedPhredMinForward, mutatedPhredMinReverse, mutNameDelimiter, constantMaxDistForward, constantMaxDistReverse, variableCollapseMaxDist, variableCollapseMinReads, variableCollapseMinRatio, umiCollapseMaxDist, filteredReadsFastqForward, filteredReadsFastqReverse, maxNReads, verbose, nThreads, chunkSize));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_mod_BKtree();

static const R_CallMethodDef CallEntries[] = {
    {"_mutscan_levenshtein_distance", (DL_FUNC) &_mutscan_levenshtein_distance, 2},
    {"_mutscan_hamming_distance", (DL_FUNC) &_mutscan_hamming_distance, 2},
    {"_mutscan_hamming_shift_distance", (DL_FUNC) &_mutscan_hamming_shift_distance, 3},
    {"_mutscan_compareCodonPositions", (DL_FUNC) &_mutscan_compareCodonPositions, 3},
    {"_mutscan_test_decomposeRead", (DL_FUNC) &_mutscan_test_decomposeRead, 12},
    {"_mutscan_test_mergeReadPairPartial", (DL_FUNC) &_mutscan_test_mergeReadPairPartial, 8},
    {"_mutscan_findClosestRefSeq", (DL_FUNC) &_mutscan_findClosestRefSeq, 4},
    {"_mutscan_findClosestRefSeqEarlyStop", (DL_FUNC) &_mutscan_findClosestRefSeqEarlyStop, 4},
    {"_mutscan_digestFastqsCpp", (DL_FUNC) &_mutscan_digestFastqsCpp, 50},
    {"_rcpp_module_boot_mod_BKtree", (DL_FUNC) &_rcpp_module_boot_mod_BKtree, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_mutscan(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
