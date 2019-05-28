#include <Rcpp.h>
#include <cerrno>
#include <zlib.h>
#include <string>
#include <vector>
#include <numeric>
#include <map>
#include <algorithm>    // std::sort

#define BUFFER_SIZE 4096

using namespace Rcpp;

// check if the last gzgets reached end of file
// return either:
// - true (reached end of file)
// - false (not yet reached end of file, but reading was ok)
// - fall back to R with an error (reading was not ok)
bool reached_end_of_file(gzFile file, char *ret) {
  static int err;
  if (ret == Z_NULL) {
    if (gzeof(file)) {
      return true;
    } else {
      const char * error_string;
      error_string = gzerror (file, &err);
      if (err) {
        stop(error_string);
      }
    }
  } else {
    return false;
  }
}

// read next four lines from gzipped file and store
// second and forth in *seq and *qual
// return either:
// - true (reached end of file, I am done)
// - false (not yet reached end of file, not done yet)
// - nothing (encounetered an error, fall back to R from reached_end_of_file())
bool get_next_seq(gzFile file, char *seq, char *qual) {
  // sequence identifier
  if (reached_end_of_file(file, gzgets(file, seq, BUFFER_SIZE))) {
    return true;
  }
  // sequence
  if (reached_end_of_file(file, gzgets(file, seq, BUFFER_SIZE))) {
    return true;
  }
  // quality identifier
  if (reached_end_of_file(file, gzgets(file, qual, BUFFER_SIZE))) {
    return true;
  }
  // quality
  if (reached_end_of_file(file, gzgets(file, qual, BUFFER_SIZE))) {
    return true;
  }
  
  return false;
}

// create the complement of a base
char complement(char n) {   
  switch(n) {   
    case 'A':
    case 'a':
      return 'T';
    case 'T':
    case 't':
      return 'A';
    case 'G':
    case 'g':
      return 'C';
    case 'C':
    case 'c':
      return 'G';
    case 'N':
    case 'n':
      return 'N';
  }   
  stop("Non-base character in sequence - aborting");
}   

// compare position of codons of the form (std::string)"x123NNN_"
// [[Rcpp::export]]
bool compareCodonPositions(std::string a, std::string b) {
  int posa = std::stoi(a.substr(1, a.length() - 5));
  int posb = std::stoi(b.substr(1, b.length() - 5));
  return (posa < posb);
}

// compare read to wiltype sequence,
// identify mutated bases/codons, filter, update counters
// and add to the name
// returns true if the read needs to be filtered out (and a counter has been incremented)
bool compareToWildtype(const std::string varSeq, const std::string wtSeq,
                       const std::vector<int> varIntQual, const double mutatedPhredMin,
                       const int nbrMutatedCodonsMax, const std::set<std::string> &forbiddenCodons,
                       const std::string codonPrefix, int &nMutQualTooLow, int &nTooManyMutCodons,
                       int &nForbiddenCodons, std::string &mutantName) {
  static std::set<std::string> mutatedCodons;
  static std::set<std::string>::iterator mutatedCodonIt;
  bool hasLowQualMutation, hasForbidden;
  
  // filter if there are too many mutated codons
  mutatedCodons.clear();
  hasLowQualMutation = false;
  for (size_t i = 0; i < varSeq.length(); i++) {
    if (varSeq[i] != wtSeq[i]) { // found mismatching base
      // record if the mutated base quality is below a threshold
      if (varIntQual[i] < mutatedPhredMin) {
        hasLowQualMutation = true;
        break;
      }
      // add codon to mutatedCodons
      mutatedCodons.insert(codonPrefix + std::to_string((int)(i / 3) * 3) +
        varSeq.substr((int)(i / 3) * 3, 3) +
        std::string("_"));
    }
  }
  if (hasLowQualMutation) {
    nMutQualTooLow++;
    return true;
  }
  // check if there are too many mutated codons
  if (mutatedCodons.size() > nbrMutatedCodonsMax) {
    nTooManyMutCodons++;
    return true;
  }
  // check if there are forbidden codons
  hasForbidden = false;
  for ( mutatedCodonIt = mutatedCodons.begin(); mutatedCodonIt != mutatedCodons.end(); mutatedCodonIt++) {
    if (forbiddenCodons.find((*mutatedCodonIt).substr((*mutatedCodonIt).length()-4, 3)) != forbiddenCodons.end()) { // found forbidden codon
      hasForbidden = true;
      break;
    }
  }
  if (hasForbidden) {
    nForbiddenCodons++;
    return true;
  }
  // create name for mutant
  std::vector<std::string> mutatedCodonsSorted(mutatedCodons.begin(), mutatedCodons.end());
  std::sort (mutatedCodonsSorted.begin(), mutatedCodonsSorted.end(), compareCodonPositions);
  for (size_t i = 0; i < mutatedCodonsSorted.size(); i++) {
    mutantName += mutatedCodonsSorted[i];
  }
  
  return false;
}

// compare constSeq to constant sequence and update
// counter by match/mismatch and by base Phred quality
// return true
bool tabulateBasesByQual(const std::string constSeq, const std::string constant,
                         const std::vector<int> constIntQual,
                         std::vector<int> &nPhredCorrect, std::vector<int> &nPhredMismatch) {
  for (size_t i = 0; i < constSeq.length(); i++) { // for each base
    if (constSeq[i] != constant[i]) { // found mismatch
      nPhredMismatch[constIntQual[i]]++;
    } else {                          // found match
      nPhredCorrect[constIntQual[i]]++;
    }
  }
  return true;
}


// stores information about retained mutants
struct mutantInfo {
  std::set<std::string> umi;   // set of unique UMIs observed for read pairs of that mutant
  int nReads;                  // number of reads
  std::string mutantName;      // name of the mutant
};

//' @title Read, filter and digest sequences from two fastq files.
//'
//' @description
//' \code{readAndDigestFastqs} reads sequences for a pair of fastq files
//' and digests them (extracts umis, constant and variable parts, filters,
//' extracts mismatch information from constant and counts the observed
//' unique variable parts).
//'
//' @details
//' more details.
//'
//' @param experimentType character(1), either "cis" or "trans". If this is set
//'   to "cis", the variable sequences from the forward and reverse reads will be
//'   consolidated into one single sequence.
//' @param fastqForward,fastqReverse character(1), paths to FASTQ files
//'   corresponding to forward and reverse reads, respectively.
//' @param skipForward,skipReverse numeric(1), the number of bases to skip in the
//'   start of each forward and reverse read, respectively.
//' @param umiLengthForward,umiLengthReverse numeric(1), the length of the
//'   barcode (UMI) sequence in the forward/reverse reads, respectively, not
//'   including the skipped bases (defined by
//'   \code{skipForward}/\code{skipReverse}).
//' @param constantLengthForward,constantLengthReverse numeric(1), the length of
//'   the constant sequence in the forward/reverse reads, respectively.
//' @param variableLengthForward,variableLengthReverse numeric(1), the length of
//'   the variable sequence in the forward/reverse reads, respectively.
//' @param adapterForward,adapterReverse character(1), the adapter sequence for
//'   forward/reverse reads, respectively. If a forward/reverse read contains the
//'   corresponding adapter sequence, the sequence pair will be filtered out.
//'   If set to \code{NULL}, no adapter filtering is performed. The number of
//'   filtered read pairs are reported in the return value.
//' @param wildTypeForward,wildTypeReverse character(1), the wild type sequence
//'   for the forward and reverse variable region.
//' @param constantForward,constantReverse character(1), the expected constant
//'   forward and reverse sequences.
//' @param avePhredMin numeric(1) Minimum average Phred score in the variable
//'   region for a read to be retained. If L contains both forward and reverse
//'   variable regions, the minimum average Phred score has to be achieved in
//'   both for a read pair to be retained.
//' @param variableNMax numeric(1) Maximum number of Ns allowed in the variable
//'   region for a read to be retained.
//' @param umiNMax numeric(1) Maximum number of Ns allowed in the UMI for a read
//'   to be retained.
//' @param nbrMutatedCodonsMax numeric(1) Maximum number of mutated codons that
//'   are allowed.
//' @param forbiddenMutatedCodons character vector. Codons (can contain ambiguous
//'   IUPAC characters, see \code{\link[Biostrings]{IUPAC_CODE_MAP}}). If a read
//'   pair contains a mutated codon matching this pattern, it will be filtered
//'   out.
//' @param mutatedPhredMin numeric(1) Minimum Phred score of a mutated base for the
//'   read to be retained. If any mutated base has a Phred score lower than
//'   \code{mutatedPhredMin}, the read will be discarded.
//' @param verbose logical(1), whether to print out progress messages.
//'
//' @return A list with ---more details---.
//'
//' @export
// [[Rcpp::export]]
List digestFastqs(std::string experimentType,
                  std::string fastqForward, std::string fastqReverse,
                  int skipForward = 1, int skipReverse = 1,
                  int umiLengthForward = 10, int umiLengthReverse = 8,
                  int constantLengthForward = 18, int constantLengthReverse = 20,
                  int variableLengthForward = 96, int variableLengthReverse = 96,
                  std::string adapterForward = "", std::string adapterReverse = "",
                  std::string wildTypeForward = "", std::string wildTypeReverse = "", 
                  std::string constantForward = "", std::string constantReverse = "", 
                  double avePhredMin = 20.0, int variableNMax = 0, int umiNMax = 0,
                  int nbrMutatedCodonsMax = 1,
                  CharacterVector forbiddenMutatedCodons = "NNW",
                  double mutatedPhredMin = 0.0,
                  bool verbose = false) {

  // --------------------------------------------------------------------------
  // pre-flight checks
  // --------------------------------------------------------------------------
  if (experimentType.compare("cis") != 0 && experimentType.compare("trans") != 0) {
    stop("experimentType' must be either 'cis' or 'trans'");
  }
  if (wildTypeForward.compare("") == 0) {
    Rcout << "skipping number-of-mutated-codons filter (missing 'wildTypeForward')";
  }
  // need more here

  // Biostrings::IUPAC_CODE_MAP
  std::map<char,std::vector<char>> IUPAC;
  IUPAC['A'] = std::vector<char>({'A'});
  IUPAC['C'] = std::vector<char>({'C'});
  IUPAC['G'] = std::vector<char>({'G'});
  IUPAC['T'] = std::vector<char>({'T'});
  IUPAC['M'] = std::vector<char>({'A','C'});
  IUPAC['R'] = std::vector<char>({'A','G'});
  IUPAC['W'] = std::vector<char>({'A','T'});
  IUPAC['S'] = std::vector<char>({'C','G'});
  IUPAC['Y'] = std::vector<char>({'C','T'});
  IUPAC['K'] = std::vector<char>({'G','T'});
  IUPAC['V'] = std::vector<char>({'A','C','G'});
  IUPAC['H'] = std::vector<char>({'A','C','T'});
  IUPAC['D'] = std::vector<char>({'A','G','T'});
  IUPAC['B'] = std::vector<char>({'C','G','T'});
  IUPAC['N'] = std::vector<char>({'A','C','G','T'});
  
  
  // --------------------------------------------------------------------------
  // open fastq files
  // --------------------------------------------------------------------------
  gzFile file1, file2;
  file1 = gzopen(fastqForward.c_str(), "rb");   
  if( !file1 ) {
    if( errno ) {
      stop("Failed to open file '", fastqForward,  "': ",
           strerror(errno), " (errno=", errno, ")");
      
    } else {
      stop("Failed to open file '", fastqForward, "': zlib out of memory");
    }
  }	    
  file2 = gzopen(fastqReverse.c_str(), "rb");   
  if( !file1 ) {
    if( errno ) {
      stop("Failed to open file '", fastqReverse,  "': ",
           strerror(errno), " (errno=", errno, ")");
      
    } else {
      stop("Failed to open file '", fastqReverse, "': zlib out of memory");
    }
  }	    

  
  // --------------------------------------------------------------------------
  // digest reads one by one
  // --------------------------------------------------------------------------
  char seq1[BUFFER_SIZE];
  char qual1[BUFFER_SIZE];
  char seq2[BUFFER_SIZE];
  char qual2[BUFFER_SIZE];
  bool done = false;
  int nTot = 0, nAdapter = 0, nAvgVarQualTooLow = 0, nTooManyNinVar = 0, nTooManyNinUMI = 0;
  int nTooManyMutCodons = 0, nForbiddenCodons = 0, nMutQualTooLow = 0, nRetain = 0;
  std::string varSeqForward, varSeqReverse, varQualForward, varQualReverse, umiSeq;
  std::string constSeqForward, constSeqReverse, constQualForward, constQualReverse;
  std::vector<int> varIntQualForward(variableLengthForward,0);
  std::vector<int> varIntQualReverse(variableLengthReverse,0);
  std::vector<int> constIntQualForward(constantLengthForward,0);
  std::vector<int> constIntQualReverse(constantLengthReverse,0);
  std::set<std::string> forbiddenCodons;
  std::string mutantName;
  std::map<std::string, mutantInfo> mutantSummary;
  std::map<std::string, mutantInfo>::iterator mutantSummaryIt;
  std::vector<int> nPhredCorrectForward(100, 0), nPhredMismatchForward(100, 0);
  std::vector<int> nPhredCorrectReverse(100, 0), nPhredMismatchReverse(100, 0);
  
  // enumerate forbidden codons based on forbiddenMutatedCodons
  if (verbose) {
    Rcout << "start enumerating forbidden codons" << std::endl;
  }
  std::string codon("NNN");
  for (size_t i = 0; i < forbiddenMutatedCodons.length(); i++) {
    std::vector<char> B1s = IUPAC[forbiddenMutatedCodons[i][0]];
    std::vector<char> B2s = IUPAC[forbiddenMutatedCodons[i][1]];
    std::vector<char> B3s = IUPAC[forbiddenMutatedCodons[i][2]];
    for (std::vector<char>::iterator b1 = B1s.begin(); b1 != B1s.end(); b1++) {
      codon[0] = (*b1);
      for (std::vector<char>::iterator b2 = B2s.begin(); b2 != B2s.end(); b2++) {
        codon[1] = (*b2);
        for (std::vector<char>::iterator b3 = B3s.begin(); b3 != B3s.end(); b3++) {
          codon[2] = (*b3);
          forbiddenCodons.insert(codon);
        }
      }
    }
  }
  if (verbose) {
    Rcout << "done enumerating forbidden codons (" << forbiddenCodons.size() << ")" << std::endl;
  }
  
  // iterate over sequences
  while (done == false) {
    mutantName = ""; // start with empty name
    
    // read sequence pair
    done = get_next_seq(file1, seq1, qual1);
    done = (done || get_next_seq(file2, seq2, qual2));
    
    // check if end-of-file was reached
    if (done) {
      break;
    }
    
    // update counters
    nTot++;
    
    // convert C char* to C++ string
    std::string sseq1(seq1);
    std::string sseq2(seq2);
    std::string squal1(qual1);
    std::string squal2(qual2);
    
    // search for adapter sequences and filter read pairs
    if ((adapterForward.compare("") != 0 && sseq1.find(adapterForward) != std::string::npos) ||
        (adapterReverse.compare("") != 0 && sseq2.find(adapterReverse) != std::string::npos)) {
      nAdapter++;
      continue;
    }
    
    // extract variable sequence and merge ("cis" only)
    varSeqForward = sseq1.substr(skipForward + umiLengthForward + constantLengthForward, variableLengthForward);
    varSeqReverse = sseq2.substr(skipReverse + umiLengthReverse + constantLengthReverse, variableLengthReverse);
    varQualForward = squal1.substr(skipForward + umiLengthForward + constantLengthForward, variableLengthForward);
    varQualReverse = squal2.substr(skipReverse + umiLengthReverse + constantLengthReverse, variableLengthReverse);
    
    // convert qualities to int
    for (size_t i = 0; i < variableLengthForward; i++) {
      varIntQualForward[i] = int(varQualForward[i]) - 33;
      varIntQualReverse[i] = int(varQualReverse[i]) - 33;
    }
    
    // for "cis" experiments, fuse forward and reverse reads
    if (experimentType.compare("cis") == 0) {
      // reverse-complement reverse read
      transform(
        begin(varSeqReverse),
        end(varSeqReverse),
        begin(varSeqReverse),
        complement);
      reverse(begin(varSeqReverse),
              end(varSeqReverse));
      reverse(begin(varIntQualReverse),
              end(varIntQualReverse));
      // merge reads (keep base with higher quality) and store in forward read
      for (size_t i = 0; i < varSeqForward.length(); i++) {
        if (varIntQualReverse[i] > varIntQualForward[i]) {
          varSeqForward[i] = varSeqReverse[i];
          varIntQualForward[i] = varIntQualReverse[i];
        }
      }
      // empty the reverse sequence and quality
      varSeqReverse.clear();
      varQualReverse.clear();
      varIntQualReverse.clear();
    }
    
    // filter if the average quality in variable region is too low
    if (std::accumulate(varIntQualForward.begin(), varIntQualForward.end(), 0.0) <
        avePhredMin * variableLengthForward ||
        (!varIntQualReverse.empty() && std::accumulate(varIntQualForward.begin(), varIntQualForward.end(), 0.0) <
          avePhredMin * variableLengthForward)) {
      nAvgVarQualTooLow++;
      continue;
    }
    
    // filter if there are too many N's in variable regions
    if (std::count(varSeqForward.begin(), varSeqForward.end(), 'N') > variableNMax ||
        std::count(varSeqReverse.begin(), varSeqReverse.end(), 'N') > variableNMax) {
      nTooManyNinVar++;
      continue;
    }
    
    // extract UMIs and filter if there are too many N's
    umiSeq = sseq1.substr(skipForward, umiLengthForward) + sseq2.substr(skipReverse, umiLengthReverse);
    if (std::count(umiSeq.begin(), umiSeq.end(), 'N') > umiNMax) {
      nTooManyNinUMI++;
      continue;
    }
    
    // if wildTypeForward is available...
    if (wildTypeForward.compare("") != 0) {
      if (compareToWildtype(varSeqForward, wildTypeForward, varIntQualForward,
                            mutatedPhredMin, nbrMutatedCodonsMax, forbiddenCodons,
                            std::string("f"), nMutQualTooLow, 
                            nTooManyMutCodons, nForbiddenCodons, mutantName)) {
        // read is to be filtered out
        continue;
      }
    }
    
    // if wildTypeReverse is available...
    if (!varSeqReverse.empty() && wildTypeReverse.compare("") != 0) {
      if (compareToWildtype(varSeqReverse, wildTypeReverse, varIntQualReverse,
                            mutatedPhredMin, nbrMutatedCodonsMax, forbiddenCodons,
                            std::string("r"), nMutQualTooLow, 
                            nTooManyMutCodons, nForbiddenCodons, mutantName)) {
        // read is to be filtered out
        continue;
      }
    }
    
    // store the read pair
    nRetain++;
    // ... create final mutant name
    if (mutantName.length() > 0) { // we have a least one mutation
      mutantName.pop_back(); // remove '_' at the end
    } else {
      if (wildTypeForward.compare("") != 0 || (!varSeqReverse.empty() && wildTypeReverse.compare("") != 0)) {
        mutantName = "WT";
      }
    }
    // ... check if mutant already exists in mutantSummary
    if (!varSeqReverse.empty()) { // "trans" experiment
      varSeqForward += (std::string("_") + varSeqReverse);
    }
    if ((mutantSummaryIt = mutantSummary.find(varSeqForward)) != mutantSummary.end()) {
      // ... ... update existing mutantInfo
      (*mutantSummaryIt).second.nReads++;
      (*mutantSummaryIt).second.umi.insert(umiSeq);
    } else {
      // ... ... create mutantInfo instance for this mutant and add it to mutantSummary
      mutantInfo newMutant;
      newMutant.nReads = 1;
      newMutant.umi.insert(umiSeq);
      newMutant.mutantName = mutantName;
      mutantSummary.insert(std::pair<std::string,mutantInfo>(varSeqForward, newMutant));
    }
    
    // for retained reads, count numbers of (mis-)matching bases by Phred quality
    if (constantForward.compare("") != 0) {
      constSeqForward = sseq1.substr(skipForward + umiLengthForward, constantLengthForward);
      constQualForward = squal1.substr(skipForward + umiLengthForward, constantLengthForward);
      for (size_t i = 0; i < constantLengthForward; i++) {
        constIntQualForward[i] = int(constQualForward[i]) - 33;
      }
      tabulateBasesByQual(constSeqForward, constantForward, constIntQualForward,
                          nPhredCorrectForward, nPhredMismatchForward);
    }
    if (constantReverse.compare("") != 0) {
      constSeqReverse = sseq2.substr(skipReverse + umiLengthReverse, constantLengthReverse);
      constQualReverse = squal2.substr(skipReverse + umiLengthReverse, constantLengthReverse);
      for (size_t i = 0; i < constantLengthReverse; i++) {
        constIntQualReverse[i] = int(constQualReverse[i]) - 33;
      }
      tabulateBasesByQual(constSeqReverse, constantReverse, constIntQualReverse,
                          nPhredCorrectReverse, nPhredMismatchReverse);
    }
    
  } // iterate over individual sequence pairs
  if (verbose) {
    Rcout << "collapsed to " << mutantSummary.size() << " elements" << std::endl;
  }

  // clean up
  gzclose(file1);
  gzclose(file2);

  // return results
  size_t dfLen = mutantSummary.size();
  std::vector<std::string> dfSeq(dfLen, ""), dfName(dfLen, "");
  std::vector<int> dfReads(dfLen, 0), dfUmis(dfLen, 0);
  int i = 0;
  for (mutantSummaryIt = mutantSummary.begin(); mutantSummaryIt != mutantSummary.end(); mutantSummaryIt++) {
    dfSeq[i] = (*mutantSummaryIt).first;
    dfName[i] = (*mutantSummaryIt).second.mutantName;
    dfReads[i] = (*mutantSummaryIt).second.nReads;
    dfUmis[i] = (*mutantSummaryIt).second.umi.size();
    i++;
  }
  DataFrame filt = DataFrame::create(Named("nbrTotal") = nTot,
                                     Named("1_nbrAdapter") = nAdapter,
                                     Named("2_nAvgVarQualTooLow") = nAvgVarQualTooLow,
                                     Named("3_nTooManyNinVar") = nTooManyNinVar,
                                     Named("4_nTooManyNinUMI") = nTooManyNinUMI,
                                     Named("5_nMutQualTooLow") = nMutQualTooLow,
                                     Named("6_nTooManyMutCodons") = nTooManyMutCodons,
                                     Named("7_nForbiddenCodons") = nForbiddenCodons,
                                     Named("nbrRetained") = nRetain,
                                     Named("check.names") = false);
  DataFrame df = DataFrame::create(Named("sequence") = dfSeq,
                                   Named("mutantName") = dfName,
                                   Named("nbrReads") = dfReads,
                                   Named("nbrUmis") = dfUmis,
                                   Named("stringsAsFactors") = false);
  DataFrame err = DataFrame::create(Named("PhredQuality") = seq_len(100) - 1,
                                    Named("nbrMatchForward") = nPhredCorrectForward,
                                    Named("nbrMismatchForward") = nPhredMismatchForward,
                                    Named("nbrMatchReverse") = nPhredCorrectReverse,
                                    Named("nbrMismatchReverse") = nPhredMismatchReverse);
  std::vector<std::string> forbiddenCodonsUsed(forbiddenCodons.begin(), forbiddenCodons.end());
  List param;
  param.push_back(experimentType, "experimentType");
  param.push_back(fastqForward, "fastqForward");
  param.push_back(fastqReverse, "fastqReverse");
  param.push_back(skipForward, "skipForward");
  param.push_back(skipReverse, "skipReverse");
  param.push_back(umiLengthForward, "umiLengthForward");
  param.push_back(umiLengthReverse, "umiLengthReverse");
  param.push_back(constantLengthForward, "constantLengthForward");
  param.push_back(constantLengthReverse, "constantLengthReverse");
  param.push_back(variableLengthForward, "variableLengthForward");
  param.push_back(variableLengthReverse, "variableLengthReverse");
  param.push_back(adapterForward, "adapterForward");
  param.push_back(adapterReverse, "adapterReverse");
  param.push_back(wildTypeForward, "wildTypeForward");
  param.push_back(wildTypeReverse, "wildTypeReverse");
  param.push_back(constantForward, "constantForward");
  param.push_back(constantReverse, "constantReverse");
  param.push_back(avePhredMin, "avePhredMin");
  param.push_back(variableNMax, "variableMax");
  param.push_back(umiNMax, "umiNMax");
  param.push_back(nbrMutatedCodonsMax, "nbrMutatedCodonsMax");
  param.push_back(forbiddenCodonsUsed, "forbiddenMutatedCodons");
  param.push_back(mutatedPhredMin, "mutatedPhredMin");
  List L = List::create(Named("parameters") = param,
                        Named("filterSummary") = filt,
                        Named("summaryTable") = df,
                        Named("errorStatistics") = err,
                        Named("experimentType") = experimentType);
  return L;
}
