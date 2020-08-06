#include <Rcpp.h>
#include <cerrno>
#include <zlib.h>
#include <string>
#include <vector>
#include <numeric>
#include <map>
#include <algorithm>    // std::sort
#include "BKtree_utils.hpp"

#define BUFFER_SIZE 4096

using namespace std::placeholders;
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
      const char *error_string;
      error_string = gzerror(file, &err);
      if (err) {
        stop(error_string);
      }
    }
  }
  return false;
}

// read next four lines from gzipped file and store
// second and fourth in *seq and *qual
// return either:
// - true (reached end of file, I am done)
// - false (not yet reached end of file, not done yet)
// - nothing (encountered an error, fall back to R from reached_end_of_file())
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
  stop("Invalid DNA base character in sequence - aborting");
}   

// initialize IUPAC code table
std::map<char,std::vector<char>> initializeIUPAC() {
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
  return IUPAC;
}

// split string by delimiter (code from https://www.fluentcpp.com/2017/04/21/how-to-split-a-string-in-c/)
std::vector<std::string> split(const std::string& s, char delimiter) {
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (std::getline(tokenStream, token, delimiter)) {
    tokens.push_back(token);
  }
  return tokens;
}

// compare position of codons of the form (std::string)"x.123.NNN_"
// [[Rcpp::export]]
bool compareCodonPositions(std::string a, std::string b, const char mutNameDelimiter) {
  int posa = std::stoi(split(a, mutNameDelimiter)[1]);
  int posb = std::stoi(split(b, mutNameDelimiter)[1]);
  return (posa < posb);
}

// compare read to wildtype sequence,
// identify mutated bases/codons, filter, update counters
// and add to the name
// returns true if the read needs to be filtered out 
// (and a counter has been incremented)
bool compareToWildtype(const std::string varSeq, const std::string wtSeq,
                       const std::vector<int> varIntQual, const double mutatedPhredMin,
                       const unsigned int nbrMutatedCodonsMax, const std::set<std::string> &forbiddenCodons,
                       const std::string codonPrefix, int &nMutQualTooLow, int &nTooManyMutCodons,
                       int &nForbiddenCodons, std::string &mutantName, const std::string mutNameDelimiter) {
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
      mutatedCodons.insert(codonPrefix + mutNameDelimiter + 
        std::to_string((int)(i / 3) + 1) + mutNameDelimiter + 
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
  for (mutatedCodonIt = mutatedCodons.begin(); mutatedCodonIt != mutatedCodons.end(); mutatedCodonIt++) {
    if (forbiddenCodons.find((*mutatedCodonIt).substr((*mutatedCodonIt).length() - 4, 3)) != forbiddenCodons.end()) { // found forbidden codon
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
  std::sort(mutatedCodonsSorted.begin(), mutatedCodonsSorted.end(), std::bind(compareCodonPositions, _1, _2, *(mutNameDelimiter.c_str())));
  for (size_t i = 0; i < mutatedCodonsSorted.size(); i++) {
    mutantName += mutatedCodonsSorted[i];
  }
  
  // if no mutant codons, name as <codonPrefix>.0.WT
  if (mutatedCodonsSorted.size() == 0) {
    mutantName += codonPrefix + mutNameDelimiter + "0" + mutNameDelimiter + "WT_";
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
  std::set<std::string> umi;      // set of unique UMIs observed for read pairs of that mutant
  int nReads;                     // number of reads
  std::set<std::string> sequence; // set of sequences for that mutant
};

// open fastq file and check if it worked
gzFile  openFastq(std::string filename) {
  gzFile file = gzopen(filename.c_str(), "rb");   
  if (!file) {
    if (errno) {
      stop("Failed to open file '", filename,  "': ",
           strerror(errno), " (errno=", errno, ")");
    } else {
      stop("Failed to open file '", filename, "': zlib out of memory");
    }
  }
  return file;
}

// given a vector of codons (potentially containing IUPAC ambiguous bases),
// enumerate all non-ambiguous codons that match to them
std::set<std::string> enumerateCodonsFromIUPAC(CharacterVector forbiddenMutatedCodons,
                                               std::map<char,std::vector<char>> IUPAC, bool verbose) {
  std::set<std::string> forbiddenCodons;
  std::string codon("NNN");
  for (int i = 0; i < forbiddenMutatedCodons.length(); i++) {
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
  return forbiddenCodons;
}

// merge forward and reverse sequences (partially overlapping)
// - don't count N's as a mismatch
// - default values for minOverlap,maxOverlap are zero, which will:
//      set minOverlap,maxOverlap to the length of the shorter read
// - first find optimal overlap without indels
//      for greedy = true, pick the longest overlap with less than maxMismatch errors
//      for greedy = false, pick the one with the highest score := number of overlap bases - number of mismatches in overlap
// - merge the reads, at each overlap position keeping the base with maximal Phred quality
// - store the results into varSeqForward and varIntQualForward
// - if no valid overlap is found within the scope of minOverlap/maxOverlap/maxMismatch, do nothing
//   (varSeqForward and varIntQualForward still correspond to the input values)
// - return false if a valid merge was found and performed, and return true otherwise
bool mergeReadPairPartial(std::string &varSeqForward, std::vector<int> &varIntQualForward,
                          std::string &varSeqReverse, std::vector<int> &varIntQualReverse,
                          size_t minOverlap = 0, size_t maxOverlap = 0, 
                          double maxFracMismatchOverlap = 0,
                          bool greedy = true) {
  // initialize overlap parameters
  size_t lenF = varSeqForward.length(), lenR = varSeqReverse.length();
  if (minOverlap == 0) {
    minOverlap = lenF;
    if (minOverlap > lenR) {
      minOverlap = lenR;
    }
  } else if (minOverlap > lenF || minOverlap > lenR) {
    return true; //  no valid overlap possible
  }
  if (maxOverlap == 0) {
    maxOverlap = lenF;
    if (maxOverlap > lenR) {
      maxOverlap = lenR;
    }
  } else if (maxOverlap < minOverlap) {
    return true; //  no valid overlap possible
  }

  // find overlap (score := number of overlap bases - number of mismatches in overlap)
  size_t o, i, j;
  int bestScore = 0, bestO = -1, score;
  double fracmm;

  for (o = maxOverlap; o >= minOverlap; o--) {
    // calculate score for overlap of o bases
    score = 0;
    for (i = lenF - o, j = 0; i < lenF; i++, j++) {
      if (varSeqForward[i] == varSeqReverse[j] || varSeqForward[i] == 'N' || varSeqReverse[j] == 'N') {
        score++;
      }
    }
    // if valid, store overlap
    fracmm = ((double)o - (double)score)/(double)o;
    if (fracmm <= maxFracMismatchOverlap && score > bestScore) {
      bestO = o;
      bestScore = score;
      // if greedy, stop checking for overlaps
      if (greedy) {
        break;
      }
    }
  }

  // if a valid overlap has been found, merge reads into varSeqForward, varIntQualForward
  if (bestO > 0) {
    // ... grow varSeqForward, varIntQualForward
    varSeqForward.resize(lenF + lenR - bestO);
    varIntQualForward.resize(lenF + lenR - bestO);
    // ... at each overlap position, keep base with higher quality) and store in forward read
    for (i = lenF - bestO, j = 0; i < lenF; i++, j++) {
      if (varIntQualReverse[j] > varIntQualForward[i]) {
        varSeqForward[i] = varSeqReverse[j];
        varIntQualForward[i] = varIntQualReverse[j];
      }
    }
    // ... add non-overlapping part of varSeqReverse, varIntQualReverse
    for (i = lenF, j = bestO; i < lenF + lenR - bestO; i++, j++) {
      varSeqForward[i] = varSeqReverse[j];
      varIntQualForward[i] = varIntQualReverse[j];
    }
    
    return false;
  
  } else {
    return true;
  }
}

// wrapper around mergeReadPairPartial used in unit testing
// (needed because mregeReadPairPartial merges in-place and returns only 'true',
//  and R passes copies of the arguments, so that the results of the merging cannot be "seen")
// [[Rcpp::export]]
List test_mergeReadPairPartial(std::string seqF, std::vector<int> qualF,
                               std::string seqR, std::vector<int> qualR,
                               size_t minOverlap = 0, size_t maxOverlap = 0, 
                               double maxFracMismatchOverlap = 0,
                               bool greedy = true) {
  mergeReadPairPartial(seqF, qualF, seqR, qualR,
                       minOverlap, maxOverlap, maxFracMismatchOverlap, greedy);
  List L = List::create(Named("mergedSeq") = seqF,
                        Named("mergedQual") = qualF);
  return L;
}


void removeEOL(std::string &seq) {
  if (seq.back() == '\n') {
    seq.pop_back();
  }
  if (seq.back() == '\r') {
    seq.pop_back();
  }
}

// Find closest wild type sequence to a variable sequence
// Here, 'closest' is defined as the sequence with the largest number of matching bases
// Assumes that the start of varSeq coincides with the start of each wtSeq
// [[Rcpp::export]]
int findClosestRefSeq(std::string varSeq, Rcpp::StringVector wtSeq) {
  // return index of most similar sequence
  int idx = 0;
  int maxsim = 0;
  int currsim;
  for (int i = 0; i < wtSeq.size(); i++) {
    currsim = 0;
    std::string currSeq = std::string(wtSeq[i]);
    for (size_t j = 0; j < currSeq.length(); j++) {
      if (currSeq[j] == varSeq[j]) {
        currsim++;
      }
    }
    if (currsim > maxsim) {
      idx = i;
      maxsim = currsim;
    }
  }
  return idx;
}

// [[Rcpp::export]]
List digestFastqsCpp(std::vector<std::string> fastqForwardVect,
                     std::vector<std::string> fastqReverseVect,
                     bool mergeForwardReverse, 
                     size_t minOverlap, size_t maxOverlap, 
                     double maxFracMismatchOverlap, bool greedyOverlap,
                     bool revComplForward, bool revComplReverse,
                     int skipForward, int skipReverse,
                     int umiLengthForward, int umiLengthReverse,
                     int constantLengthForward,
                     int constantLengthReverse,
                     int variableLengthForward,
                     int variableLengthReverse,
                     std::string adapterForward, std::string adapterReverse,
                     std::string primerForward, std::string primerReverse,
                     Rcpp::StringVector wildTypeForward, 
                     Rcpp::StringVector wildTypeReverse, 
                     std::string constantForward, std::string constantReverse, 
                     double avePhredMinForward = 20.0, double avePhredMinReverse = 20.0,
                     int variableNMaxForward = 0, int variableNMaxReverse = 0, 
                     int umiNMax = 0,
                     unsigned int nbrMutatedCodonsMaxForward = 1,
                     unsigned int nbrMutatedCodonsMaxReverse = 1,
                     CharacterVector forbiddenMutatedCodonsForward = "NNW",
                     CharacterVector forbiddenMutatedCodonsReverse = "NNW",
                     double mutatedPhredMinForward = 0.0,
                     double mutatedPhredMinReverse = 0.0,
                     std::string mutNameDelimiter = ".",
                     double variableCollapseMaxDist = 0.0,
                     double umiCollapseMaxDist = 0.0,
                     int maxNReads = -1, bool verbose = false) {

  
  // Biostrings::IUPAC_CODE_MAP
  std::map<char,std::vector<char>> IUPAC = initializeIUPAC();

  // --------------------------------------------------------------------------
  // digest reads one by one
  // --------------------------------------------------------------------------
  char seq1[BUFFER_SIZE];
  char qual1[BUFFER_SIZE];
  char seq2[BUFFER_SIZE];
  char qual2[BUFFER_SIZE];
  bool noReverse;
  int nTot = 0, nAdapter = 0, nNoPrimer = 0, nReadTooShort = 0, nNoValidOverlap = 0, nAvgVarQualTooLow = 0, nTooManyNinVar = 0, nTooManyNinUMI = 0;
  int nTooManyMutCodons = 0, nForbiddenCodons = 0, nMutQualTooLow = 0, nRetain = 0;
  unsigned int primerPosForward, primerPosReverse;
  std::string varSeqForward, varSeqReverse, varQualForward, varQualReverse, umiSeq;
  std::string constSeqForward, constSeqReverse, constQualForward, constQualReverse;
  std::string mutantName;
  std::map<std::string, mutantInfo> mutantSummary;
  std::map<std::string, mutantInfo>::iterator mutantSummaryIt;
  std::vector<int> nPhredCorrectForward(100, 0), nPhredMismatchForward(100, 0);
  std::vector<int> nPhredCorrectReverse(100, 0), nPhredMismatchReverse(100, 0);
  
  // enumerate forbidden codons based on forbiddenMutatedCodons
  std::set<std::string> forbiddenCodonsForward = enumerateCodonsFromIUPAC(forbiddenMutatedCodonsForward, IUPAC, verbose);
  std::set<std::string> forbiddenCodonsReverse = enumerateCodonsFromIUPAC(forbiddenMutatedCodonsReverse, IUPAC, verbose);

  // --------------------------------------------------------------------------
  // iterate over fastq files
  // --------------------------------------------------------------------------
  for (size_t f = 0; f < fastqForwardVect.size(); f++) {
    bool done = false;
    std::string fastqForward = fastqForwardVect[f];
    std::string fastqReverse = fastqReverseVect[f];
    
    // --------------------------------------------------------------------------
    // open fastq files
    // --------------------------------------------------------------------------
    gzFile file1 = openFastq(fastqForward);
    gzFile file2 = NULL;
    if (fastqReverse.compare("") != 0) {
      file2 = openFastq(fastqReverse);
    }
    
    // iterate over sequences
    if (verbose) {
      Rcout << "start reading sequences for file " <<
               (fastqReverse.compare("") != 0 ? "pair " : "") << (f + 1) <<
               " of " << fastqForwardVect.size() << "..." << std::endl;
    }
    while (done == false) {
      mutantName = ""; // start with empty name
      
      // read sequence pair
      done = get_next_seq(file1, seq1, qual1);
      if (fastqReverse.compare("") != 0) {
        done = (done || get_next_seq(file2, seq2, qual2));
      }
      
      // if maxNReads has been reached, break
      if (maxNReads != (-1) && nTot >= maxNReads) {
        done = true;
      }
      
      // check if end-of-file was reached
      if (done) {
        break;
      }
      
      // update counters
      nTot++;
      if (nTot % 200000 == 0) { // every 200,000 reads (every ~1.6 seconds)
        Rcpp::checkUserInterrupt(); // ... check for user interrupt
        // ... and give an update
        if (verbose && nTot % 1000000 == 0) {
          Rcout << "    " << nTot << " read pairs read ("
                << std::setprecision(3) << (100*((double)nRetain/nTot)) << "% retained)" << std::endl;
        }
      }

      // convert C char* to C++ string
      std::string sseq1(seq1);
      std::string squal1(qual1);
      std::string sseq2, squal2;
      if (fastqReverse.compare("") != 0) {
        sseq2 = seq2;
        squal2 = qual2;
      }
      
      // search for adapter sequences and filter read pairs
      if ((adapterForward.compare("") != 0 && sseq1.find(adapterForward) != std::string::npos) ||
          (fastqReverse.compare("") != 0 && adapterReverse.compare("") != 0 && sseq2.find(adapterReverse) != std::string::npos)) {
        nAdapter++;
        continue;
      }
      
      // search for primer sequences, exclude read pairs if not both primers are found
      if ((primerForward.compare("") != 0 && sseq1.find(primerForward) == std::string::npos) ||
          (fastqReverse.compare("") != 0 && primerReverse.compare("") != 0 && sseq2.find(primerReverse) == std::string::npos)) {
        nNoPrimer++;
        continue;
      }
      
      // extract variable sequence
      // forward
      if (skipForward != (-1) && umiLengthForward != (-1) && constantLengthForward != (-1)) {
        // if skipForward, umiLengthForward, constantLengthForward are all not -1, extract by position
        if (variableLengthForward != (-1)) {
          // variable length given
          varSeqForward = sseq1.substr(skipForward + umiLengthForward + constantLengthForward, variableLengthForward);
          varQualForward = squal1.substr(skipForward + umiLengthForward + constantLengthForward, variableLengthForward);
        } else {
          // variable length not given - take the rest of the read
          varSeqForward = sseq1.substr(skipForward + umiLengthForward + constantLengthForward, sseq1.length());
          varQualForward = squal1.substr(skipForward + umiLengthForward + constantLengthForward, squal1.length());
        }
      } else if (primerForward.compare("") != 0) {
        // otherwise, primer sequence should be given, and will be used as a trigger
        primerPosForward = sseq1.find(primerForward);
        if (variableLengthForward != (-1)) {
          // variable length given
          varSeqForward = sseq1.substr(primerPosForward + primerForward.length(), variableLengthForward);
          varQualForward = squal1.substr(primerPosForward + primerForward.length(), variableLengthForward);
        } else {
          // variable length not given - take the rest of the read
          varSeqForward = sseq1.substr(primerPosForward + primerForward.length(), sseq1.length());
          varQualForward = squal1.substr(primerPosForward + primerForward.length(), squal1.length());
        }
      } else {
        stop("Either expected sequence lengths or a primer sequence must be provided (forward)");
      }
      // reverse
      if (fastqReverse.compare("") != 0) {
        if (skipReverse != (-1) && umiLengthReverse != (-1) && constantLengthReverse != (-1)) {
          // if skipReverse, umiLengthReverse, constantLengthReverse are all not -1, extract by position
          if (variableLengthReverse != (-1)) {
            // variable length given
            varSeqReverse = sseq2.substr(skipReverse + umiLengthReverse + constantLengthReverse, variableLengthReverse);
            varQualReverse = squal2.substr(skipReverse + umiLengthReverse + constantLengthReverse, variableLengthReverse);
          } else {
            // variable length not given - take the rest of the read
            varSeqReverse = sseq2.substr(skipReverse + umiLengthReverse + constantLengthReverse, sseq2.length());
            varQualReverse = squal2.substr(skipReverse + umiLengthReverse + constantLengthReverse, squal2.length());
          }
        } else if (primerReverse.compare("") != 0) {
          // otherwise, primer sequence should be given, and will be used as a trigger
          primerPosReverse = sseq2.find(primerReverse);
          if (variableLengthReverse != (-1)) {
            // variable length given
            varSeqReverse = sseq2.substr(primerPosReverse + primerReverse.length(), variableLengthReverse);
            varQualReverse = squal2.substr(primerPosReverse + primerReverse.length(), variableLengthReverse);
          } else {
            // variable length not given - take the rest of the read
            varSeqReverse = sseq2.substr(primerPosReverse + primerReverse.length(), sseq2.length());
            varQualReverse = squal2.substr(primerPosReverse + primerReverse.length(), squal2.length());
          }
        } else {
          stop("Either expected sequence lengths or a primer sequence must be provided (reverse)");
        }
      }
      
      // check if the last character(s) are new line, if so remove them
      removeEOL(varSeqForward);
      removeEOL(varQualForward);
      if (fastqReverse.compare("") != 0) {
        removeEOL(varSeqReverse);
        removeEOL(varQualReverse);
      }
      
      // check that extracted sequences and qualities are of the right length 
      // (if the read is too short, substr() will just read until the end of it)
      // if the read sequence is too short, discard the read pair
      // don't raise an error, since this could cause unintended problems when extracting sequence parts based on primers
      // (there may be a 'random' primer match in the middle of the read)
      if ((variableLengthForward != (-1) && varSeqForward.length() != (size_t) variableLengthForward) || 
          (fastqReverse.compare("") != 0 && variableLengthReverse != (-1) && varSeqReverse.length() != (size_t) variableLengthReverse) || 
          (variableLengthForward != (-1) && varQualForward.length() != (size_t) variableLengthForward) || 
          (fastqReverse.compare("") != 0 && variableLengthReverse != (-1) && varQualReverse.length() != (size_t) variableLengthReverse)) {
        nReadTooShort++;
        continue;
      }
  
      // reverse complement if requested
      if (revComplForward) {
        transform(
          begin(varSeqForward), end(varSeqForward),
          begin(varSeqForward), complement);
        reverse(begin(varSeqForward), end(varSeqForward));
        reverse(begin(varQualForward), end(varQualForward));
      }
      if (fastqReverse.compare("") != 0 && revComplReverse) {
        transform(
          begin(varSeqReverse), end(varSeqReverse),
          begin(varSeqReverse), complement);
        reverse(begin(varSeqReverse), end(varSeqReverse));
        reverse(begin(varQualReverse), end(varQualReverse));
      }
      
      // convert qualities to int
      std::vector<int> varIntQualForward(varSeqForward.length(), 0);
      for (size_t i = 0; i < varSeqForward.length(); i++) {
        varIntQualForward[i] = int(varQualForward[i]) - 33;
      }
      
      std::vector<int> varIntQualReverse;
      if (fastqReverse.compare("") != 0) {
        varIntQualReverse.resize(varSeqReverse.length()); // set number of elements
        std::fill(varIntQualReverse.begin(), varIntQualReverse.end(), 0); // fill with zero
        for (size_t i = 0; i < varSeqReverse.length(); i++) {
          varIntQualReverse[i] = int(varQualReverse[i]) - 33;
        }
      }
      
      // if requested, fuse forward and reverse reads
      if (mergeForwardReverse) {
        if (mergeReadPairPartial(varSeqForward, varIntQualForward,
                                 varSeqReverse, varIntQualReverse,
                                 minOverlap, maxOverlap, maxFracMismatchOverlap,
                                 greedyOverlap)) {
          // read should be filtered out - no valid overlap found
          nNoValidOverlap++;
          continue;
        }
      }
      // if no reverse sequence was provided, hereafter it is identical to the situation 
      // where forward and reverse reads were merged
      noReverse = mergeForwardReverse || fastqReverse.compare("") == 0;
      
      // filter if the average quality in variable region is too low
      if (std::accumulate(varIntQualForward.begin(), varIntQualForward.end(), 0.0) <
          avePhredMinForward * varSeqForward.length() ||
          (!noReverse && std::accumulate(varIntQualReverse.begin(), varIntQualReverse.end(), 0.0) <
            avePhredMinReverse * varSeqReverse.length())) {
        nAvgVarQualTooLow++;
        continue;
      }
      
      // filter if there are too many N's in variable regions
      if (std::count(varSeqForward.begin(), varSeqForward.end(), 'N') > variableNMaxForward ||
          (!noReverse && std::count(varSeqReverse.begin(), varSeqReverse.end(), 'N') > variableNMaxReverse)) {
        nTooManyNinVar++;
        continue;
      }
      
      // TODO: Allow UMI in only one of the reads?
      // extract UMIs and filter if there are too many N's
      if (umiLengthForward != (-1)) {
        if (fastqReverse.compare("") != 0 && umiLengthReverse != (-1)) {
          umiSeq = sseq1.substr(skipForward, umiLengthForward) + sseq2.substr(skipReverse, umiLengthReverse);
        } else {
          umiSeq = sseq1.substr(skipForward, umiLengthForward);
        }
        if (std::count(umiSeq.begin(), umiSeq.end(), 'N') > umiNMax) {
          nTooManyNinUMI++;
          continue;
        }
      } else {
        umiSeq = "";
      }
      
      // if wildTypeForward is available...
      if (std::string(wildTypeForward[0]).compare("") != 0) {
        std::vector<std::string> refNamesForward = wildTypeForward.attr("names");
        int idxForward = findClosestRefSeq(varSeqForward, wildTypeForward);
        std::string wtForward = std::string(wildTypeForward[idxForward]);
        std::string wtNameForward = std::string(refNamesForward[idxForward]);
        if (compareToWildtype(varSeqForward, wtForward, varIntQualForward,
                              mutatedPhredMinForward, nbrMutatedCodonsMaxForward, forbiddenCodonsForward,
                              wtNameForward, nMutQualTooLow, 
                              nTooManyMutCodons, nForbiddenCodons, mutantName, mutNameDelimiter)) {
          // read is to be filtered out
          continue;
        }
      } else if (varSeqForward.length() > 0) { // variable seq, but no reference -> add variable seq to mutantName
        mutantName += (varSeqForward + std::string("_"));
      }
      
      // if wildTypeReverse is available...
      if (!noReverse && std::string(wildTypeReverse[0]).compare("") != 0) {
        std::vector<std::string> refNamesReverse = wildTypeReverse.attr("names");
        int idxReverse = findClosestRefSeq(varSeqReverse, wildTypeReverse);
        std::string wtReverse = std::string(wildTypeReverse[idxReverse]);
        std::string wtNameReverse = std::string(refNamesReverse[idxReverse]);
        if (compareToWildtype(varSeqReverse, wtReverse, varIntQualReverse,
                              mutatedPhredMinReverse, nbrMutatedCodonsMaxReverse, forbiddenCodonsReverse,
                              wtNameReverse, nMutQualTooLow, 
                              nTooManyMutCodons, nForbiddenCodons, mutantName, mutNameDelimiter)) {
          // read is to be filtered out
          continue;
        }
      } else if (!noReverse && varSeqReverse.length() > 0) { // variable seq, but no reference -> add variable seq to mutantName
        mutantName += (varSeqReverse + std::string("_"));
      }
      
      // store the read pair
      nRetain++;
      // ... create final mutant name
      if (mutantName.length() > 0) { // we have a least one mutation, or sequence-based name
        mutantName.pop_back(); // remove '_' at the end
      } else {
        if (std::string(wildTypeForward[0]).compare("") != 0 || 
            (!noReverse && std::string(wildTypeReverse[0]).compare("") != 0)) {
          mutantName = "WT";
        }
      }
      
      if (!noReverse) { // "trans" experiment
        varSeqForward += (std::string("_") + varSeqReverse);
      }
      // ... check if mutant already exists in mutantSummary
      if ((mutantSummaryIt = mutantSummary.find(mutantName)) != mutantSummary.end()) {
        // ... ... update existing mutantInfo
        (*mutantSummaryIt).second.nReads++;
        (*mutantSummaryIt).second.umi.insert(umiSeq);
        (*mutantSummaryIt).second.sequence.insert(varSeqForward);
      } else {
        // ... ... create mutantInfo instance for this mutant and add it to mutantSummary
        mutantInfo newMutant;
        newMutant.nReads = 1;
        newMutant.umi.insert(umiSeq);
        newMutant.sequence.insert(varSeqForward);
        mutantSummary.insert(std::pair<std::string,mutantInfo>(mutantName, newMutant));
      }
      
      // for retained reads, count numbers of (mis-)matching bases by Phred quality
      if (constantForward.compare("") != 0) {
        constSeqForward = sseq1.substr(skipForward + umiLengthForward, constantLengthForward);
        constQualForward = squal1.substr(skipForward + umiLengthForward, constantLengthForward);
        
        // check that extracted sequences have the right length
        if (constSeqForward.length() != (size_t) constantLengthForward || 
            constQualForward.length() != (size_t) constantLengthForward) {
          stop("The read is not long enough to extract a forward constant sequence of the indicated length");
        }
        
        // reverse complement if requested
        if (revComplForward) {
          transform(begin(constSeqForward), end(constSeqForward),
                    begin(constSeqForward), complement);
          reverse(constSeqForward.begin(), constSeqForward.end());
          reverse(constQualForward.begin(), constQualForward.end());
        }
        
        // populate an integer vector of base qualities
        std::vector<int> constIntQualForward(constantLengthForward, 0);
        for (size_t i = 0; i < (size_t) constantLengthForward; i++) {
          constIntQualForward[i] = int(constQualForward[i]) - 33;
        }
        tabulateBasesByQual(constSeqForward, constantForward, constIntQualForward,
                            nPhredCorrectForward, nPhredMismatchForward);
      }
      
      if (fastqReverse.compare("") != 0 && constantReverse.compare("") != 0) {
        constSeqReverse = sseq2.substr(skipReverse + umiLengthReverse, constantLengthReverse);
        constQualReverse = squal2.substr(skipReverse + umiLengthReverse, constantLengthReverse);
        
        // check that extracted sequences have the right length
        if (constSeqReverse.length() != (size_t) constantLengthReverse || 
            constQualReverse.length() != (size_t) constantLengthReverse) {
          stop("The read is not long enough to extract a reverse constant sequence of the indicated length");
        }
        
        // reverse (complement) sequence and quality string
        if (revComplReverse) {
          transform(begin(constSeqReverse), end(constSeqReverse),
                    begin(constSeqReverse), complement);
          reverse(constSeqReverse.begin(), constSeqReverse.end());
          reverse(constQualReverse.begin(), constQualReverse.end());
        }
        
        // populate an integer vector of base qualities
        std::vector<int> constIntQualReverse(constantLengthReverse,0);
        for (size_t i = 0; i < (size_t) constantLengthReverse; i++) {
          constIntQualReverse[i] = int(constQualReverse[i]) - 33;
        }
        tabulateBasesByQual(constSeqReverse, constantReverse, constIntQualReverse,
                            nPhredCorrectReverse, nPhredMismatchReverse);
      }
    } // iterate over individual sequence pairs

    // clean up
    gzclose(file1);
    if (fastqReverse.compare("") != 0) {
      gzclose(file2);
    }
    
    if (verbose) {
      Rcout << "done reading sequences" << std::endl;
    }
  } // iterate over fastq files
  if (verbose) {
    Rcout << "retained " << mutantSummary.size() << " unique features" << std::endl;
  }

  // collapse similar variable sequences in mutantSummary
  if (variableCollapseMaxDist > 0.0) {
    if (std::string(wildTypeForward[0]).compare("") != 0 ||
        std::string(wildTypeReverse[0]).compare("") != 0) {
      warning("Skipping variable sequence collapsing because wildtype reference sequence(s) are given");
    } else {
      // get sequence length
      mutantSummaryIt = mutantSummary.begin();
      size_t seqlen = (*mutantSummaryIt).first.length();

      // calculate Hamming distance tolerance
      int tol;
      if (variableCollapseMaxDist >= 1.0) {
        tol = (int)variableCollapseMaxDist;
      } else {
        tol = (int)(variableCollapseMaxDist *
          ((*mutantSummaryIt).first.find("_") != std::string::npos ? seqlen-1 : seqlen));
      }

      if (verbose) {
        Rcout << "start collapsing variable sequences (tolerance: " << tol << ")...";
      }

      // sort mutantSummary decreasingly by read count
      // ... create an empty intermediate vector
      std::vector<std::pair<std::string,mutantInfo>> vec;
      std::vector<std::pair<std::string,mutantInfo>>::iterator vecIt;
      // copy key-value pairs from mutantSummary to vec
      std::copy(mutantSummary.begin(), mutantSummary.end(),
                std::back_inserter<std::vector<std::pair<std::string,mutantInfo>>>(vec));
      // ... sort vec by decreasing order of pair.second.nReads
      //     (if second values are equal, order by the pair's first value)
      std::sort(vec.begin(), vec.end(),
                [](const std::pair<std::string,mutantInfo>& l,
                   const std::pair<std::string,mutantInfo>& r) {
                  if (l.second.nReads != r.second.nReads)
                    return l.second.nReads > r.second.nReads;
                  return l.first < r.first;
                  });

      // store sequences (from names) in BK tree
      BKtree tree;
      for (vecIt = vec.begin(); vecIt != vec.end(); vecIt++) {
        if ((*vecIt).first.length() != seqlen) {
          warning("Skipping variable sequence collapsing because reads are not all of the same length");
          tree.remove_all();
          break;
        } else {
          tree.insert((*vecIt).first);
        }
      }
      vec.clear(); // remove temporary vector

      if (tree.size > 0) {
        std::string querySeq, collapsedName;
        std::vector<std::string> simSeqs;
        std::map<std::string, std::string> single2collapsed;
        
        // start querying in the order of tree.items (ordered decreasingly by nReads)
        while (tree.size > 0) {
          querySeq = tree.first();
          simSeqs = tree.search(querySeq, tol);
          for (size_t i = 0; i < simSeqs.size(); i++) {
            single2collapsed[simSeqs[i]] = querySeq;
            tree.remove(simSeqs[i]);
          }
        }

        // group into sets of similar sequences
        std::map<std::string, mutantInfo> collapsedMutantSummary;
        std::map<std::string, mutantInfo>::iterator collapsedMutantSummaryIt;
        for (mutantSummaryIt = mutantSummary.begin(); mutantSummaryIt != mutantSummary.end(); mutantSummaryIt++) {
          collapsedName = single2collapsed[(*mutantSummaryIt).first];
          if ((collapsedMutantSummaryIt = collapsedMutantSummary.find(collapsedName)) != collapsedMutantSummary.end()) {
            // ... fuse with existing mutantInfo
            (*collapsedMutantSummaryIt).second.nReads += (*mutantSummaryIt).second.nReads;
            (*collapsedMutantSummaryIt).second.umi.insert((*mutantSummaryIt).second.umi.begin(),
                                                          (*mutantSummaryIt).second.umi.end());
            (*collapsedMutantSummaryIt).second.sequence.insert((*mutantSummaryIt).second.sequence.begin(),
                                                               (*mutantSummaryIt).second.sequence.end());
          } else {
            // ... insert first mutantInfo
            collapsedMutantSummary.insert(std::pair<std::string,mutantInfo>(collapsedName, (*mutantSummaryIt).second));
          }
        }
        if (verbose) {
          Rcout << "done (reduced from " << mutantSummary.size() << " to " << collapsedMutantSummary.size() << ")" << std::endl;
        }
        mutantSummary = collapsedMutantSummary;
      }
    }
  }
  
  // collapse similar UMI sequences in each variable sequence mutantSummary
  if (umiCollapseMaxDist > 0.0) {
    // calculate Hamming distance tolerance
    int tol;
    if (umiCollapseMaxDist >= 1.0) {
      tol = (int)umiCollapseMaxDist;
    } else {
      tol = (int)(umiCollapseMaxDist * (*(*mutantSummary.begin()).second.umi.begin()).length());
    }
    
    if (verbose) {
      Rcout << "start collapsing UMIs (tolerance: " << tol << ")...";
    }

    // store sequences (from names) in BK tree
    BKtree tree;
    std::set<std::string>::iterator umiIt;
    for (mutantSummaryIt = mutantSummary.begin(); mutantSummaryIt != mutantSummary.end(); mutantSummaryIt++) {
      if ((*mutantSummaryIt).second.umi.size() == 1) {
        continue;
      } else {
        tree.remove_all();
        for (umiIt = (*mutantSummaryIt).second.umi.begin(); umiIt != (*mutantSummaryIt).second.umi.end(); umiIt++) {
          tree.insert((*umiIt));
        }

        std::vector<std::string> simSeqs;
        std::set<std::string> collapsedUmis;
        while (tree.size > 0) {
          simSeqs = tree.search(tree.first(), tol);
          collapsedUmis.insert(tree.first());
          for (size_t i = 0; i < simSeqs.size(); i++) {
            tree.remove(simSeqs[i]);
          }
        }
        
        (*mutantSummaryIt).second.umi = collapsedUmis;
      }
    }
    
    if (verbose) {
      Rcout << "done" << std::endl;
    }
  }

  // return results
  size_t dfLen = mutantSummary.size();
  std::vector<std::string> dfSeq(dfLen, ""), dfName(dfLen, "");
  std::vector<int> dfReads(dfLen, 0), dfUmis(dfLen, 0);
  int i = 0;
  for (mutantSummaryIt = mutantSummary.begin(); mutantSummaryIt != mutantSummary.end(); mutantSummaryIt++) {
    // collapse all sequences associated with the mutant
    std::vector<std::string> sequenceVector((*mutantSummaryIt).second.sequence.begin(), 
                                            (*mutantSummaryIt).second.sequence.end());
    std::string collapsedSequence = "";
    for (size_t i = 0; i < sequenceVector.size(); i++) {
      collapsedSequence += sequenceVector[i] + ",";
    }
    collapsedSequence.pop_back(); // remove final ","
    dfName[i] = (*mutantSummaryIt).first;
    dfSeq[i] = collapsedSequence;
    dfReads[i] = (*mutantSummaryIt).second.nReads;
    dfUmis[i] = (*mutantSummaryIt).second.umi.size();
    i++;
  }
  DataFrame filt = DataFrame::create(Named("nbrTotal") = nTot,
                                     Named("f1_nbrAdapter") = nAdapter,
                                     Named("f2_nbrNoPrimer") = nNoPrimer,
                                     Named("f3_nbrReadTooShort") = nReadTooShort,
                                     Named("f4_nbrNoValidOverlap") = nNoValidOverlap,
                                     Named("f5_nbrAvgVarQualTooLow") = nAvgVarQualTooLow,
                                     Named("f6_nbrTooManyNinVar") = nTooManyNinVar,
                                     Named("f7_nbrTooManyNinUMI") = nTooManyNinUMI,
                                     Named("f8_nbrMutQualTooLow") = nMutQualTooLow,
                                     Named("f9_nbrTooManyMutCodons") = nTooManyMutCodons,
                                     Named("f10_nbrForbiddenCodons") = nForbiddenCodons,
                                     Named("nbrRetained") = nRetain);
  DataFrame df = DataFrame::create(Named("mutantName") = dfName,
                                   Named("sequence") = dfSeq,
                                   Named("nbrReads") = dfReads,
                                   Named("nbrUmis") = dfUmis,
                                   Named("stringsAsFactors") = false);
  DataFrame err = DataFrame::create(Named("PhredQuality") = seq_len(100) - 1,
                                    Named("nbrMatchForward") = nPhredCorrectForward,
                                    Named("nbrMismatchForward") = nPhredMismatchForward,
                                    Named("nbrMatchReverse") = nPhredCorrectReverse,
                                    Named("nbrMismatchReverse") = nPhredMismatchReverse);
  std::vector<std::string> forbiddenCodonsUsedForward(forbiddenCodonsForward.begin(), forbiddenCodonsForward.end());
  std::vector<std::string> forbiddenCodonsUsedReverse(forbiddenCodonsReverse.begin(), forbiddenCodonsReverse.end());
  List param;
  param.push_back(fastqForwardVect, "fastqForward");
  param.push_back(fastqReverseVect, "fastqReverse");
  param.push_back(mergeForwardReverse, "mergeForwardReverse");
  param.push_back(minOverlap, "minOverlap");
  param.push_back(maxOverlap, "maxOverlap");
  param.push_back(maxFracMismatchOverlap, "maxFracMismatchOverlap");
  param.push_back(greedyOverlap, "greedyOverlap");
  param.push_back(revComplForward, "revComplForward");
  param.push_back(revComplReverse, "revComplReverse");
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
  param.push_back(primerForward, "primerForward");
  param.push_back(primerReverse, "primerReverse");
  param.push_back(wildTypeForward, "wildTypeForward");
  param.push_back(wildTypeReverse, "wildTypeReverse");
  param.push_back(constantForward, "constantForward");
  param.push_back(constantReverse, "constantReverse");
  param.push_back(avePhredMinForward, "avePhredMinForward");
  param.push_back(avePhredMinReverse, "avePhredMinReverse");
  param.push_back(variableNMaxForward, "variableNMaxForward");
  param.push_back(variableNMaxReverse, "variableNMaxReverse");
  param.push_back(umiNMax, "umiNMax");
  param.push_back(nbrMutatedCodonsMaxForward, "nbrMutatedCodonsMaxForward");
  param.push_back(nbrMutatedCodonsMaxReverse, "nbrMutatedCodonsMaxReverse");
  param.push_back(forbiddenCodonsUsedForward, "forbiddenMutatedCodonsForward");
  param.push_back(forbiddenCodonsUsedReverse, "forbiddenMutatedCodonsReverse");
  param.push_back(mutatedPhredMinForward, "mutatedPhredMinForward");
  param.push_back(mutatedPhredMinReverse, "mutatedPhredMinReverse");
  param.push_back(mutNameDelimiter, "mutNameDelimiter");
  param.push_back(variableCollapseMaxDist, "variableCollapseMaxDist");
  param.push_back(umiCollapseMaxDist, "umiCollapseMaxDist");
  param.push_back(maxNReads, "maxNReads");
  List L = List::create(Named("parameters") = param,
                        Named("filterSummary") = filt,
                        Named("summaryTable") = df,
                        Named("errorStatistics") = err);
  return L;
}
