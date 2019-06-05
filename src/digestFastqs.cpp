#include <Rcpp.h>
#include <cerrno>
#include <zlib.h>
#include <string>
#include <vector>
#include <numeric>
#include <map>
#include <algorithm>    // std::sort

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

// compare position of codons of the form (std::string)"x123NNN_"
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
// enumerate all non-ambiguious codons that match to them
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

// merge forward and reverse sequences
// (keeping the base with maximal Phred quality)
bool mergeReadPair(std::string &varSeqForward, std::vector<int> &varIntQualForward,
                   std::string &varSeqReverse, std::vector<int> &varIntQualReverse) {
  // merge reads (keep base with higher quality) and store in forward read
  for (size_t i = 0; i < varSeqForward.length(); i++) {
    if (varIntQualReverse[i] > varIntQualForward[i]) {
      varSeqForward[i] = varSeqReverse[i];
      varIntQualForward[i] = varIntQualReverse[i];
    }
  }

  return true;
}

// Find closest wild type sequence to a variable sequence
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
List digestFastqsCpp(std::string fastqForward, std::string fastqReverse,
                     bool mergeForwardReverse, bool revComplForward, bool revComplReverse,
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
                     bool verbose = false) {

  // Biostrings::IUPAC_CODE_MAP
  std::map<char,std::vector<char>> IUPAC = initializeIUPAC();

  
  // --------------------------------------------------------------------------
  // open fastq files
  // --------------------------------------------------------------------------
  gzFile file1 = openFastq(fastqForward);
  gzFile file2 = openFastq(fastqReverse);

  
  // --------------------------------------------------------------------------
  // digest reads one by one
  // --------------------------------------------------------------------------
  char seq1[BUFFER_SIZE];
  char qual1[BUFFER_SIZE];
  char seq2[BUFFER_SIZE];
  char qual2[BUFFER_SIZE];
  bool done = false;
  int nTot = 0, nAdapter = 0, nNoPrimer = 0, nAvgVarQualTooLow = 0, nTooManyNinVar = 0, nTooManyNinUMI = 0;
  int nTooManyMutCodons = 0, nForbiddenCodons = 0, nMutQualTooLow = 0, nRetain = 0;
  unsigned int primerPosForward, primerPosReverse;
  std::string varSeqForward, varSeqReverse, varQualForward, varQualReverse, umiSeq;
  std::string constSeqForward, constSeqReverse, constQualForward, constQualReverse;
  std::vector<int> constIntQualForward(constantLengthForward,0);
  std::vector<int> constIntQualReverse(constantLengthReverse,0);
  std::string mutantName;
  std::map<std::string, mutantInfo> mutantSummary;
  std::map<std::string, mutantInfo>::iterator mutantSummaryIt;
  std::vector<int> nPhredCorrectForward(100, 0), nPhredMismatchForward(100, 0);
  std::vector<int> nPhredCorrectReverse(100, 0), nPhredMismatchReverse(100, 0);
  
  // enumerate forbidden codons based on forbiddenMutatedCodons
  std::set<std::string> forbiddenCodonsForward = enumerateCodonsFromIUPAC(forbiddenMutatedCodonsForward, IUPAC, verbose);
  std::set<std::string> forbiddenCodonsReverse = enumerateCodonsFromIUPAC(forbiddenMutatedCodonsReverse, IUPAC, verbose);
  
  // iterate over sequences
  if (verbose) {
    Rcout << "start reading sequences..." << std::endl;
  }
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
    if (nTot % 200000 == 0) { // every 200,000 reads (every ~1.6 seconds)
      Rcpp::checkUserInterrupt(); // ... check for user interrupt
      // ... and give an update
      if (verbose && nTot % 1000000 == 0) {
        Rcout << "    " << nTot << " read pairs read (" 
              << std::setprecision(1) << (100*nRetain/nTot) << " retained)" << std::endl;
      }
    }

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
    
    // search for primer sequences, exclude read pairs if not both primers are found
    if ((primerForward.compare("") != 0 && sseq1.find(primerForward) == std::string::npos) ||
        (primerReverse.compare("") != 0 && sseq2.find(primerReverse) == std::string::npos)) {
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
        varSeqForward.pop_back();
        varQualForward.pop_back();
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
        varSeqForward.pop_back();
        varQualForward.pop_back();
      }
    } else {
      stop("Either expected sequence lengths or a primer sequence must be provided (forward)");
    }
    //reverse
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
        varSeqReverse.pop_back();
        varQualReverse.pop_back();
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
        varSeqReverse.pop_back();
        varQualReverse.pop_back();
      }
    } else {
      stop("Either expected sequence lengths or a primer sequence must be provided (reverse)");
    }

    // check that extracted sequences and qualities are of the right length 
    // (if the read is too short, substr() will just read until the end of it)
    if ((variableLengthForward != (-1) && varSeqForward.length() != variableLengthForward) || 
        (variableLengthReverse != (-1) && varSeqReverse.length() != variableLengthReverse) || 
        (variableLengthForward != (-1) && varQualForward.length() != variableLengthForward) || 
        (variableLengthReverse != (-1) && varQualReverse.length() != variableLengthReverse)) {
      stop("The read is not long enough to extract a variable sequence of the indicated length");
    }

    // reverse complement if requested
    if (revComplForward) {
      transform(
        begin(varSeqForward), end(varSeqForward),
        begin(varSeqForward), complement);
      reverse(begin(varSeqForward), end(varSeqForward));
      reverse(begin(varQualForward), end(varQualForward));
    }
    if (revComplReverse) {
      transform(
        begin(varSeqReverse), end(varSeqReverse),
        begin(varSeqReverse), complement);
      reverse(begin(varSeqReverse), end(varSeqReverse));
      reverse(begin(varQualReverse), end(varQualReverse));
    }
    
    // convert qualities to int
    std::vector<int> varIntQualForward(varSeqForward.length(), 0);
    std::vector<int> varIntQualReverse(varSeqReverse.length(), 0);
    for (size_t i = 0; i < varSeqForward.length(); i++) {
      varIntQualForward[i] = int(varQualForward[i]) - 33;
    }
    for (size_t i = 0; i < varSeqReverse.length(); i++) {
      varIntQualReverse[i] = int(varQualReverse[i]) - 33;
    }

    // if requested, fuse forward and reverse reads
    if (mergeForwardReverse) {
      mergeReadPair(varSeqForward, varIntQualForward,
                    varSeqReverse, varIntQualReverse);
    }
    
    // filter if the average quality in variable region is too low
    if (std::accumulate(varIntQualForward.begin(), varIntQualForward.end(), 0.0) <
        avePhredMinForward * varSeqForward.length() ||
        (!mergeForwardReverse && std::accumulate(varIntQualReverse.begin(), varIntQualReverse.end(), 0.0) <
          avePhredMinReverse * varSeqReverse.length())) {
      nAvgVarQualTooLow++;
      continue;
    }
    
    // filter if there are too many N's in variable regions
    if (std::count(varSeqForward.begin(), varSeqForward.end(), 'N') > variableNMaxForward ||
        (!mergeForwardReverse && std::count(varSeqReverse.begin(), varSeqReverse.end(), 'N') > variableNMaxReverse)) {
      nTooManyNinVar++;
      continue;
    }
    
    // TODO: Allow UMI in only one of the reads?
    // extract UMIs and filter if there are too many N's
    if (umiLengthForward != (-1) && umiLengthReverse != (-1)) {
      umiSeq = sseq1.substr(skipForward, umiLengthForward) + sseq2.substr(skipReverse, umiLengthReverse);
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
    }
    
    // if wildTypeReverse is available...
    if (!mergeForwardReverse && std::string(wildTypeReverse[0]).compare("") != 0) {
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
    }
    
    // store the read pair
    nRetain++;
    // ... create final mutant name
    if (mutantName.length() > 0) { // we have a least one mutation
      mutantName.pop_back(); // remove '_' at the end
    } else {
      if (std::string(wildTypeForward[0]).compare("") != 0 || 
          (!mergeForwardReverse && std::string(wildTypeReverse[0]).compare("") != 0)) {
        mutantName = "WT";
      }
    }
    // ... check if mutant already exists in mutantSummary
    if (!mergeForwardReverse) { // "trans" experiment
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
      
      // check that extracted sequences have the right length
      if (constSeqForward.length() != constantLengthForward || constQualForward.length() != constantLengthForward) {
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
      for (size_t i = 0; i < constantLengthForward; i++) {
        constIntQualForward[i] = int(constQualForward[i]) - 33;
      }
      tabulateBasesByQual(constSeqForward, constantForward, constIntQualForward,
                          nPhredCorrectForward, nPhredMismatchForward);
    }
    if (constantReverse.compare("") != 0) {
      constSeqReverse = sseq2.substr(skipReverse + umiLengthReverse, constantLengthReverse);
      constQualReverse = squal2.substr(skipReverse + umiLengthReverse, constantLengthReverse);
      
      // check that extracted sequences have the right length
      if (constSeqReverse.length() != constantLengthReverse || constQualReverse.length() != constantLengthReverse) {
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
      for (size_t i = 0; i < constantLengthReverse; i++) {
        constIntQualReverse[i] = int(constQualReverse[i]) - 33;
      }
      tabulateBasesByQual(constSeqReverse, constantReverse, constIntQualReverse,
                          nPhredCorrectReverse, nPhredMismatchReverse);
    }
    
  } // iterate over individual sequence pairs
  if (verbose) {
    Rcout << "done reading sequences" << std::endl << "    collapsed to " <<
             mutantSummary.size() << " elements" << std::endl;
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
                                     Named("f1_nbrAdapter") = nAdapter,
                                     Named("f1.5_nNoPrimer") = nNoPrimer,
                                     Named("f2_nAvgVarQualTooLow") = nAvgVarQualTooLow,
                                     Named("f3_nTooManyNinVar") = nTooManyNinVar,
                                     Named("f4_nTooManyNinUMI") = nTooManyNinUMI,
                                     Named("f5_nMutQualTooLow") = nMutQualTooLow,
                                     Named("f6_nTooManyMutCodons") = nTooManyMutCodons,
                                     Named("f7_nForbiddenCodons") = nForbiddenCodons,
                                     Named("nbrRetained") = nRetain);
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
  std::vector<std::string> forbiddenCodonsUsedForward(forbiddenCodonsForward.begin(), forbiddenCodonsForward.end());
  std::vector<std::string> forbiddenCodonsUsedReverse(forbiddenCodonsReverse.begin(), forbiddenCodonsReverse.end());
  List param;
  param.push_back(fastqForward, "fastqForward");
  param.push_back(fastqReverse, "fastqReverse");
  param.push_back(mergeForwardReverse, "mergeForwardReverse");
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
  List L = List::create(Named("parameters") = param,
                        Named("filterSummary") = filt,
                        Named("summaryTable") = df,
                        Named("errorStatistics") = err);
  return L;
}
