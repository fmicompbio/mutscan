#include <Rcpp.h>
#include <cerrno>
#include <zlib.h>
#include <string>
#include <vector>
#include <numeric>
#include <map>
#include <algorithm>    // std::sort
#include "BKtree_utils.h"

#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

// #define BUFFER_SIZE 65536 // maximum length of a single read + 2

// define constants that are used below
#define NO_SIMILAR_REF    -1 // no similar enough wildtype sequence was found
#define TOO_MANY_BEST_REF -2 // too many equally good hits among the WT sequences
#define QUALITY_OFFSET    33

// FastqBuffer_utils needs the BUFFER_SIZE to be defined
#include "FastqBuffer_utils.h"

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
  // Check if we have read until a newline character. Otherwise, the read is 
  // too long -> break
  if (std::string(ret).back() != '\n') {
    stop("Encountered a read exceeding the maximal allowed length");
  }
  return false;
}

// read next four lines from gzipped file and store
// second and fourth in *seq and *qual
// return either:
// - true (reached end of file, I am done)
// - false (not yet reached end of file, not done yet)
// - nothing (encountered an error, fall back to R from reached_end_of_file())
bool get_next_seq(gzFile file, char *seq, char *qual, size_t BUFFER_SIZE) {
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

// initialize one- to three-letter amino acid conversion table
std::map<char,std::string> initializeThreeAA() {
    std::map<char,std::string> threeAA;
    threeAA['A'] = std::string("Ala");
    threeAA['R'] = std::string("Arg");
    threeAA['N'] = std::string("Asn");
    threeAA['D'] = std::string("Asp");
    threeAA['C'] = std::string("Cys");
    threeAA['Q'] = std::string("Gln");
    threeAA['E'] = std::string("Glu");
    threeAA['G'] = std::string("Gly");
    threeAA['H'] = std::string("His");
    threeAA['I'] = std::string("Ile");
    threeAA['L'] = std::string("Leu");
    threeAA['K'] = std::string("Lys");
    threeAA['M'] = std::string("Met");
    threeAA['F'] = std::string("Phe");
    threeAA['P'] = std::string("Pro");
    threeAA['S'] = std::string("Ser");
    threeAA['T'] = std::string("Thr");
    threeAA['W'] = std::string("Trp");
    threeAA['Y'] = std::string("Tyr");
    threeAA['V'] = std::string("Val");
    threeAA['*'] = std::string("*");
    return threeAA;
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

// translate sequence
// [[Rcpp::export]]
std::string translateString(const std::string& s) {
  const char* tr = "KQE*TPASRRG*ILVLNHDYTPASSRGCILVFKQE*TPASRRGWMLVLNHDYTPASSRGCILVF";
  std::string aa = "";
  size_t i = 0, j = 0;
  const int pow4[] = {1, 4, 16};
  int codonVal = 0;
  
  for (i = 0; i < s.size(); i++) {
    if (s[i] == '_') {
      j = 0;
      aa.push_back('_');
    } else {
      switch(s[i]) {
        case 'A':
        case 'a':
          // codonVal += 0 * pow4[j];
          break;
        
        case 'C':
        case 'c':
          codonVal += 1 * pow4[j];
          break;
        
        case 'G':
        case 'g':
          codonVal += 2 * pow4[j];
          break;
        
        case 'T':
        case 't':
        case 'u':
        case 'U':
          codonVal += 3 * pow4[j];
          break;
          
        default:
          codonVal += 64;
          break;
      }
      
      if (j == 2) {
        if (codonVal > 63) {
          aa.push_back('X');
        } else {
          aa.push_back(tr[codonVal]);
        }
        j = 0;
        codonVal = 0;
      } else {
        j++;
      }
    
    }
    
  }
  return aa;
}

// check if string ends with another string
// see https://stackoverflow.com/questions/874134/find-out-if-string-ends-with-another-string-in-c
bool hasEnding (std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

// Helper function to create the ID for a single mutation group (delins or 
// single mutation)
std::string makeSingleBaseHGVSid(int posMin, int posMax, const std::string wtSeq, 
                                 const std::string varSeq) {
    std::string id = "";
    if (posMax > posMin) {
        // delins
        id += std::to_string(posMin) + "_" + std::to_string(posMax) + 
            "delins" + 
            varSeq.substr((size_t)(posMin - 1), (size_t)(posMax - posMin + 1)) + ";";
    } else {
        // single base mutation
        id += std::to_string(posMin) + wtSeq[posMin - 1] + 
            ">" + varSeq[posMin - 1] + ";";
    }
    return id;
}

// [[Rcpp::export]]
std::string makeBaseHGVS(const std::vector<std::string> mutationsSorted, 
                         const std::string mutNameDelimiter, 
                         const std::string wtSeq, const std::string varSeq) {
  int prevPosMin = -5, prevPosMax = -5; // something far enough from the first mutation position to not cause a delins
  bool moreThanOne = false;
  int pos;
  std::string mutId = "";
    
  if (mutationsSorted.size() == 0) {
    // no mutations - return only _ (all variants should end with this character,
    // as it will be stripped off later)
    return std::string("_");
  } else {
    for (size_t i = 0; i < mutationsSorted.size(); i++) {
      // get the mutated position
      pos = std::stoi(split(mutationsSorted[i], *(mutNameDelimiter.c_str()))[1]);
      if (pos - prevPosMax > 2) {
        // need to start a new mutation group
        // first finish the previous one and add to mutId
        if (prevPosMin > 0) {
          // an actual mutation (the initialization value was negative)
          if (mutId != "") {
            // not the first mutation we add
            moreThanOne = true;
          }
          mutId += makeSingleBaseHGVSid(prevPosMin, prevPosMax, wtSeq, varSeq);
        }
        prevPosMin = pos;
        prevPosMax = pos;
      } else {
        prevPosMax = pos;
      }
    }
    // process last mutation
    if (mutId != "") {
      // not the first mutation we add
      moreThanOne = true;
    }
    mutId += makeSingleBaseHGVSid(prevPosMin, prevPosMax, wtSeq, varSeq);
    mutId.pop_back(); // remove final ;
    if (moreThanOne) {
      mutId = "[" + mutId + "]_";
    } else {
      mutId += "_";
    }
  }
  return mutId;
}

std::string makeAAHGVS(const std::vector<std::string> mutationsSorted, 
                       const std::string mutNameDelimiter, 
                       const std::string wtSeq, 
                       std::map<char, std::string> &threeAA) {
  int pos;
  std::string mutAA;
  char wtAA;
  std::string mutId = "";
    
  if (mutationsSorted.size() == 0) {
    // no mutations - return only _ (all variants should end with this character,
    // as it will be stripped off later)
    return std::string("_");
  } else {
    if (mutationsSorted.size() > 1) {
      mutId += "[";
    }
    for (size_t i = 0; i < mutationsSorted.size(); i++) {
      // get the mutated position and the mutated amino acid
      pos = std::stoi(split(mutationsSorted[i], *(mutNameDelimiter.c_str()))[1]);
      mutAA = split(mutationsSorted[i], *(mutNameDelimiter.c_str()))[2];
      // get the corresponding wt amino acid
      wtAA = wtSeq[pos - 1];
      // assemble mutation identifier
      mutId += "(" + threeAA[wtAA] + std::to_string(pos) + threeAA[*(mutAA.c_str())] + ");";
    }
    mutId.pop_back(); // remove final ;
    if (mutationsSorted.size() > 1) {
      mutId += "]";
    }
    mutId += "_";
    return mutId;
  }
}

// [[Rcpp::export]]
std::string test_makeAAHGVS(const std::vector<std::string> mutationsSorted, 
                            const std::string mutNameDelimiter, 
                            const std::string wtSeq) {
    std::map<char,std::string> threeAA = initializeThreeAA();
    return makeAAHGVS(mutationsSorted, mutNameDelimiter, wtSeq, threeAA);
} 
    
// compare read to wildtype sequence,
// identify mutated bases/codons, filter, update counters
// and add to the name (mutantName* may be empty if no mutations are found,
// except mutantName*HGVS will be set to the nearest wildtype name := `codonPrefix`)
// returns true if the read needs to be filtered out
// (and a counter has been incremented)
// TODO:
// - instead of conditionally creating mutantName,
//   build mutantNameBase and mutantNameCodon in parallel and
//   conditionally set to mutantName base at the end
// - at the end, take mutantNameBase and mutantNameAA and
//   create mutantNameBaseHGVS and mutantNameAAHGVS from it (helper function
//   that takes name and wt/mutant sequences as input):
//   * always using `codonPrefix` as the reference sequence name
//   * if it contains no ':c', add it for mutantNameBaseHGVS
//   * replace ':c' by ':p' for mutantNameAAHGVS (or add if there is no ':c')
//   * use `threeAA` for mapping one- to three-letter aa
//   * follow future standard: encode base substitutions individually, combine
//     all that are spaced less than 2nt appart into a delins
//   * make sure all names (including *HGVS) always end with a final '_'
//     (except non-HGVS names if there are no mutations)
// - fix tests
// - add tests for *HGVS
bool compareToWildtype(const std::string varSeq, const std::string wtSeq,
                       const std::vector<int> varIntQual, const double mutatedPhredMin,
                       const int nbrMutatedCodonsMax, const std::set<std::string> &forbiddenCodons,
                       const std::string codonPrefix, const int nbrMutatedBasesMax,
                       int &nMutQualTooLow, int &nTooManyMutCodons, int &nForbiddenCodons,
                       int &nTooManyMutBases, std::string &mutantName, 
                       std::string &mutantNameBase, std::string &mutantNameCodon,
                       std::string &mutantNameBaseHGVS, std::string &mutantNameAA,
                       std::string &mutantNameAAHGVS, int &nMutBases,
                       int &nMutCodons, int &nMutAAs, std::set<std::string> &mutationTypes, 
                       const std::string mutNameDelimiter,
                       const bool collapseToWT,
                       std::map<char, std::string> &threeAA) {
  // exactly one of nbrMutatedCodonsMax or nbrMutatedBasesMax should be -1 (checked in the R code).
  // the one that is not -1 will be used for filtering and naming the mutant

  std::set<std::string> mutatedBases, mutatedCodons, mutatedAAs;
  std::set<std::string>::iterator mutatedCodonIt;
  bool hasLowQualMutation, hasForbidden;
  std::string varCodon, wtCodon, varAA, wtAA;

  // create prefixes for HGVS identifiers
  // if codonPrefix doesn't end with :c, add it -> basePrefixHGVS
  // replace :c by :p in basePrefixHGVS -> aaPrefixHGVS
  std::string basePrefixHGVS = codonPrefix, aaPrefixHGVS = codonPrefix;
  if (!hasEnding(codonPrefix, ":c")) {
      basePrefixHGVS = basePrefixHGVS + ":c";
      aaPrefixHGVS = aaPrefixHGVS + ":p";
  } else {
      aaPrefixHGVS.pop_back();
      aaPrefixHGVS = aaPrefixHGVS + "p";
  }
  
  // filter if there are too many mutated codons
  // mutatedCodons.clear();
  hasLowQualMutation = false;
  for (size_t i = 0; i < varSeq.length(); i++) {
    if (varSeq[i] != wtSeq[i]) { // found mismatching base
      // record if the mutated base quality is below a threshold
      if (varIntQual[i] < mutatedPhredMin) {
        hasLowQualMutation = true;
        break;
      }
      // add codon or base (depending on which to consider) to mutatedCodonsOrBases
      // if (nbrMutatedCodonsMax != (-1)) {
      // consider codons
      varCodon = varSeq.substr((int)(i / 3) * 3, 3);
      wtCodon = wtSeq.substr((int)(i / 3) * 3, 3);
      mutatedCodons.insert(codonPrefix + mutNameDelimiter +
        std::to_string((int)(i / 3) + 1) + mutNameDelimiter +
        varCodon +
        std::string("_"));
      // } else {
        // nbrMutatedBasesMax != -1
        // consider individual bases
      mutatedBases.insert(codonPrefix + mutNameDelimiter +
        std::to_string((int)(i) + 1) + mutNameDelimiter +
        varSeq.substr((int)(i), 1) +
        std::string("_"));
      // }
      varAA = translateString(varCodon);
      wtAA = translateString(wtCodon);
      if (varAA != wtAA) {
        // add to mutatedAA
        mutatedAAs.insert(codonPrefix + mutNameDelimiter + 
          std::to_string((int)(i / 3) + 1) + mutNameDelimiter +
          varAA + 
          std::string("_"));
        // nonsynonymous or stop
        if (varAA == "*") {
          mutationTypes.insert("stop");
        } else {
          mutationTypes.insert("nonsynonymous");
        }
      } else {
        // synonymous
        mutationTypes.insert("silent");
      }
    }
  }
  if (hasLowQualMutation) {
#ifdef _OPENMP
    #pragma omp atomic
#endif
    nMutQualTooLow++;
    return true;
  }
  // check if there are too many mutated codons/bases
  if (nbrMutatedCodonsMax != (-1)) {
    // consider codons
    if ((int)mutatedCodons.size() > nbrMutatedCodonsMax) {
#ifdef _OPENMP
      #pragma omp atomic
#endif
      nTooManyMutCodons++;
      return true;
    }
  } else {
    //consider bases
    if ((int)mutatedBases.size() > nbrMutatedBasesMax) {
#ifdef _OPENMP
      #pragma omp atomic
#endif
      nTooManyMutBases++;
      return true;
    }
  }

  // check if there are forbidden codons
  // if (nbrMutatedCodonsMax != (-1)) {
  hasForbidden = false;
  for (mutatedCodonIt = mutatedCodons.begin(); mutatedCodonIt != mutatedCodons.end(); mutatedCodonIt++) {
    if (forbiddenCodons.find((*mutatedCodonIt).substr((*mutatedCodonIt).length() - 4, 3)) != forbiddenCodons.end()) { // found forbidden codon
      hasForbidden = true;
      break;
    }
  }
  if (hasForbidden) {
#ifdef _OPENMP
    #pragma omp atomic
#endif
    nForbiddenCodons++;
    return true;
  }
  // }

  nMutBases += (int)mutatedBases.size();
  nMutCodons += (int)mutatedCodons.size();
  nMutAAs += (int)mutatedAAs.size();

  // create name for mutant
  if (collapseToWT) {
    mutantNameBase += codonPrefix + "_";
    mutantNameCodon += codonPrefix + "_";
    mutantName += codonPrefix + "_";
    mutantNameAA += codonPrefix + "_";
    mutantNameBaseHGVS += basePrefixHGVS + "_";
    mutantNameAAHGVS += aaPrefixHGVS + "_";
  } else {
    std::vector<std::string> mutatedBasesSorted, mutatedCodonsSorted, mutatedAAsSorted;
    
    mutatedBasesSorted.assign(mutatedBases.begin(), mutatedBases.end());
    std::sort(mutatedBasesSorted.begin(), mutatedBasesSorted.end(),
              std::bind(compareCodonPositions, _1, _2, *(mutNameDelimiter.c_str())));
    
    mutatedCodonsSorted.assign(mutatedCodons.begin(), mutatedCodons.end());
    std::sort(mutatedCodonsSorted.begin(), mutatedCodonsSorted.end(),
              std::bind(compareCodonPositions, _1, _2, *(mutNameDelimiter.c_str())));
    
    mutatedAAsSorted.assign(mutatedAAs.begin(), mutatedAAs.end());
    std::sort(mutatedAAsSorted.begin(), mutatedAAsSorted.end(),
              std::bind(compareCodonPositions, _1, _2, *(mutNameDelimiter.c_str())));
    
    for (size_t i = 0; i < mutatedBasesSorted.size(); i++) {
      mutantNameBase += mutatedBasesSorted[i];
      if (nbrMutatedCodonsMax == (-1)) {
        mutantName += mutatedBasesSorted[i];
      }
    }
    for (size_t i = 0; i < mutatedCodonsSorted.size(); i++) {
      mutantNameCodon += mutatedCodonsSorted[i];
      if (nbrMutatedCodonsMax != (-1)) {
        mutantName += mutatedCodonsSorted[i];
      }
    }
    for (size_t i = 0; i < mutatedAAsSorted.size(); i++) {
      mutantNameAA += mutatedAAsSorted[i];
    }
    // if no mutant codons, name as <codonPrefix>.0.WT
    if (mutatedBasesSorted.size() == 0) {
      mutantNameBase += codonPrefix + mutNameDelimiter + "0" + mutNameDelimiter + "WT_";
      mutantNameCodon += codonPrefix + mutNameDelimiter + "0" + mutNameDelimiter + "WT_";
      mutantName += codonPrefix + mutNameDelimiter + "0" + mutNameDelimiter + "WT_";
      mutantNameBaseHGVS += basePrefixHGVS + "_";
    } else {
      mutantNameBaseHGVS += basePrefixHGVS + "." + 
          makeBaseHGVS(mutatedBasesSorted, mutNameDelimiter, wtSeq, varSeq);
    }
    if (mutatedAAsSorted.size() == 0) {
      mutantNameAA += codonPrefix + mutNameDelimiter + "0" + mutNameDelimiter + "WT_";
      mutantNameAAHGVS += aaPrefixHGVS + "_";
    } else {
      mutantNameAAHGVS += aaPrefixHGVS + "." + 
          makeAAHGVS(mutatedAAsSorted, mutNameDelimiter, 
                     translateString(wtSeq), threeAA);
    }
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
#ifdef _OPENMP
      #pragma omp atomic
#endif
      nPhredMismatch[constIntQual[i]]++;
    } else {                          // found match
#ifdef _OPENMP
      #pragma omp atomic
#endif
      nPhredCorrect[constIntQual[i]]++;
    }
  }
  return true;
}


// stores information about retained mutants
struct mutantInfo {
  std::set<std::string> umi;      // set of unique UMIs observed for read pairs of that mutant
  int nReads;                     // number of reads
  int maxNReads;                  // maximum number of reads for an individual mutant, in case of collapsing multiple
                                  // mutants into one entry
  std::set<int> nMutBases;
  std::set<int> nMutCodons;
  std::set<int> nMutAAs;
  std::set<std::string> mutationTypes;
  std::set<std::string> sequence; // set of sequences for that mutant
  std::set<std::string> mutantNameBase;
  std::set<std::string> mutantNameCodon;
  std::set<std::string> mutantNameBaseHGVS;
  std::set<std::string> mutantNameAA;
  std::set<std::string> mutantNameAAHGVS;
  std::set<std::string> sequenceAA;
  std::string varLengths;         // lengths of individual variable segments (e.g. "10,20" or "10,20_20,10")
};

// open fastq file and check if it worked
gzFile  openFastq(std::string filename, const char* mode = "rb") {
  gzFile file = gzopen(filename.c_str(), mode);
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

// given a read and information about its 'composition', extract UMI, constant sequence,
// variable sequence
// S - skip
// U - umi
// V - variable
// C - constant
// P - primer
bool decomposeRead(const std::string sseq,
                   const std::string squal,
                   const std::string elements,
                   const std::vector<int> elementLengths,
                   const std::vector<std::string> primerSeqs,
                   std::string &umiSeq,
                   std::string &varSeq, std::string &varQual,
                   std::vector<int> &varLengths,
                   std::string &constSeq, std::string &constQual,
                   int &nNoPrimer, int &nReadWrongLength) {
  // check if elements contains a P (primer)
  size_t pPos = elements.find('P');  // location of P in the elements string
  if (pPos != std::string::npos) {

    // search for the primer
    bool foundprimer = false;
    for (size_t j = 0; j < primerSeqs.size(); j++) {
      size_t primerPos = sseq.find(primerSeqs[j]); // primer position in read
      if (primerPos != std::string::npos) {
        foundprimer = true;
        // primer found - split read in parts before/after primer and check separately
        // none of them can have a P in the elements substring, since only one P is allowed in
        // the original elements string
        // part before the primer
        std::vector<int> elementsPre = std::vector<int>(elementLengths.begin(),
                                                        elementLengths.begin() + pPos);
        bool pre = decomposeRead(sseq.substr(0, primerPos), squal.substr(0, primerPos),
                                 elements.substr(0, pPos), elementsPre, primerSeqs,
                                 umiSeq, varSeq, varQual, varLengths, constSeq,
                                 constQual, nNoPrimer, nReadWrongLength);
        if (!pre) {
          return false;
        }
        // part after the primer
        std::vector<int> elementsPost = std::vector<int>(elementLengths.begin() + pPos + 1,
                                                         elementLengths.end());
        bool post = decomposeRead(sseq.substr(primerPos + primerSeqs[j].size(), sseq.size()),
                                  squal.substr(primerPos + primerSeqs[j].size(), squal.size()),
                                  elements.substr(pPos + 1, elements.size()),
                                  elementsPost, primerSeqs,
                                  umiSeq, varSeq, varQual, varLengths, constSeq,
                                  constQual, nNoPrimer, nReadWrongLength);
        if (!post) {
          return false;
        }
        break;
      }
    }
    // no primer found - read should not be considered further
    if (!foundprimer) {
#ifdef _OPENMP
      #pragma omp atomic
#endif
      nNoPrimer++;
      return false;
    }
  } else {
    // no primer, extract read parts from given lengths
    // iterate through elements

    // sum of the element lengths (including -1), will be
    // used to derive the length if not specified
    int totalElementLength = std::accumulate(elementLengths.begin(),
                                             elementLengths.end(), 0);
    size_t currentPos = 0; // current position in the string
    size_t currentLength; // length of the current element

    // check that the length is ok
    if (elementLengths.size() > 0) { // don't check if there are no elements,
                                     // e.g. before the primer
      if (std::find(elementLengths.begin(), elementLengths.end(), -1) == elementLengths.end()) {
        // no -1's - sum of element lengths must be equal to the sequence length
        if (totalElementLength != (int)sseq.size()) {
#ifdef _OPENMP
          #pragma omp atomic
#endif
          nReadWrongLength++;
          return false;
        }
      } else {
        // -1 present - sum of element lengths must be < sequence length
        if (totalElementLength >= (int)sseq.size()) {
#ifdef _OPENMP
          #pragma omp atomic
#endif
          nReadWrongLength++;
          return false;
        }
      }
    }

    for (size_t i = 0; i < elements.size(); i++) { // go through elements
      if (elementLengths[i] == -1) {
        // get element length as the total read length minus the lengths of the specified parts
        currentLength = sseq.size() - totalElementLength - 1;
      } else {
        currentLength = elementLengths[i];
      }
      if (elements[i] == 'U') { // add to umiSeq
        umiSeq += sseq.substr(currentPos, currentLength);
      } else if (elements[i] == 'C') { // add to constant seq
        constSeq += sseq.substr(currentPos, currentLength);
        constQual += squal.substr(currentPos, currentLength);
      } else if (elements[i] == 'V') { // add to variable seq
        varLengths.push_back((int)currentLength);
        varSeq += sseq.substr(currentPos, currentLength);
        varQual += squal.substr(currentPos, currentLength);
      }
      currentPos += currentLength;
    }
  }

  return true; // keep read
}

// wrapper around decomposeRead used in unit testing
// (needed because decomposeRead merges in-place and returns only 'true' or 'false',
//  and R passes copies of the arguments, so that the results of the extraction cannot be "seen")
// [[Rcpp::export]]
List test_decomposeRead(const std::string sseq,
                        const std::string squal,
                        const std::string elements,
                        const std::vector<int> elementLengths,
                        const std::vector<std::string> primerSeqs,
                        std::string umiSeq,
                        std::string varSeq, std::string varQual,
                        std::vector<int> varLengths,
                        std::string constSeq, std::string constQual,
                        int nNoPrimer, int nReadWrongLength) {
  decomposeRead(sseq, squal, elements, elementLengths, primerSeqs,
                umiSeq, varSeq, varQual, varLengths, constSeq, constQual,
                nNoPrimer, nReadWrongLength);
  List L = List::create(Named("umiSeq") = umiSeq,
                        Named("varSeq") = varSeq,
                        Named("varQual") = varQual,
                        Named("varLengths") = varLengths,
                        Named("constSeq") = constSeq,
                        Named("constQual") = constQual,
                        Named("nNoPrimer") = nNoPrimer,
                        Named("nReadWrongLength") = nReadWrongLength);
  return L;
}


// merge forward and reverse sequences (partially overlapping)
// - don't count N's as a mismatch
// - default values for minOverlap,maxOverlap are zero, which will:
//      set minOverlap,maxOverlap to the length of the shorter read
// - first find optimal overlap without indels
//      for greedy = true, pick the longest overlap with less than maxMismatch errors
//      for greedy = false, pick the one with the highest score := number of overlap bases - number of mismatches in overlap
// - merge the reads, at each overlap position keeping the base with maximal Phred quality
// - check if the merging position (based on comparing sequences) is compatible with
//      the variable segment lengths (break points are aligned) and update varLengthsForward accordingly
// - store the results into varSeqForward and varIntQualForward
// - if no valid overlap is found within the scope of minOverlap/maxOverlap/maxMismatch, do nothing
//   (varSeqForward and varIntQualForward still correspond to the input values)
// - return false if a valid merge was found and performed, and return true otherwise
bool mergeReadPairPartial(std::string &varSeqForward, std::vector<int> &varIntQualForward,
                          std::string &varSeqReverse, std::vector<int> &varIntQualReverse,
                          std::vector<int> &varLengthsForward, std::vector<int> &varLengthsReverse,
                          size_t minOverlap = 0, size_t maxOverlap = 0,
                          size_t minMergedLength = 0, size_t maxMergedLength = 0,
                          double maxFracMismatchOverlap = 0,
                          bool greedy = true) {
  // initialize overlap parameters
  size_t lenF = varSeqForward.length(), lenR = varSeqReverse.length();
  if (minOverlap == 0) {
    minOverlap = lenF;
    if (minOverlap > lenR) {
      minOverlap = lenR;
    }
  }
  // if maxMergedLength is 0, set it to the sum of lenF and lenR
  // (the merged sequence can't be longer than that)
  if (maxMergedLength == 0) {
    maxMergedLength = lenF + lenR;
  }
  // if maxMergedLength is specified, adjust allowed minOverlap
  if (minOverlap < (lenF + lenR - maxMergedLength)) {
    minOverlap = (lenF + lenR - maxMergedLength);
  }
  if (minOverlap > lenF || minOverlap > lenR) {
    return true; //  no valid overlap possible
  }

  if (maxOverlap == 0) {
    maxOverlap = lenF;
    if (maxOverlap > lenR) {
      maxOverlap = lenR;
    }
  }
  // if minMergedLength is 0, set it to the larger of lenF and lenR
  // (the merged sequence can't be shorter than that)
  if (minMergedLength == 0) {
    minMergedLength = lenF;
    if (minMergedLength < lenR) {
      minMergedLength = lenR;
    }
  }
  // if minMergedLength is specified, adjust allowed maxOverlap
  if (maxOverlap > (lenF + lenR - minMergedLength)) {
    maxOverlap = (lenF + lenR - minMergedLength);
  }
  if (maxOverlap < minOverlap) {
    return true; //  no valid overlap possible
  }

  // find overlap (score := number of overlap bases - number of mismatches in overlap)
  size_t o, i, j;
  size_t cumpos = 0, bestO = 0, matchedlen = 0;
  int bestScore = -1, score;
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
    // ... check if variable segments lengths' breaks are aligned
    o = 0;
    for (i = 0; i < varLengthsForward.size(); i++) {
      cumpos += varLengthsForward[i];
      if (cumpos > lenF - bestO) { // break within overlap
        if (cumpos < lenF &&
            (int)(cumpos - (lenF - bestO) - matchedlen) != varLengthsReverse[i - o]) {
          return true; // no valid overlap possible
        } else if (cumpos < lenF) {
          matchedlen += varLengthsReverse[i - o];
        }
      } else { // before overlap -> increase break offset for reverse (forward[i] <=> reverse[i - o])
        o++;
      }
    }
    // ... ... is the next reverse read break still in the overlap region?
    i = varLengthsForward.size() - 1;
    if (varLengthsReverse.size() > i - o && matchedlen + varLengthsReverse[i - o] < bestO) {
      return true;
    }
    // ... ... all variable lengths breaks are aligned -> create segment lengths vector for the merged
    varLengthsForward.resize(o); // keep the first o element from forward (before the overlap)
    varLengthsForward.push_back(varLengthsReverse[0] +
        (lenF - bestO - std::accumulate(varLengthsForward.begin(),varLengthsForward.end(), 0)));
    for (i = 1; i < varLengthsReverse.size(); i++) {
      varLengthsForward.push_back(varLengthsReverse[i]);
    }
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
                               std::vector<int> lenF, std::vector<int> lenR,
                               size_t minOverlap = 0, size_t maxOverlap = 0,
                               size_t minMergedLength = 0, size_t maxMergedLength = 0,
                               double maxFracMismatchOverlap = 0,
                               bool greedy = true) {
  bool result = mergeReadPairPartial(seqF, qualF, seqR, qualR, lenF, lenR,
                                     minOverlap, maxOverlap,
                                     minMergedLength, maxMergedLength,
                                     maxFracMismatchOverlap, greedy);
  List L = List::create(Named("return") = result,
                        Named("mergedSeq") = seqF,
                        Named("mergedQual") = qualF,
                        Named("mergedLengths") = lenF);
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
int findClosestRefSeq(std::string &varSeq, std::vector<std::string> &wtSeq,
                      size_t upperBoundMismatch, int &sim) {
  // return index of most similar sequence
  int idx = NO_SIMILAR_REF;
  int maxsim = 0;
  int nbrbesthits = 0;
  int currsim;
  for (size_t i = 0; i < wtSeq.size(); i++) {
    currsim = 0;
    std::string currSeq = wtSeq[i];
    for (size_t j = 0; j < currSeq.length(); j++) {
      currsim += (currSeq[j] == varSeq[j]);
    }
    if (((int)varSeq.size() - currsim <= (int)upperBoundMismatch)) {
      if (currsim == maxsim) {
        nbrbesthits++;
      } else if (currsim > maxsim) {
        nbrbesthits = 1;
        idx = i;
        maxsim = currsim;
      }
    }
  }
  sim = maxsim;
  if (nbrbesthits > 1) {
    return TOO_MANY_BEST_REF;
  } else {
    return idx;
  }
}

// Find closest wild type sequence to a variable sequence
// Similar to findClosestRefSeq, but implements an early stopping
// criterion if it is clear that the similarity is not going
// be large enough
// [[Rcpp::export]]
int findClosestRefSeqEarlyStop(std::string &varSeq, std::vector<std::string> &wtSeq,
                               size_t upperBoundMismatch, int &sim) {
  // return index of most similar sequence
  int idx = NO_SIMILAR_REF;
  int maxsim = -1;
  int nbrbesthits = 0;
  int currsim;
  size_t minl;
  for (size_t i = 0; i < wtSeq.size(); i++) {
    currsim = 0;
    std::string currSeq = wtSeq[i];
    minl = std::min(varSeq.size(), currSeq.size());
    for (size_t j = 0; j < minl; j++) {
      // if currsim + (minl - j) >= varSeq.size() - upperBoundMismatch, we can still continue looking
      // if varSeq is longer than minl, the additional positions are considered mismatches
      if (currsim < (int)(j - minl + varSeq.size() - upperBoundMismatch)) {
        // no chance to reach the minimal similarity - break
        break;
      }
      currsim += (currSeq[j] == varSeq[j]);
    }
    if (((int)varSeq.size() - currsim <= (int)upperBoundMismatch)) {
      if (currsim == maxsim) {
        nbrbesthits++;
      } else if (currsim > maxsim) {
        nbrbesthits = 1;
        idx = (int)i;
        maxsim = currsim;
      }
    }
  }
  sim = maxsim;
  if (nbrbesthits > 1) {
    return TOO_MANY_BEST_REF;
  } else {
    return idx;
  }
}

// Find closest wild type sequence to a variable sequence
// Here, 'closest' is defined as the sequence with the largest number of matching bases
// Assumes that the start of varSeq coincides with the start of each wtSeq
// Input: a variable sequence, a BKtree and a map (mapping entries in the tree to
// their indexes in the original vector)
int findClosestRefSeqTree(std::string &varSeq, BKtree &wtTree,
                          std::map<std::string, int> &wtTreeIdx,
                          size_t upperBoundMismatch, int &sim) {
  // return index of most similar sequence
  std::string idx;
  std::vector<std::string> searchResults;
  for (size_t i = 0; i < std::min(upperBoundMismatch + 1, varSeq.size()); i++) {
    // find all sequences in the tree within distance i of varSeq
    searchResults = wtTree.search(varSeq, i);
    if (searchResults.size() > 0) {
      sim = varSeq.size() - i;
      // if we find multiple best hits, skip the read (return TOO_MANY_BEST_REF)
      if (searchResults.size() > 1) {
        return TOO_MANY_BEST_REF;
      } else {
        idx = searchResults[0];
        return wtTreeIdx[idx];
      }
    }
  }
  return NO_SIMILAR_REF;
}

// Convert mutant summary to DataFrame
DataFrame mutantSummaryToDataFrame(std::map<std::string, mutantInfo> mutantSummary) {
    std::map<std::string, mutantInfo>::iterator mutantSummaryIt;
    
    size_t dfLen = mutantSummary.size();
    std::vector<std::string> dfSeq(dfLen, ""), dfName(dfLen, "");
    std::vector<int> dfReads(dfLen, 0), dfUmis(dfLen, 0), dfMaxReads(dfLen, 0);
    std::vector<std::string> dfMutBases(dfLen, ""), dfMutCodons(dfLen, ""); 
    std::vector<std::string> dfMutAAs(dfLen, ""), dfVarLengths(dfLen, "");
    std::vector<std::string> dfMutantNameBase(dfLen, ""), dfMutantNameCodon(dfLen, "");
    std::vector<std::string> dfMutantNameBaseHGVS(dfLen, ""), dfMutantNameAA(dfLen, "");
    std::vector<std::string> dfMutantNameAAHGVS(dfLen, ""), dfSeqAA(dfLen, "");
    std::vector<std::string> dfMutationTypes(dfLen, "");
    int j = 0;
    for (mutantSummaryIt = mutantSummary.begin(); mutantSummaryIt != mutantSummary.end(); mutantSummaryIt++) {
        // collapse all sequences associated with the mutant
        std::vector<std::string> sequenceVector((*mutantSummaryIt).second.sequence.begin(),
                                                (*mutantSummaryIt).second.sequence.end());
        std::string collapsedSequence = "";
        for (size_t i = 0; i < sequenceVector.size(); i++) {
            collapsedSequence += sequenceVector[i] + ",";
        }
        if (!collapsedSequence.empty()) {
            collapsedSequence.pop_back(); // remove final ","
        }
        
        // mutantNameBase
        std::vector<std::string> mutantNameBaseVector((*mutantSummaryIt).second.mutantNameBase.begin(),
                                                      (*mutantSummaryIt).second.mutantNameBase.end());
        std::string collapsedMutantNameBase = "";
        for (size_t i = 0; i < mutantNameBaseVector.size(); i++) {
            collapsedMutantNameBase += mutantNameBaseVector[i] + ",";
        }
        if (!collapsedMutantNameBase.empty()) {
            collapsedMutantNameBase.pop_back(); // remove final ","
        }
        
        // mutantNameCodon
        std::vector<std::string> mutantNameCodonVector((*mutantSummaryIt).second.mutantNameCodon.begin(),
                                                       (*mutantSummaryIt).second.mutantNameCodon.end());
        std::string collapsedMutantNameCodon = "";
        for (size_t i = 0; i < mutantNameCodonVector.size(); i++) {
            collapsedMutantNameCodon += mutantNameCodonVector[i] + ",";
        }
        if (!collapsedMutantNameCodon.empty()) {
            collapsedMutantNameCodon.pop_back(); // remove final ","
        }
        
        // mutantNameBaseHGVS
        std::vector<std::string> mutantNameBaseHGVSVector((*mutantSummaryIt).second.mutantNameBaseHGVS.begin(),
                                                          (*mutantSummaryIt).second.mutantNameBaseHGVS.end());
        std::string collapsedMutantNameBaseHGVS = "";
        for (size_t i = 0; i < mutantNameBaseHGVSVector.size(); i++) {
            collapsedMutantNameBaseHGVS += mutantNameBaseHGVSVector[i] + ",";
        }
        if (!collapsedMutantNameBaseHGVS.empty()) {
            collapsedMutantNameBaseHGVS.pop_back(); // remove final ","
        }
        
        // mutantNameAA
        std::vector<std::string> mutantNameAAVector((*mutantSummaryIt).second.mutantNameAA.begin(),
                                                    (*mutantSummaryIt).second.mutantNameAA.end());
        std::string collapsedMutantNameAA = "";
        for (size_t i = 0; i < mutantNameAAVector.size(); i++) {
            collapsedMutantNameAA += mutantNameAAVector[i] + ",";
        }
        if (!collapsedMutantNameAA.empty()) {
            collapsedMutantNameAA.pop_back(); // remove final ","
        }
        
        // mutantNameAAHGVS
        std::vector<std::string> mutantNameAAHGVSVector((*mutantSummaryIt).second.mutantNameAAHGVS.begin(),
                                                        (*mutantSummaryIt).second.mutantNameAAHGVS.end());
        std::string collapsedMutantNameAAHGVS = "";
        for (size_t i = 0; i < mutantNameAAHGVSVector.size(); i++) {
            collapsedMutantNameAAHGVS += mutantNameAAHGVSVector[i] + ",";
        }
        if (!collapsedMutantNameAAHGVS.empty()) {
            collapsedMutantNameAAHGVS.pop_back(); // remove final ","
        }
        
        // mutationTypes
        std::vector<std::string> mutationTypesVector((*mutantSummaryIt).second.mutationTypes.begin(),
                                                     (*mutantSummaryIt).second.mutationTypes.end());
        std::string collapsedMutationTypes = "";
        for (size_t i = 0; i < mutationTypesVector.size(); i++) {
            collapsedMutationTypes += mutationTypesVector[i] + ",";
        }
        if (!collapsedMutationTypes.empty()) {
            collapsedMutationTypes.pop_back(); // remove final ","
        }
        
        // collapse all aa sequences associated with the mutant
        std::vector<std::string> sequenceAAVector((*mutantSummaryIt).second.sequenceAA.begin(),
                                                  (*mutantSummaryIt).second.sequenceAA.end());
        std::string collapsedSequenceAA = "";
        for (size_t i = 0; i < sequenceAAVector.size(); i++) {
            collapsedSequenceAA += sequenceAAVector[i] + ",";
        }
        if (!collapsedSequenceAA.empty()) {
            collapsedSequenceAA.pop_back(); // remove final ","
        }
        
        // collapse all observed nMutBases/nMutCodons
        std::vector<int> nMutBasesVector((*mutantSummaryIt).second.nMutBases.begin(),
                                         (*mutantSummaryIt).second.nMutBases.end());
        std::string collapsedNMutBases = "";
        for (size_t i = 0; i < nMutBasesVector.size(); i++) {
            collapsedNMutBases += (std::to_string(nMutBasesVector[i]) + ",");
        }
        if (!collapsedNMutBases.empty()) {
            collapsedNMutBases.pop_back();
        }
        
        // codons
        std::vector<int> nMutCodonsVector((*mutantSummaryIt).second.nMutCodons.begin(),
                                          (*mutantSummaryIt).second.nMutCodons.end());
        std::string collapsedNMutCodons = "";
        for (size_t i = 0; i < nMutCodonsVector.size(); i++) {
            collapsedNMutCodons += std::to_string(nMutCodonsVector[i]) + ",";
        }
        if (!collapsedNMutCodons.empty()) {
            collapsedNMutCodons.pop_back();
        }
        
        // AAs
        std::vector<int> nMutAAsVector((*mutantSummaryIt).second.nMutAAs.begin(),
                                       (*mutantSummaryIt).second.nMutAAs.end());
        std::string collapsedNMutAAs = "";
        for (size_t i = 0; i < nMutAAsVector.size(); i++) {
            collapsedNMutAAs += std::to_string(nMutAAsVector[i]) + ",";
        }
        if (!collapsedNMutAAs.empty()) {
            collapsedNMutAAs.pop_back();
        }
        
        dfName[j] = (*mutantSummaryIt).first;
        dfSeq[j] = collapsedSequence;
        dfReads[j] = (*mutantSummaryIt).second.nReads;
        dfMaxReads[j] = (*mutantSummaryIt).second.maxNReads;
        dfUmis[j] = (*mutantSummaryIt).second.umi.size();
        dfMutBases[j] = collapsedNMutBases;
        dfMutCodons[j] = collapsedNMutCodons;
        dfMutAAs[j] = collapsedNMutAAs;
        dfMutantNameBase[j] = collapsedMutantNameBase;
        dfMutantNameCodon[j] = collapsedMutantNameCodon;
        dfMutantNameBaseHGVS[j] = collapsedMutantNameBaseHGVS;
        dfMutantNameAA[j] = collapsedMutantNameAA;
        dfMutantNameAAHGVS[j] = collapsedMutantNameAAHGVS;
        dfMutationTypes[j] = collapsedMutationTypes;
        dfSeqAA[j] = collapsedSequenceAA;
        dfVarLengths[j] = (*mutantSummaryIt).second.varLengths;
        j++;
    }
    
    DataFrame df = DataFrame::create(Named("mutantName") = dfName,
                                     Named("sequence") = dfSeq,
                                     Named("nbrReads") = dfReads,
                                     Named("maxNbrReads") = dfMaxReads,
                                     Named("nbrUmis") = dfUmis,
                                     Named("nbrMutBases") = dfMutBases,
                                     Named("nbrMutCodons") = dfMutCodons,
                                     Named("nbrMutAAs") = dfMutAAs,
                                     Named("varLengths") = dfVarLengths,
                                     Named("mutantNameBase") = dfMutantNameBase,
                                     Named("mutantNameCodon") = dfMutantNameCodon,
                                     Named("mutantNameBaseHGVS") = dfMutantNameBaseHGVS,
                                     Named("mutantNameAA") = dfMutantNameAA,
                                     Named("mutantNameAAHGVS") = dfMutantNameAAHGVS,
                                     Named("mutationTypes") = dfMutationTypes,
                                     Named("sequenceAA") = dfSeqAA,
                                     Named("stringsAsFactors") = false);
    
    return df;
}

// Collapse similar sequences
//   `seqs`: unique valid nucleotide sequences (or pairs of sequences  concatenated with "_") of equal lengths
//   `scores`: scores corresponding to `seqs`
// [[Rcpp::export]]
Rcpp::DataFrame groupSimilarSequences(std::vector<std::string> seqs,
                                      std::vector<double> scores, 
                                      double variableCollapseMaxDist, 
                                      int variableCollapseMinReads,
                                      double variableCollapseMinRatio,
                                      bool verbose) {
    
    // combine seqs and scores
    std::map<std::string, double> seqsScores;
    std::map<std::string, double>::iterator seqsScoresIt, seqsScoresSimIt;
    for (size_t i = 0; i < seqs.size(); i++) {
        seqsScores.insert(std::pair<std::string, double>(seqs[i], scores[i]));
    }
    
    // get sequence length
    size_t seqlen = seqs[0].length();

    // calculate Hamming distance tolerance
    int tol;
    if (variableCollapseMaxDist >= 1.0) {
        tol = (int)variableCollapseMaxDist;
    } else {
        tol = (int)(variableCollapseMaxDist *
            (seqs[0].find("_") != std::string::npos ? seqlen-1 : seqlen));
    }
    
    if (verbose) {
        Rcout << "start collapsing variable sequences (tolerance: " << tol << ")...";
    }
    
    // sort seqsScores decreasingly by scores
    // ... create an empty intermediate vector
    std::vector<std::pair<std::string, double>> vec;
    std::vector<std::pair<std::string, double>>::iterator vecIt;
    // copy key-value pairs from seqsScores to vec
    std::copy(seqsScores.begin(), seqsScores.end(),
              std::back_inserter<std::vector<std::pair<std::string, double>>>(vec));
    // ... sort vec by decreasing order of pair.second
    //     (if second values are equal, order by the pair's first value)
    std::sort(vec.begin(), vec.end(),
              [](const std::pair<std::string, double>& l,
                 const std::pair<std::string, double>& r) {
                  if (l.second != r.second)
                      return l.second > r.second;
                  return l.first < r.first;
              });
    
    // store sequences in BK tree
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
    
    std::string querySeq, representativeSeq;
    std::vector<std::string> simSeqs;
    std::map<std::string, std::string> single2representative;
    
    // start querying in the order of tree.items (ordered decreasingly by scores)
    size_t start_size = (double)tree.size;
    while (tree.size > 0) {
        querySeq = tree.first();
        // check in seqsScores if score for querySeq is < variableCollapseMinReads
        seqsScoresIt = seqsScores.find(querySeq);
        if (variableCollapseMinReads > 0 &&
            seqsScoresIt != seqsScores.end() &&
            (*seqsScoresIt).second < variableCollapseMinReads) {
            // in that case, score < variableCollapseMinReads for all other
            // sequences in the tree as well
            // if so, extract all remaining sequences in the tree and
            // add them to single2representative, each mapping to itself
            std::vector<std::string> rest = tree.get_all();
            for (size_t i = 0; i < rest.size(); i++) {
                single2representative[rest[i]] = rest[i];
            }
            tree.remove_all();
            // after that, tree.size = 0 so the loop will stop
        } else {
            simSeqs = tree.search(querySeq, tol);
            for (size_t i = 0; i < simSeqs.size(); i++) {
                // check that the score for querySeq is high enough compared to that of simSeqs[i]
                // if not, don't collapse simSeqs[i] with querySeq
                // must check explicitly if simSeqs[i] == querySeq, since the ratio in
                // that case will always be 1, and the function will loop indefinitely
                // if the querySeq is not removed
                if (((seqsScoresSimIt = seqsScores.find(simSeqs[i])) != seqsScores.end() &&
                    (*seqsScoresIt).second >= variableCollapseMinRatio * (*seqsScoresSimIt).second) ||
                    querySeq == simSeqs[i]) {
                    single2representative[simSeqs[i]] = querySeq;
                    tree.remove(simSeqs[i]);
                }
            }
            
            // check for user interruption and print progress
            if ((start_size - tree.size) % 2000 == 0) { // every 2,000 queries (every ~1-2s)
                Rcpp::checkUserInterrupt(); // ... check for user interrupt
                // ... and give an update
                if (verbose && (start_size - tree.size) % 2000 == 0) {
                    Rcout << "    " << std::setprecision(4) <<
                        (100.0 * ((double)(start_size - tree.size) / (double)start_size)) <<
                            "% done" << std::endl;
                }
            }
        }
    }
    
    // return grouping results as a DataFrame
    std::vector<std::string> repseqs(seqs.size());
    std::unordered_set<std::string> uniqueRepseqs;
    for (size_t i = 0; i < seqs.size(); i++) {
        repseqs[i] = single2representative[seqs[i]];
        uniqueRepseqs.insert(repseqs[i]);
    }
    Rcpp::DataFrame df = DataFrame::create(Named("sequence") = seqs,
                                           Named("representative") = repseqs);

    if (verbose) {
        Rcout << "done (reduced from " << seqs.size() << " to " << uniqueRepseqs.size() << ")" << std::endl;
    }
    uniqueRepseqs.clear();

    return df;
}



// [[Rcpp::export]]
List digestFastqsCpp(std::vector<std::string> fastqForwardVect,
                     std::vector<std::string> fastqReverseVect,
                     bool mergeForwardReverse,
                     size_t minOverlap, size_t maxOverlap,
                     size_t minMergedLength, size_t maxMergedLength,
                     double maxFracMismatchOverlap, bool greedyOverlap,
                     bool revComplForward, bool revComplReverse,
                     std::string elementsForward,
                     std::vector<int> elementLengthsForward,
                     std::string elementsReverse,
                     std::vector<int> elementLengthsReverse,
                     std::string adapterForward, std::string adapterReverse,
                     std::vector<std::string> primerForward,
                     std::vector<std::string> primerReverse,
                     std::vector<std::string> wildTypeForward,
                     std::vector<std::string> wildTypeForwardNames,
                     std::vector<std::string> wildTypeReverse,
                     std::vector<std::string> wildTypeReverseNames,
                     std::vector<std::string> constantForward,
                     std::vector<std::string> constantReverse,
                     double avePhredMinForward = 20.0, double avePhredMinReverse = 20.0,
                     int variableNMaxForward = 0, int variableNMaxReverse = 0,
                     int umiNMax = 0,
                     int nbrMutatedCodonsMaxForward = 1,
                     int nbrMutatedCodonsMaxReverse = 1,
                     int nbrMutatedBasesMaxForward = -1,
                     int nbrMutatedBasesMaxReverse = -1,
                     CharacterVector forbiddenMutatedCodonsForward = "NNW",
                     CharacterVector forbiddenMutatedCodonsReverse = "NNW",
                     bool useTreeWTmatch = false,
                     bool collapseToWTForward = false,
                     bool collapseToWTReverse = false,
                     double mutatedPhredMinForward = 0.0,
                     double mutatedPhredMinReverse = 0.0,
                     std::string mutNameDelimiter = ".",
                     int constantMaxDistForward = -1,
                     int constantMaxDistReverse = -1,
                     double umiCollapseMaxDist = 0.0,
                     std::string filteredReadsFastqForward = "",
                     std::string filteredReadsFastqReverse = "",
                     int maxNReads = -1, bool verbose = false,
                     int nThreads = 1, int chunkSize = 100000,
                     size_t maxReadLength = 1024) {

  // See https://github.com/Rdatatable/data.table/issues/3300#issuecomment-457017735 for
  // a discussion on limiting the number of threads available to OpenMP
#ifdef _OPENMP
  if (nThreads > std::min(omp_get_max_threads(), omp_get_thread_limit())) {
    warning("'nThreads' is larger than the number of available threads. Reducing 'nThreads' to %i", std::min(omp_get_max_threads(), omp_get_thread_limit()));
    nThreads = std::min(omp_get_max_threads(), omp_get_thread_limit());
  }
#else
  if (nThreads > 1) {
    warning("OpenMP parallelization not available. Ignoring 'nThreads'.");
  }
#endif

  // Biostrings::IUPAC_CODE_MAP
  std::map<char,std::vector<char>> IUPAC = initializeIUPAC();
  std::map<char,std::string> threeAA = initializeThreeAA();

  // --------------------------------------------------------------------------
  // declare variables
  // --------------------------------------------------------------------------
  size_t BUFFER_SIZE = maxReadLength + 2;
  FastqBuffer *chunkBuffer = new FastqBuffer((size_t)chunkSize, BUFFER_SIZE, fastqReverseVect[0].compare("") != 0);

  int nTot = 0, nAdapter = 0, nNoPrimer = 0, nReadWrongLength = 0, nTooManyMutConstant = 0;
  int nNoValidOverlap = 0, nAvgVarQualTooLow = 0, nTooManyNinVar = 0, nTooManyNinUMI = 0;
  int nTooManyMutCodons = 0, nTooManyMutBases = 0, nForbiddenCodons = 0, nMutQualTooLow = 0;
  int nTooManyBestWTHits = 0, nTooManyBestConstantHits = 0, nRetain = 0;

  std::map<std::string, mutantInfo> mutantSummary;
  std::map<std::string, mutantInfo>::iterator mutantSummaryIt, mutantSummarySimIt;
  std::vector<int> nPhredCorrectForward(100, 0), nPhredMismatchForward(100, 0);
  std::vector<int> nPhredCorrectReverse(100, 0), nPhredMismatchReverse(100, 0);

  // enumerate forbidden codons based on forbiddenMutatedCodons
  std::set<std::string> forbiddenCodonsForward = enumerateCodonsFromIUPAC(forbiddenMutatedCodonsForward, IUPAC, verbose);
  std::set<std::string> forbiddenCodonsReverse = enumerateCodonsFromIUPAC(forbiddenMutatedCodonsReverse, IUPAC, verbose);

  // max number of mismatches for constant sequences
  size_t upperBoundMismatchConstForward = 0, upperBoundMismatchConstReverse = 0;

  // build BKtree(s) for wildtype sequences
  bool useTreeWTmatchForward = useTreeWTmatch, useTreeWTmatchReverse = useTreeWTmatch;
  BKtree wtTreeForward, wtTreeReverse;
  std::map<std::string, int> wtTreeForwardIdx, wtTreeReverseIdx;
  size_t upperBoundMismatchForward = 0, upperBoundMismatchReverse = 0;
  // forward
  if (wildTypeForward[0].compare("") != 0 && useTreeWTmatch) {
    std::vector<std::string>::iterator wtTreeForwardVecIt;
    size_t wtForwardLen = wildTypeForward[0].length();
    for (wtTreeForwardVecIt = wildTypeForward.begin();
         wtTreeForwardVecIt != wildTypeForward.end(); wtTreeForwardVecIt++) {
      if ((*wtTreeForwardVecIt).length() != wtForwardLen) {
        // not all WT sequences are of the same length ->
        // don't use the tree for finding the most similar WT sequence
        wtTreeForward.remove_all();
        useTreeWTmatchForward = false;
        break;
      } else {
        wtTreeForward.insert(*wtTreeForwardVecIt);
      }
    }
    for (size_t i = 0; i < wildTypeForward.size(); i++) {
      wtTreeForwardIdx[wildTypeForward[i]] = i;
    }
  }
  // determine the upper bound for the number of base mismatches in the tree search
  if (nbrMutatedCodonsMaxForward != (-1)) {
    upperBoundMismatchForward = 3 * nbrMutatedCodonsMaxForward;
  } else {
    upperBoundMismatchForward = nbrMutatedBasesMaxForward;
  }

  if (wildTypeReverse[0].compare("") != 0 && useTreeWTmatch) {
    std::vector<std::string>::iterator wtTreeReverseVecIt;
    size_t wtReverseLen = wildTypeReverse[0].length();
    for (wtTreeReverseVecIt = wildTypeReverse.begin();
         wtTreeReverseVecIt != wildTypeReverse.end(); wtTreeReverseVecIt++) {
      if ((*wtTreeReverseVecIt).length() != wtReverseLen) {
        // not all WT sequences are of the same length ->
        // don't use the tree for finding the most similar WT sequence
        wtTreeReverse.remove_all();
        useTreeWTmatchReverse = false;
        break;
      } else {
        wtTreeReverse.insert(*wtTreeReverseVecIt);
      }
    }
    for (size_t i = 0; i < wildTypeReverse.size(); i++) {
      wtTreeReverseIdx[wildTypeReverse[i]] = i;
    }
  }
  // determine the upper bound for the number of base mismatches in the tree search
  if (nbrMutatedCodonsMaxForward != (-1)) {
    upperBoundMismatchReverse = 3 * nbrMutatedCodonsMaxReverse;
  } else {
    upperBoundMismatchReverse = nbrMutatedBasesMaxReverse;
  }

  // open files for outputting filtered reads
  gzFile outfile1 = NULL, outfile2 = NULL;
  if (filteredReadsFastqForward.compare("") != 0) {
    outfile1 = openFastq(filteredReadsFastqForward, "wb");
    if (filteredReadsFastqReverse.compare("") != 0) {
      outfile2 = openFastq(filteredReadsFastqReverse, "wb");
    }
  }

  // --------------------------------------------------------------------------
  // iterate over fastq files
  // --------------------------------------------------------------------------
  for (size_t f = 0; f < fastqForwardVect.size(); f++) {
    // if maxNReads has been reached in the previous file, break
    bool done;
    if (maxNReads != (-1) && nTot >= maxNReads) {
      done = true;
    } else {
      done = false;
    }
    std::string fastqForward = fastqForwardVect[f];
    std::string fastqReverse = fastqReverseVect[f];

    // --------------------------------------------------------------------------
    // open fastq file(s)
    // --------------------------------------------------------------------------
    gzFile file1 = openFastq(fastqForward, "rb");
    gzFile file2 = NULL;
    if (fastqReverse.compare("") != 0) {
      file2 = openFastq(fastqReverse, "rb");
    }

    // iterate over sequences
    if (verbose) {
      Rcout << "start reading sequences for file " <<
               (fastqReverse.compare("") != 0 ? "pair " : "") << (f + 1) <<
               " of " << fastqForwardVect.size() << "..." << std::endl;
    }

    // counter for reads being added to the chunk vector
    size_t iChunk = 0;
    while (done == false) {
      done = get_next_seq(file1,
                          chunkBuffer->seq1 + iChunk * BUFFER_SIZE,
                          chunkBuffer->qual1 + iChunk * BUFFER_SIZE,
                          BUFFER_SIZE);
      if (fastqReverse.compare("") != 0) {
        done = (done || get_next_seq(file2,
                                     chunkBuffer->seq2 + iChunk * BUFFER_SIZE,
                                     chunkBuffer->qual2 + iChunk * BUFFER_SIZE,
                                     BUFFER_SIZE));
      }
      if (done == false) {
        iChunk++;
        nTot++;
      }
      if (nTot % 200000 == 0) { // every 200,000 reads (every ~1.6 seconds)
        Rcpp::checkUserInterrupt(); // ... check for user interrupt
      }

      // if maxNReads has been reached, break
      if (maxNReads != (-1) && nTot >= maxNReads) {
        done = true;
      }

      // process reads in the chunk
      if ((int)iChunk == chunkSize || done) {
#ifdef _OPENMP
#pragma omp parallel for schedule(guided) num_threads(nThreads)
#endif
        for (size_t ci = 0; ci < iChunk; ci++) { // for each read in chunk
          // while (done == false) {
          // initialize read parts as empty strings (since they will be added onto
          // when extracting the respective parts from the reads)
          std::string varSeqForward = "", varSeqReverse = "", varQualForward = "";
          std::string varQualReverse = "", umiSeq = "", constSeqForward = "";
          std::string constSeqReverse = "", constQualForward = "", constQualReverse = "";
          std::string mutantName = "", mutantNameBase = "", mutantNameCodon = "",
              mutantNameBaseHGVS = "", mutantNameAA = "", mutantNameAAHGVS = "";
          std::vector<int> varLengthsForward, varLengthsReverse;
          std::string varLengthsForwardStr = "", varLengthsReverseStr = "";
          int nMutBases = 0, nMutCodons = 0, nMutAAs = 0;
          std::set<std::string> mutationTypes;
          int maxSim = 0;
          int constantLengthForward, constantLengthReverse;
          std::map<std::string, mutantInfo>::iterator mutantSummaryParIt;

          // convert C char* to C++ string
          std::string sseq1(chunkBuffer->seq1 + ci * BUFFER_SIZE);
          std::string squal1(chunkBuffer->qual1 + ci * BUFFER_SIZE);
          std::string sseq2, squal2;
          if (fastqReverse.compare("") != 0) {
            sseq2 = chunkBuffer->seq2 + ci * BUFFER_SIZE;
            squal2 = chunkBuffer->qual2 + ci * BUFFER_SIZE;
          }

          // search for adapter sequences and filter read pairs
          if ((adapterForward.compare("") != 0 && sseq1.find(adapterForward) != std::string::npos) ||
              (fastqReverse.compare("") != 0 && adapterReverse.compare("") != 0 && sseq2.find(adapterReverse) != std::string::npos)) {
#ifdef _OPENMP
            #pragma omp atomic
#endif
            nAdapter++;
            chunkBuffer->write_seq(ci, outfile1, outfile2, nTot-(int)iChunk+(int)ci, "adapter");
            continue;
          }

          // check if the last character(s) are new line, if so remove them
          removeEOL(sseq1);
          removeEOL(squal1);
          if (fastqReverse.compare("") != 0) {
            removeEOL(sseq2);
            removeEOL(squal2);
          }

          // Extract the components of the reads
          // UMI, constant sequence, variable sequence, forward/reverse
          // if processing of either forward or reverse read signals 'break', go to the next read pair
          if (!decomposeRead(sseq1, squal1, elementsForward, elementLengthsForward,
                             primerForward, umiSeq, varSeqForward, varQualForward,
                             varLengthsForward, constSeqForward, constQualForward,
                             nNoPrimer, nReadWrongLength)) {
            chunkBuffer->write_seq(ci, outfile1, outfile2, nTot-(int)iChunk+(int)ci, "noPrimer_readWrongLength");
            continue;
          }
          if (fastqReverse.compare("") != 0 &&
              !decomposeRead(sseq2, squal2, elementsReverse, elementLengthsReverse,
                             primerReverse, umiSeq, varSeqReverse, varQualReverse,
                             varLengthsReverse, constSeqReverse, constQualReverse,
                             nNoPrimer, nReadWrongLength)) {
            chunkBuffer->write_seq(ci, outfile1, outfile2, nTot-(int)iChunk+(int)ci, "noPrimer_readWrongLength");
            continue;
          }

          // reverse complement if requested
          if (revComplForward) {
            transform(
              begin(varSeqForward), end(varSeqForward),
              begin(varSeqForward), complement);
            reverse(begin(varSeqForward), end(varSeqForward));
            reverse(begin(varQualForward), end(varQualForward));
            std::reverse(varLengthsForward.begin(), varLengthsForward.end());
          }
          if (fastqReverse.compare("") != 0 && revComplReverse) {
            transform(
              begin(varSeqReverse), end(varSeqReverse),
              begin(varSeqReverse), complement);
            reverse(begin(varSeqReverse), end(varSeqReverse));
            reverse(begin(varQualReverse), end(varQualReverse));
            std::reverse(varLengthsReverse.begin(), varLengthsReverse.end());
          }

          // convert qualities to int
          std::vector<int> varIntQualForward(varSeqForward.length(), 0);
          for (size_t i = 0; i < varSeqForward.length(); i++) {
            varIntQualForward[i] = int(varQualForward[i]) - QUALITY_OFFSET;
          }

          std::vector<int> varIntQualReverse;
          if (fastqReverse.compare("") != 0) {
            varIntQualReverse.resize(varSeqReverse.length()); // set number of elements
            std::fill(varIntQualReverse.begin(), varIntQualReverse.end(), 0); // fill with zero
            for (size_t i = 0; i < varSeqReverse.length(); i++) {
              varIntQualReverse[i] = int(varQualReverse[i]) - QUALITY_OFFSET;
            }
          }

          // if requested, fuse forward and reverse reads
          if (mergeForwardReverse) {
            if (mergeReadPairPartial(varSeqForward, varIntQualForward,
                                     varSeqReverse, varIntQualReverse,
                                     varLengthsForward, varLengthsReverse,
                                     minOverlap, maxOverlap,
                                     minMergedLength, maxMergedLength,
                                     maxFracMismatchOverlap, greedyOverlap)) {
              // read should be filtered out - no valid overlap found
#ifdef _OPENMP
              #pragma omp atomic
#endif
              nNoValidOverlap++;
              chunkBuffer->write_seq(ci, outfile1, outfile2, nTot-(int)iChunk+(int)ci, "noValidOverlap");
              continue;
            }
          }
          // if no reverse sequence was provided, hereafter it is identical to the situation
          // where forward and reverse reads were merged
          bool noReverse = mergeForwardReverse || fastqReverse.compare("") == 0;

          // filter if the average quality in variable region is too low
          if (std::accumulate(varIntQualForward.begin(), varIntQualForward.end(), 0.0) <
            avePhredMinForward * varSeqForward.length() ||
            (!noReverse && std::accumulate(varIntQualReverse.begin(), varIntQualReverse.end(), 0.0) <
              avePhredMinReverse * varSeqReverse.length())) {
#ifdef _OPENMP
            #pragma omp atomic
#endif
            nAvgVarQualTooLow++;
            chunkBuffer->write_seq(ci, outfile1, outfile2, nTot-(int)iChunk+(int)ci, "avgVarQualTooLow");
            continue;
          }

          // filter if there are too many N's in variable regions
          if (std::count(varSeqForward.begin(), varSeqForward.end(), 'N') > variableNMaxForward ||
              (!noReverse && std::count(varSeqReverse.begin(), varSeqReverse.end(), 'N') > variableNMaxReverse)) {
#ifdef _OPENMP
            #pragma omp atomic
#endif
            nTooManyNinVar++;
            chunkBuffer->write_seq(ci, outfile1, outfile2, nTot-(int)iChunk+(int)ci, "tooManyNinVar");
            continue;
          }

          if (umiSeq != "" && std::count(umiSeq.begin(), umiSeq.end(), 'N') > umiNMax) {
#ifdef _OPENMP
            #pragma omp atomic
#endif
            nTooManyNinUMI++;
            chunkBuffer->write_seq(ci, outfile1, outfile2, nTot-(int)iChunk+(int)ci, "tooManyNinUMI");
            continue;
          }

          // if wildTypeForward is available...
          if (wildTypeForward[0].compare("") != 0) {
            int idxForward;
            if (useTreeWTmatchForward) {
              idxForward = findClosestRefSeqTree(varSeqForward, wtTreeForward,
                                                 wtTreeForwardIdx, upperBoundMismatchForward, maxSim);
            } else {
              idxForward = findClosestRefSeqEarlyStop(varSeqForward, wildTypeForward,
                                                      upperBoundMismatchForward, maxSim);
            }
            // if idxForward = NO_SIMILAR_REF (defined above), no similar enough wildtype sequence was found - skip the read
            if (idxForward == NO_SIMILAR_REF) {
              if (nbrMutatedCodonsMaxForward != (-1)) {
#ifdef _OPENMP
                #pragma omp atomic
#endif
                nTooManyMutCodons++;
              } else {
#ifdef _OPENMP
                #pragma omp atomic
#endif
                nTooManyMutBases++;
              }
              chunkBuffer->write_seq(ci, outfile1, outfile2, nTot-(int)iChunk+(int)ci, "mutQualTooLow_tooManyMutCodons_forbiddenCodons");
              continue;
            }
            // if idxForward = TOO_MANY_BEST_REF (defined above), too many equally good hits among the WT sequences - skip the read
            if (idxForward == TOO_MANY_BEST_REF) {
#ifdef _OPENMP
              #pragma omp atomic
#endif
              nTooManyBestWTHits++;
              chunkBuffer->write_seq(ci, outfile1, outfile2, nTot-(int)iChunk+(int)ci, "tooManyBestWTHits");
              continue;
            }
            std::string wtForward = wildTypeForward[idxForward];
            std::string wtNameForward = wildTypeForwardNames[idxForward];
            
            // if the variable region is longer than the (best matching) WT sequence, 
            // skip the read (consider it to have the wrong length)
            if (varSeqForward.length() > wtForward.length()) {
#ifdef _OPENMP
              #pragma omp atomic
#endif
              nReadWrongLength++;
              chunkBuffer->write_seq(ci, outfile1, outfile2, nTot-(int)iChunk+(int)ci, "noPrimer_readWrongLength");
              continue;  
            }
            
            if (compareToWildtype(varSeqForward, wtForward, varIntQualForward,
                                  mutatedPhredMinForward, nbrMutatedCodonsMaxForward, forbiddenCodonsForward,
                                  wtNameForward, nbrMutatedBasesMaxForward, nMutQualTooLow,
                                  nTooManyMutCodons, nForbiddenCodons, nTooManyMutBases, 
                                  mutantName, mutantNameBase, mutantNameCodon,
                                  mutantNameBaseHGVS, mutantNameAA, mutantNameAAHGVS,
                                  nMutBases, nMutCodons, 
                                  nMutAAs, mutationTypes, mutNameDelimiter, collapseToWTForward,
                                  threeAA)) {
              // read is to be filtered out
              chunkBuffer->write_seq(ci, outfile1, outfile2, nTot-(int)iChunk+(int)ci, "mutQualTooLow_tooManyMutCodons_forbiddenCodons");
              continue;
            }
          } else if (varSeqForward.length() > 0) { // variable seq, but no reference -> add variable seq to mutantName
            mutantName += (varSeqForward + std::string("_"));
            mutantNameBase += (varSeqForward + std::string("_"));
            mutantNameCodon += (varSeqForward + std::string("_"));
            mutantNameAA += (translateString(varSeqForward) + std::string("_"));
            // don't set mutantNameBaseHGVS, mutantNameAAHGVS if there is no wt
          }
          // convert varLengthsForward to string
          if (varLengthsForward.size() > 0) {
            for (size_t i = 0; i < varLengthsForward.size(); i++) {
                varLengthsForwardStr += (std::to_string(varLengthsForward[i]) + ",");
            }
            varLengthsForwardStr.pop_back();
          }

          // if wildTypeReverse is available...
          if (!noReverse && wildTypeReverse[0].compare("") != 0) {
            int idxReverse;
            if (useTreeWTmatchReverse) {
              idxReverse = findClosestRefSeqTree(varSeqReverse, wtTreeReverse,
                                                 wtTreeReverseIdx, upperBoundMismatchReverse, maxSim);
            } else {
              idxReverse = findClosestRefSeqEarlyStop(varSeqReverse, wildTypeReverse,
                                                      upperBoundMismatchReverse, maxSim);
            }
            // if idxReverse = NO_SIMILAR_REF (defined above), no similar enough wildtype sequence was found - skip the read
            if (idxReverse == NO_SIMILAR_REF) {
              if (nbrMutatedCodonsMaxReverse != (-1)) {
#ifdef _OPENMP
                #pragma omp atomic
#endif
                nTooManyMutCodons++;
              } else {
#ifdef _OPENMP
                #pragma omp atomic
#endif
                nTooManyMutBases++;
              }
              chunkBuffer->write_seq(ci, outfile1, outfile2, nTot-(int)iChunk+(int)ci, "mutQualTooLow_tooManyMutCodons_forbiddenCodons");
              continue;
            }
            // if idxReverse = TOO_MANY_BEST_REF (defined above), too many equally good hits among the WT sequences - skip the read
            if (idxReverse == TOO_MANY_BEST_REF) {
#ifdef _OPENMP
              #pragma omp atomic
#endif
              nTooManyBestWTHits++;
              chunkBuffer->write_seq(ci, outfile1, outfile2, nTot-(int)iChunk+(int)ci, "tooManyBestWTHits");
              continue;
            }
            std::string wtReverse = wildTypeReverse[idxReverse];
            std::string wtNameReverse = wildTypeReverseNames[idxReverse];
            
            // if the variable region is longer than the (best matching) WT sequence, 
            // skip the read (consider it to have the wrong length)
            if (varSeqReverse.length() > wtReverse.length()) {
#ifdef _OPENMP
              #pragma omp atomic
#endif
              nReadWrongLength++;
              chunkBuffer->write_seq(ci, outfile1, outfile2, nTot-(int)iChunk+(int)ci, "noPrimer_readWrongLength");
              continue;  
            }
            
            if (compareToWildtype(varSeqReverse, wtReverse, varIntQualReverse,
                                  mutatedPhredMinReverse, nbrMutatedCodonsMaxReverse, forbiddenCodonsReverse,
                                  wtNameReverse, nbrMutatedBasesMaxReverse, nMutQualTooLow,
                                  nTooManyMutCodons, nForbiddenCodons, nTooManyMutBases, 
                                  mutantName, mutantNameBase, mutantNameCodon,
                                  mutantNameBaseHGVS, mutantNameAA, mutantNameAAHGVS,
                                  nMutBases, nMutCodons, 
                                  nMutAAs, mutationTypes, mutNameDelimiter, collapseToWTReverse,
                                  threeAA)) {
              // read is to be filtered out
              chunkBuffer->write_seq(ci, outfile1, outfile2, nTot-(int)iChunk+(int)ci, "mutQualTooLow_tooManyMutCodons_forbiddenCodons");
              continue;
            }
          } else if (!noReverse && varSeqReverse.length() > 0) { // variable seq, but no reference -> add variable seq to mutantName
            mutantName += (varSeqReverse + std::string("_"));
            mutantNameBase += (varSeqReverse + std::string("_"));
            mutantNameCodon += (varSeqReverse + std::string("_"));
            mutantNameAA += (translateString(varSeqReverse) + std::string("_"));
            // don't set mutantNameBaseHGVS, mutantNameAAHGVS if there is no wt
          }
          // convert varLengthsReverse to string
          if (varLengthsReverse.size() > 0) {
              for (size_t i = 0; i < varLengthsReverse.size(); i++) {
                  varLengthsReverseStr += (std::to_string(varLengthsReverse[i]) + ",");
              }
              varLengthsReverseStr.pop_back();
          }

          // for retained reads, count numbers of (mis-)matching bases in constant seq by Phred quality
          constantLengthForward = (int) constSeqForward.length();
          constantLengthReverse = (int) constSeqReverse.length();
          std::vector<int> constIntQualForward(constantLengthForward, 0);
          std::vector<int> constIntQualReverse(constantLengthReverse,0);
          int idxConstForward = 0;
          int idxConstReverse = 0;

          if (constantForward[0].compare("") != 0) {
            // reverse complement if requested
            if (revComplForward) {
              transform(begin(constSeqForward), end(constSeqForward),
                        begin(constSeqForward), complement);
              reverse(constSeqForward.begin(), constSeqForward.end());
              reverse(constQualForward.begin(), constQualForward.end());
            }

            // populate an integer vector of base qualities
            for (size_t i = 0; i < (size_t) constantLengthForward; i++) {
              constIntQualForward[i] = int(constQualForward[i]) - QUALITY_OFFSET;
            }

            // find closest constant sequence and tabulate mismatches by quality
            maxSim = 0;
            // define upperBoundMismatch for this read
            if (constantMaxDistForward != (-1)) {
              upperBoundMismatchConstForward = (size_t)constantMaxDistForward;
            } else {
              upperBoundMismatchConstForward = constSeqForward.size();
            }
            idxConstForward = findClosestRefSeqEarlyStop(constSeqForward, constantForward,
                                                         upperBoundMismatchConstForward, maxSim);
            if (idxConstForward == NO_SIMILAR_REF) {
#ifdef _OPENMP
              #pragma omp atomic
#endif
              nTooManyMutConstant++;
              chunkBuffer->write_seq(ci, outfile1, outfile2, nTot-(int)iChunk+(int)ci, "tooManyMutConstant");
              continue;
            }
            // more than one equally good best hit among the constant sequences - skip read
            if (idxConstForward == TOO_MANY_BEST_REF) {
#ifdef _OPENMP
              #pragma omp atomic
#endif
              nTooManyBestConstantHits++;
              chunkBuffer->write_seq(ci, outfile1, outfile2, nTot-(int)iChunk+(int)ci, "tooManyBestConstantHits");
              continue;
            }
          }

          if (fastqReverse.compare("") != 0 &&
              constantReverse[0].compare("") != 0) {
            // reverse (complement) sequence and quality string
            if (revComplReverse) {
              transform(begin(constSeqReverse), end(constSeqReverse),
                        begin(constSeqReverse), complement);
              reverse(constSeqReverse.begin(), constSeqReverse.end());
              reverse(constQualReverse.begin(), constQualReverse.end());
            }

            // populate an integer vector of base qualities
            for (size_t i = 0; i < (size_t) constantLengthReverse; i++) {
              constIntQualReverse[i] = int(constQualReverse[i]) - QUALITY_OFFSET;
            }

            // find closest constant sequence and tabulate mismatches by quality
            maxSim = 0;
            // define upperBoundMismatch for this read
            if (constantMaxDistReverse != (-1)) {
              upperBoundMismatchConstReverse = (size_t)constantMaxDistReverse;
            } else {
              upperBoundMismatchConstReverse = constSeqReverse.size();
            }
            idxConstReverse = findClosestRefSeqEarlyStop(constSeqReverse, constantReverse,
                                                         upperBoundMismatchConstReverse, maxSim);
            if (idxConstReverse == NO_SIMILAR_REF) {
#ifdef _OPENMP
              #pragma omp atomic
#endif
              nTooManyMutConstant++;
              chunkBuffer->write_seq(ci, outfile1, outfile2, nTot-(int)iChunk+(int)ci, "tooManyMutConstant");
              continue;
            }
            // more than one equally good best hit among the constant sequences - skip read
            if (idxConstReverse == TOO_MANY_BEST_REF) {
#ifdef _OPENMP
              #pragma omp atomic
#endif
              nTooManyBestConstantHits++;
              chunkBuffer->write_seq(ci, outfile1, outfile2, nTot-(int)iChunk+(int)ci, "tooManyBestConstantHits");
              continue;
            }
          }

          if (constantForward[0].compare("") != 0) {
            tabulateBasesByQual(constSeqForward, constantForward[idxConstForward],
                                constIntQualForward,
                                nPhredCorrectForward, nPhredMismatchForward);
          }
          if (fastqReverse.compare("") != 0 &&
              constantReverse[0].compare("") != 0) {
            tabulateBasesByQual(constSeqReverse, constantReverse[idxConstReverse],
                                constIntQualReverse,
                                nPhredCorrectReverse, nPhredMismatchReverse);
          }

          // store the read pair
#ifdef _OPENMP
          #pragma omp atomic
#endif
          nRetain++;
          // ... create final mutant name
          if (mutantName.length() > 0) { // we have a least one mutation, or sequence-based name
            mutantName.pop_back(); // remove '_' at the end
            mutantNameBase.pop_back(); // remove '_' at the end
            mutantNameCodon.pop_back(); // remove '_' at the end
          } else {
            // will we ever go in here?
            if (wildTypeForward[0].compare("") != 0 ||
                (!noReverse && wildTypeReverse[0].compare("") != 0)) {
              mutantName = "WT";
              mutantNameBase = "WT";
              mutantNameCodon = "WT";
            }
          }
          if (mutantNameBaseHGVS.length() > 0) { // we have a (closest) wildtype name
            mutantNameBaseHGVS.pop_back(); // remove '_' at the end
            mutantNameAAHGVS.pop_back(); // remove '_' at the end
          }
          if (mutantNameAA.length() > 0) { // we have a least one mutation, or sequence-based name
            mutantNameAA.pop_back(); // remove '_' at the end
          } else {
            // will we ever go in here?
            if (wildTypeForward[0].compare("") != 0 ||
                (!noReverse && wildTypeReverse[0].compare("") != 0)) {
              mutantNameAA = "WT";
            }
          }
          if (!noReverse) { // "trans" experiment
            varSeqForward += (std::string("_") + varSeqReverse);
            varLengthsForwardStr += (std::string("_") + varLengthsReverseStr);
          }
          // ... check if mutant already exists in mutantSummary
#ifdef _OPENMP
#pragma omp critical
#endif
{
          if ((mutantSummaryParIt = mutantSummary.find(mutantName)) != mutantSummary.end()) {
            // ... ... update existing mutantInfo
            (*mutantSummaryParIt).second.nReads++;
            (*mutantSummaryParIt).second.maxNReads++;
            if (umiSeq != "") {
              (*mutantSummaryParIt).second.umi.insert(umiSeq);
            }
            (*mutantSummaryParIt).second.sequence.insert(varSeqForward);
            (*mutantSummaryParIt).second.nMutBases.insert(nMutBases);
            (*mutantSummaryParIt).second.nMutCodons.insert(nMutCodons);
            (*mutantSummaryParIt).second.nMutAAs.insert(nMutAAs);
            (*mutantSummaryParIt).second.mutationTypes.insert(mutationTypes.begin(), mutationTypes.end());
            (*mutantSummaryParIt).second.mutantNameBase.insert(mutantNameBase);
            (*mutantSummaryParIt).second.mutantNameCodon.insert(mutantNameCodon);
            (*mutantSummaryParIt).second.mutantNameBaseHGVS.insert(mutantNameBaseHGVS);
            (*mutantSummaryParIt).second.mutantNameAA.insert(mutantNameAA);
            (*mutantSummaryParIt).second.mutantNameAAHGVS.insert(mutantNameAAHGVS);
            (*mutantSummaryParIt).second.sequenceAA.insert(translateString(varSeqForward));
          } else {
            // ... ... create mutantInfo instance for this mutant and add it to mutantSummary
            mutantInfo newMutant;
            newMutant.nReads = 1;
            newMutant.maxNReads = 1;
            newMutant.nMutBases.insert(nMutBases);
            newMutant.nMutCodons.insert(nMutCodons);
            newMutant.nMutAAs.insert(nMutAAs);
            if (umiSeq != "") {
              newMutant.umi.insert(umiSeq);
            }
            newMutant.sequence.insert(varSeqForward);
            newMutant.varLengths = varLengthsForwardStr;
            newMutant.mutationTypes.insert(mutationTypes.begin(), mutationTypes.end());
            newMutant.mutantNameBase.insert(mutantNameBase);
            newMutant.mutantNameCodon.insert(mutantNameCodon);
            newMutant.mutantNameBaseHGVS.insert(mutantNameBaseHGVS);
            newMutant.mutantNameAA.insert(mutantNameAA);
            newMutant.mutantNameAAHGVS.insert(mutantNameAAHGVS);
            newMutant.sequenceAA.insert(translateString(varSeqForward));
            mutantSummary.insert(std::pair<std::string,mutantInfo>(mutantName, newMutant));
          }
}
        } // iterate over individual sequence pairs
        iChunk = 0;

        // ... and give an update
        if (verbose) {
          Rcout << "    " << nTot << " read pairs processed ("
                << std::setprecision(3) << (100*((double)nRetain/nTot)) << "% retained)" << std::endl;
        }
      }
      // check if end-of-file was reached
      if (done) {
        break;
      }
    }

    // clean up
    gzclose(file1);
    if (fastqReverse.compare("") != 0) {
      gzclose(file2);
    }
  } // iterate over fastq files

  // close output files
  if (outfile1 != NULL) {
    gzclose(outfile1);
  }
  if (outfile2 != NULL) {
    gzclose(outfile2);
  }

  if (verbose) {
    Rcout << "done reading sequences" << std::endl;
    Rcout << "retained " << mutantSummary.size() << " unique features" << std::endl;
  }

  // reading of reads is done - don't need the chunkBuffer anymore
  delete chunkBuffer;

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
    int mutCounter = 0;
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
      mutCounter++;
      // check for user interruption and print progress
      if (mutCounter % 200 == 0) { // every 200 queries
        Rcpp::checkUserInterrupt(); // ... check for user interrupt
      }
    }

    if (verbose) {
      Rcout << "done" << std::endl;
    }
  }

  // return results
  // ... parameters
  // put back wildtype sequences and names in Rcpp::StringVector
  Rcpp::StringVector wildTypeForwardRcpp(wildTypeForward.size());
  wildTypeForwardRcpp = wildTypeForward;
  wildTypeForwardRcpp.attr("names") = wildTypeForwardNames;
  Rcpp::StringVector wildTypeReverseRcpp(wildTypeReverse.size());
  wildTypeReverseRcpp = wildTypeReverse;
  wildTypeReverseRcpp.attr("names") = wildTypeReverseNames;
  
  std::vector<std::string> forbiddenCodonsUsedForward(forbiddenCodonsForward.begin(), forbiddenCodonsForward.end());
  std::vector<std::string> forbiddenCodonsUsedReverse(forbiddenCodonsReverse.begin(), forbiddenCodonsReverse.end());
  List param;
  param.push_back(fastqForwardVect, "fastqForward");
  param.push_back(fastqReverseVect, "fastqReverse");
  param.push_back(mergeForwardReverse, "mergeForwardReverse");
  param.push_back(minOverlap, "minOverlap");
  param.push_back(maxOverlap, "maxOverlap");
  param.push_back(minMergedLength, "minMergedLength");
  param.push_back(maxMergedLength, "maxMergedLength");
  param.push_back(maxFracMismatchOverlap, "maxFracMismatchOverlap");
  param.push_back(greedyOverlap, "greedyOverlap");
  param.push_back(revComplForward, "revComplForward");
  param.push_back(revComplReverse, "revComplReverse");
  param.push_back(elementsForward, "elementsForward");
  param.push_back(elementLengthsForward, "elementLengthsForward");
  param.push_back(elementsReverse, "elementsReverse");
  param.push_back(elementLengthsReverse, "elementLengthsReverse");
  param.push_back(adapterForward, "adapterForward");
  param.push_back(adapterReverse, "adapterReverse");
  param.push_back(primerForward, "primerForward");
  param.push_back(primerReverse, "primerReverse");
  param.push_back(wildTypeForwardRcpp, "wildTypeForward");
  param.push_back(wildTypeReverseRcpp, "wildTypeReverse");
  param.push_back(constantForward, "constantForward");
  param.push_back(constantReverse, "constantReverse");
  param.push_back(avePhredMinForward, "avePhredMinForward");
  param.push_back(avePhredMinReverse, "avePhredMinReverse");
  param.push_back(variableNMaxForward, "variableNMaxForward");
  param.push_back(variableNMaxReverse, "variableNMaxReverse");
  param.push_back(umiNMax, "umiNMax");
  param.push_back(nbrMutatedCodonsMaxForward, "nbrMutatedCodonsMaxForward");
  param.push_back(nbrMutatedCodonsMaxReverse, "nbrMutatedCodonsMaxReverse");
  param.push_back(nbrMutatedBasesMaxForward, "nbrMutatedBasesMaxForward");
  param.push_back(nbrMutatedBasesMaxReverse, "nbrMutatedBasesMaxReverse");
  param.push_back(forbiddenCodonsUsedForward, "forbiddenMutatedCodonsForward");
  param.push_back(forbiddenCodonsUsedReverse, "forbiddenMutatedCodonsReverse");
  param.push_back(useTreeWTmatch, "useTreeWTmatch");
  param.push_back(collapseToWTForward, "collapseToWTForward");
  param.push_back(collapseToWTReverse, "collapseToWTReverse");
  param.push_back(mutatedPhredMinForward, "mutatedPhredMinForward");
  param.push_back(mutatedPhredMinReverse, "mutatedPhredMinReverse");
  param.push_back(mutNameDelimiter, "mutNameDelimiter");
  param.push_back(constantMaxDistForward, "constantMaxDistForward");
  param.push_back(constantMaxDistReverse, "constantMaxDistReverse");
  param.push_back(umiCollapseMaxDist, "umiCollapseMaxDist");
  param.push_back(filteredReadsFastqForward, "filteredReadsFastqForward");
  param.push_back(filteredReadsFastqReverse, "filteredReadsFastqReverse");
  param.push_back(maxNReads, "maxNReads");
  param.push_back(nThreads, "nThreads");
  param.push_back(chunkSize, "chunkSize");
  param.push_back(maxReadLength, "maxReadLength");

  // ... error statistics
  DataFrame err = DataFrame::create(Named("PhredQuality") = seq_len(100) - 1,
                                    Named("nbrMatchForward") = nPhredCorrectForward,
                                    Named("nbrMismatchForward") = nPhredMismatchForward,
                                    Named("nbrMatchReverse") = nPhredCorrectReverse,
                                    Named("nbrMismatchReverse") = nPhredMismatchReverse);
  
  // ... filter statistics
  DataFrame filt = DataFrame::create(Named("nbrTotal") = nTot,
                                     Named("f1_nbrAdapter") = nAdapter,
                                     Named("f2_nbrNoPrimer") = nNoPrimer,
                                     Named("f3_nbrReadWrongLength") = nReadWrongLength,
                                     Named("f4_nbrNoValidOverlap") = nNoValidOverlap,
                                     Named("f5_nbrAvgVarQualTooLow") = nAvgVarQualTooLow,
                                     Named("f6_nbrTooManyNinVar") = nTooManyNinVar,
                                     Named("f7_nbrTooManyNinUMI") = nTooManyNinUMI,
                                     Named("f8_nbrTooManyBestWTHits") = nTooManyBestWTHits,
                                     Named("f9_nbrMutQualTooLow") = nMutQualTooLow,
                                     Named("f10a_nbrTooManyMutCodons") = nTooManyMutCodons,
                                     Named("f10b_nbrTooManyMutBases") = nTooManyMutBases,
                                     Named("f11_nbrForbiddenCodons") = nForbiddenCodons,
                                     Named("f12_nbrTooManyMutConstant") = nTooManyMutConstant,
                                     Named("f13_nbrTooManyBestConstantHits") = nTooManyBestConstantHits,
                                     Named("nbrRetained") = nRetain);
  
  // ... main data frame
  DataFrame df = mutantSummaryToDataFrame(mutantSummary);
  
  // ... pack into list
  List L = List::create(Named("parameters") = param,
                        Named("filterSummary") = filt,
                        Named("summaryTable") = df,
                        Named("errorStatistics") = err);
  return L;
}
