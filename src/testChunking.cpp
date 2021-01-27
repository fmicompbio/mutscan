#include <Rcpp.h>
#include <cerrno>
#include <zlib.h>
#include <string>
#include <vector>
#include <numeric>
#include <map>
#include <algorithm>    // std::sort
#include <chrono>
#include <thread>

#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

#ifndef BUFFER_SIZE 
#define BUFFER_SIZE 4096
// define constants that are used below
#define NO_SIMILAR_REF    -1 // no similar enough wildtype sequence was found
#define TOO_MANY_BEST_REF -2 // too many equally good hits among the WT sequences
#endif

#include "FastqEntry_utils.hpp"

using namespace std::placeholders;
using namespace Rcpp;

bool reached_end_of_file(gzFile, char*);
bool get_next_seq(gzFile, char*, char*);
gzFile openFastq(std::string, const char*);
int findClosestRefSeqEarlyStop(std::string&, Rcpp::StringVector&, size_t, int&);
void removeEOL(std::string&);

// [[Rcpp::export]]
std::vector<int> testChunking(std::vector<std::string> fastqForwardVect,
                  std::vector<std::string> fastqReverseVect,
                  Rcpp::StringVector wildTypeForward, 
                  Rcpp::StringVector wildTypeReverse, 
                  size_t nbrMutatedBasesMaxForward,
                  size_t nbrMutatedBasesMaxReverse,
                  int maxNReads = -1, int nthreads = 1, int chunkSize = 1) {
  // Rcout << "starting" << std::endl << std::flush;

  // std::vector<std::string> wildTypeForwardStd(wildTypeForward.size());
  // for (size_t j = 0; j < wildTypeForward.size(); j++) {
  //   wildTypeForwardStd[j] = wildTypeForward(j);
  // }
  
  int nTot = 0;
  FastqEntry** chunkVector = new FastqEntry*[chunkSize];
  // Rcout << "initialized chunkVector" << std::endl << std::flush;
  const size_t upperBoundMismatchForward = nbrMutatedBasesMaxForward;
  const size_t upperBoundMismatchReverse = nbrMutatedBasesMaxReverse;
  std::vector<int> matchResults;
  size_t i;
  
  for (int i = 0; i < chunkSize; i++) {
    chunkVector[i] = new FastqEntry;
  }
  
  for (size_t f = 0; f < fastqForwardVect.size(); f++) {
    bool done = false;
    std::string fastqForward = fastqForwardVect[f];
    std::string fastqReverse = fastqReverseVect[f];
    // Rcout << "fastqForwardVect item " << f << std::endl << std::flush;
    // --------------------------------------------------------------------------
    // open fastq files
    // --------------------------------------------------------------------------
    gzFile file1 = openFastq(fastqForward, "rb");
    gzFile file2 = NULL;
    if (fastqReverse.compare("") != 0) {
      file2 = openFastq(fastqReverse, "rb");
    }
    // Rcout << "opened the fastq files" << std::endl << std::flush;
    // counter for reads being added to the chunk vector
    size_t iChunk = 0;
    while (done == false) {
      done = get_next_seq(file1, chunkVector[iChunk]->seq1, chunkVector[iChunk]->qual1);
      if (fastqReverse.compare("") != 0) {
        done = (done || get_next_seq(file2, chunkVector[iChunk]->seq2, chunkVector[iChunk]->qual2));
      }
      // Rcout << "adding read " << iChunk << " " << std::string(chunkVector[iChunk]->seq1) << std::endl << std::flush;
      iChunk++;
      // update counters
      nTot++;
      // if maxNReads has been reached, break
      if (maxNReads != (-1) && nTot >= maxNReads) {
        done = true;
      }
      
      // Rcout << "nTot: " << nTot << std::endl << std::flush;
      // process reads in the chunk
      if ((int)iChunk == chunkSize || done) {
        // Rcout << "starting to process a chunk" << std::endl << std::flush;
#ifdef _OPENMP
#pragma omp parallel for schedule(guided) num_threads(nthreads) private(i) firstprivate(wildTypeForward, chunkVector) shared(matchResults)
#endif
        for (i = 0; i < iChunk; i++) {
          int maxSim = 0;
          int idxForward = -1;
          // convert C char* to C++ string
          std::string varSeqForward(chunkVector[i]->seq1);
          std::string varQualForward(chunkVector[i]->qual1);
          std::string varSeqReverse, varQualReverse;
          if (fastqReverse.compare("") != 0) {
            varSeqReverse = chunkVector[i]->seq2;
            varQualReverse = chunkVector[i]->qual2;
          }
          
          // check if the last character(s) are new line, if so remove them
          removeEOL(varSeqForward);
          removeEOL(varQualForward);
          if (fastqReverse.compare("") != 0) {
            removeEOL(varSeqReverse);
            removeEOL(varQualReverse);
          }

          idxForward = findClosestRefSeqEarlyStop(varSeqForward, wildTypeForward,
                                                  upperBoundMismatchForward, maxSim);
#ifdef _OPENMP
          #pragma omp critical
#endif
          {
            matchResults.push_back(idxForward);
          }
        } // end processing for one read
        iChunk = 0; 
      }
      
      // check if end-of-file was reached
      if (done) {
        break;
      }
    } // end while loop - all reads are processed
  } // end for each input fastq file
  
  for (int i = 0; i < chunkSize; i++) {
    delete chunkVector[i];
  }
  
  delete[] chunkVector;
  
  return matchResults;
}
  
  
  
  
  
  