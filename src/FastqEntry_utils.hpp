#include <Rcpp.h>
#include <string.h>

using namespace Rcpp;

class FastqEntry {
public:
  char *seq1, *qual1, *seq2, *qual2;
  
  // constructor
  FastqEntry() {
    seq1 = new char[BUFFER_SIZE];
    qual1 = new char[BUFFER_SIZE];
    seq2 = new char[BUFFER_SIZE];
    qual2 = new char[BUFFER_SIZE];
  }
  
  // copy constructor
  FastqEntry(const FastqEntry& fqe) {
    seq1 = new char[BUFFER_SIZE];
    qual1 = new char[BUFFER_SIZE];
    seq2 = new char[BUFFER_SIZE];
    qual2 = new char[BUFFER_SIZE];
    strncpy(seq1, fqe.seq1, BUFFER_SIZE);
    strncpy(qual1, fqe.qual1, BUFFER_SIZE);
    strncpy(seq2, fqe.seq2, BUFFER_SIZE);
    strncpy(qual2, fqe.qual2, BUFFER_SIZE);
  }
  
  // destructor
  ~FastqEntry() {
    delete[] seq1;
    delete[] qual1;
    delete[] seq2;
    delete[] qual2;
  }
  
};