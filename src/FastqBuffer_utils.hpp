#include <Rcpp.h>
#include <string.h>
#include <zlib.h>

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
  
  // Write one (pair of) filtered reads to output fastq file(s)
  bool write_seq(gzFile file1, gzFile file2, const int n, const char* label) {
    bool success = true;
#ifdef _OPENMP
#pragma omp critical
#endif
{
  if (file1 != NULL) {
    std::string read_id = "@S" + std::to_string(n) + "_" + label + " 1\n";
    if (success && gzputs(file1, read_id.c_str()) == (-1)) {
      success = false;
    }
    if (success && gzputs(file1, seq1) == (-1)) {
      success = false;
    }
    if (success && gzputs(file1, "+\n") == (-1)) {
      success = false;
    }
    if (success && gzputs(file1, qual1) == (-1)) {
      success = false;
    }
  }
  
  if (file2 != NULL) {
    std::string read_id = "@S" + std::to_string(n) + "_" + label + " 2\n";
    if (success && gzputs(file2, read_id.c_str()) == (-1)) {
      success = false;
    }
    if (success && gzputs(file2, seq2) == (-1)) {
      success = false;
    }
    if (success && gzputs(file2, "+\n") == (-1)) {
      success = false;
    }
    if (success && gzputs(file2, qual2) == (-1)) {
      success = false;
    }
  }
}
return success;
  }
  
};
