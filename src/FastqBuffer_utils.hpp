#include <Rcpp.h>
#include <string.h>
#include <zlib.h>

using namespace Rcpp;

class FastqBuffer {
private:
  int nentries; // number of fastq entries allocated (BUFFER_SIZE-1 max length)
  bool paired;  // paired read entries?
  char *buffer; // single block of memory
  //  each entry e has 2-times (4-times for paired=true) BUFFER_SIZE bytes available
  //  entries e1, e2, ... are stored in buffer in the following order:
  //   e1.seq1, e2.seq1, ..., e1.qual1, e2.qual2, ..., e1.seq2, e2.seq2, ..., e1.qual2, e2.qual2, ...
  //   (space for seq2 and qual2 is only allocated for paired=true)
public:
  char *seq1, *qual1, *seq2, *qual2;
  
  // constructor
  FastqBuffer(int n, bool p = true) {
    nentries = n;
    paired = p;
    buffer = new char[nentries * BUFFER_SIZE * (paired ? 4 : 2)];
    seq1 = buffer;
    qual1 = buffer + (nentries * BUFFER_SIZE);
    if (paired) {
      seq2 = buffer + (2 * nentries * BUFFER_SIZE);
      qual2 = buffer + (3 * nentries * BUFFER_SIZE);
    } else {
      seq2 = NULL;
      qual2 = NULL;
    }
  }
  
  // copy constructor
  FastqBuffer(const FastqBuffer& fqb) {
    paired = fqb.paired;
    nentries = fqb.nentries;
    buffer = new char[nentries * BUFFER_SIZE * (paired ? 4 : 2)];
    strncpy(buffer, fqb.buffer, nentries * BUFFER_SIZE * (paired ? 4 : 2));
    seq1 = buffer;
    qual1 = buffer + (nentries * BUFFER_SIZE);
    if (paired) {
      seq2 = buffer + (2 * nentries * BUFFER_SIZE);
      qual2 = buffer + (3 * nentries * BUFFER_SIZE);
    } else {
      seq2 = NULL;
      qual2 = NULL;
    }
  }
  
  // destructor
  ~FastqBuffer() {
    delete[] buffer;
  }
  
  // Write one (pair of) filtered reads to output fastq file(s)
  bool write_seq(int i, gzFile file1, gzFile file2, const int n, const char* label) {
    bool success = (i < nentries);
#ifdef _OPENMP
#pragma omp critical
#endif
{
    if (success && file1 != NULL) {
      std::string read_id = "@S" + std::to_string(n) + "_" + label + " 1\n";
      if (success && gzputs(file1, read_id.c_str()) == (-1)) {
        success = false;
      }
      if (success && gzputs(file1, seq1 + (i * BUFFER_SIZE)) == (-1)) {
        success = false;
      }
      if (success && gzputs(file1, "+\n") == (-1)) {
        success = false;
      }
      if (success && gzputs(file1, qual1 + (i * BUFFER_SIZE)) == (-1)) {
        success = false;
      }
    }
    
    if (success && file2 != NULL) {
      std::string read_id = "@S" + std::to_string(n) + "_" + label + " 2\n";
      if (success && gzputs(file2, read_id.c_str()) == (-1)) {
        success = false;
      }
      if (success && gzputs(file2, seq2 + (i * BUFFER_SIZE)) == (-1)) {
        success = false;
      }
      if (success && gzputs(file2, "+\n") == (-1)) {
        success = false;
      }
      if (success && gzputs(file2, qual2 + (i * BUFFER_SIZE)) == (-1)) {
        success = false;
      }
    }
}
    return success;
  }

};