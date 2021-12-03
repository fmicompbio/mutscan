#ifndef STRINGDIST_HPP
#define STRINGDIST_HPP

#include <Rcpp.h>
#include <string>

// prototypes for string distance functions, e.g. used by the BKtree class
int levenshtein_distance(const std::string &, const std::string &, int);
int hamming_distance(const std::string &, const std::string &, int);
int hamming_shift_distance(const std::string &, const std::string &, int);


#endif