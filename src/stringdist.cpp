#include "stringdist.hpp"

// distance metrices used by the BKtree class:
// calculate levenshtein distance between pair of strings
// [[Rcpp::export]]
int levenshtein_distance(const std::string &str1, const std::string &str2,
                         int ignored_variable = -1){
    size_t m = str1.length(), n = str2.length();
    
    int** dn = new int*[m + 1];
    for(size_t i = 0; i < m + 1; ++i) {
        dn[i] = new int[n + 1];
    }
    
    for (size_t i = 0; i < m + 1; ++i) {
        for (size_t j = 0; j < n + 1; ++j) {
            if (i == 0) {
                dn[0][j] = j; // deletions at the start of str2
                
            } else if (j == 0) {
                dn[i][0] = i; // deletions at the start of str1
                
            } else if (str1[i-1] == str2[j-1]) { // match of str1[i-1] and str2[j-1]
                dn[i][j] = dn[i-1][j-1];
                
            } else { // mismatch between str1[i-1] and str2[j-1] -> find minimal source
                dn[i][j] = 1 + std::min(
                    dn[i-1][j],          // deletion in str1
                           std::min(dn[i][j-1], // deletion in str2
                                    dn[i-1][j-1])        // mismatch
                );
            }
        }
    }
    
    int d = dn[m][n];
    for (size_t i = 0; i < m + 1; ++i) {
        delete [] dn[i];
    }
    delete [] dn;
    
    return d;
}

// calculate hamming distance between pair of strings of equal lengths
// int hamming_distance(const std::string &str1, const std::string &str2,
//                      int ignored_variable = -1){
//   int d = 0;
//   
//   for (size_t i = 0; i <= str1.length(); i++) {
//     d += (str1[i] != str2[i]);
//   }
//   
//   return d;
// }
// [[Rcpp::export]]
int hamming_distance(const std::string &str1, const std::string &str2,
                     int ignored_variable = -1){
    const char *a = str1.data(), *b = str2.data();
    
    return std::inner_product(a, a + str1.length(), b, 0,
                              std::plus<int>(), std::not_equal_to<int>());
}

// calculate hamming distance + shifts between pair of strings of equal lengths
// int hamming_shift_distance(const std::string &str1, const std::string &str2,
//                            int max_abs_shift = -1){
//   int d = str1.size(), ds = 0;
//   if (max_abs_shift < 0) {
//     max_abs_shift = (int)str1.size() - 1;
//   }
// 
//   d = str1.size();
//   for (int s = -max_abs_shift; s <= max_abs_shift; s++) {
//     ds = std::abs(s) + std::abs(s); // overhangs are treated as mismatches
//     for (size_t i = std::max(0, s), j = std::max(0, -s);
//          i < str1.length() + std::min(0, s); i++, j++) {
//       ds += (str1[i] != str2[j]);
//     }
//     if (ds < d) {
//       d = ds;
//     }
//   }
// 
//   return d;
// }
// [[Rcpp::export]]
int hamming_shift_distance(const std::string &str1, const std::string &str2,
                           int max_abs_shift = -1){
    int d = str1.size(), ds = 0;
    const char *a = str1.data(), *b = str2.data();
    if (max_abs_shift < 0) {
        max_abs_shift = (int)str1.size() - 1;
    }
    
    for (int s = -max_abs_shift; s <= max_abs_shift; s++) {
        ds = std::abs(s) + std::abs(s) + // overhangs are treated as mismatches
            std::inner_product(a + std::max(0, s),
                               a + str1.length() + std::min(0, s),
                               b + std::max(0, -s), 0,
                               std::plus<int>(), std::not_equal_to<int>());
        if (ds < d) {
            d = ds;
        }
    }
    
    return d;
}

