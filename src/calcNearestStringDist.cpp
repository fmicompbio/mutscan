#include <Rcpp.h>
#include <string>
#include <vector>
#include <climits>
#include "stringdist.hpp"

#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#endif

using namespace Rcpp;

//' Calculate distances to the nearest string
//' 
//' Given a character vector, calculate the distance for each element
//' to the nearest neighbor amongst all the other elements.
//' 
//' @param x A character vector.
//' 
//' @return An integer vector of the same length as \code{x}.
//' @export
// [[Rcpp::export]]
IntegerVector calcNearestStringDist(std::vector<std::string> x) {
    size_t i, j;
    int dist1;
    IntegerVector dists(x.size(), INT_MAX);
    
    for (i = 0; i < x.size() - 1; i++) {
        for (j = i + 1; j < x.size(); j++) {
            dist1 = hamming_distance(x[i], x[j], -1);
            if (dists[i] > dist1)
                dists[i] = dist1;
            if (dists[j] > dist1)
                dists[j] = dist1;
        }
    }
    
    return dists;
}
