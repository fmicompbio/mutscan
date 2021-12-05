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
//' @param metric A character scalar defining the string distance metric. One
//'   of \code{"hamming"} (default), \code{"hamming_shift"} or
//'   \code{"levenshtein"}.
//' @param nThreads numeric(1), number of threads to use for parallel processing.
//' 
//' @return An integer vector of the same length as \code{x}.
//' 
//' @example
//' calcNearestStringDist(c("lazy", "hazy", "crazy"))
//' calcNearestStringDist(c("lazy", "hazy", "crazy"), metric = "hamming_shift")
//' calcNearestStringDist(c("lazy", "hazy", "crazy"), metric = "levenshtein")
//' 
//' @export
// [[Rcpp::export]]
IntegerVector calcNearestStringDist(std::vector<std::string> x,
                                    std::string metric = "hamming",
                                    int nThreads = 1) {
    // declare variables
    size_t i, j, n = x.size();
    int dist1;
    int (*distance)(const std::string&, const std::string&, int); // pointer to function to calcluate string distance
    IntegerVector dists(n, INT_MAX);
    
    // set distance function pointer
    if (metric == "hamming") {
        distance = &hamming_distance;
        
    } else if (metric == "hamming_shift") {
        distance = &hamming_shift_distance;
        
    } else if (metric == "levenshtein") {
        distance = &levenshtein_distance;
        
    } else {
        stop("unknown distance metric '%s'", metric);
    }
    
    // choose algorithm depending on number of threads
    if (nThreads < 3) { // serial version: O(n*log(n))
        for (i = 0; i < n - 1; i++) {
            for (j = i + 1; j < n; j++) {
                dist1 = distance(x[i], x[j], -1);
                if (dists[i] > dist1)
                    dists[i] = dist1;
                if (dists[j] > dist1)
                    dists[j] = dist1;
            }
        }
    } else {             // parallel version (O(n^2))
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nThreads) private(i,j,dist1) shared(dists)
#endif
        for (i = 0; i < n; i++) {
            int mindist = INT_MAX;
            for (j = 0; j < n; j++) {
                if (i == j)
                    continue;
                dist1 = distance(x[i], x[j], -1);
                if (mindist > dist1)
                    mindist = dist1;
            }
#ifdef _OPENMP
            #pragma omp critical
#endif
            dists[i] = mindist;
        }
    }
    
    return dists;
}
