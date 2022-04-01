#include <string>
#include <vector>
#include <numeric>
#include <map>
#include <Rcpp.h>
using namespace Rcpp;

// Split a string at a delimiter and return a set
std::set<std::string> splitSet(const std::string& s, char delimiter) {
    std::set<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.insert(token);
    }
    return tokens;
}

// Merge entries (e.g., sequences) corresponding to the same mutant
// [[Rcpp::export]]
DataFrame mergeValues(std::vector<std::string> mutNamesIn, std::vector<std::string> valuesIn,
                      char delimiter = ',') {
    std::map<std::string, std::set<std::string>> valueSet;
    std::map<std::string, std::set<std::string>>::iterator valueSetIt;

    for (size_t i=0; i<mutNamesIn.size(); i++) {
        if ((valueSetIt = valueSet.find(mutNamesIn[i])) != valueSet.end()) {
            // mutant already present
            std::set<std::string> sst = splitSet(valuesIn[i], delimiter);
            (*valueSetIt).second.insert(sst.begin(), sst.end());
        } else {
            // mutant not yet present
            valueSet.insert(std::pair<std::string,std::set<std::string>>(mutNamesIn[i], 
                                                                         splitSet(valuesIn[i], delimiter)));
        }
    }
    
    size_t dfLen = valueSet.size();
    std::vector<std::string> dfValue(dfLen, ""), dfName(dfLen, "");
    
    int j = 0;
    for (valueSetIt = valueSet.begin(); valueSetIt != valueSet.end(); valueSetIt++) {
        std::vector<std::string> valueVector((*valueSetIt).second.begin(),
                                             (*valueSetIt).second.end());
        std::string collapsedValue = "";
        for (size_t i = 0; i < valueVector.size(); i++) {
            collapsedValue += valueVector[i] + std::string(1, delimiter);
        }
        if (!collapsedValue.empty()) {
            collapsedValue.pop_back(); // remove final ","
        }
        dfName[j] = (*valueSetIt).first;
        dfValue[j] = collapsedValue;
        j++;
    }
    
    DataFrame df = DataFrame::create(Named("mutantNameColl") = dfName,
                                     Named("valueColl") = dfValue);
    
    return df;
} 
