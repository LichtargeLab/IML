#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// This snippet code will create design matrix from vcf file format.
// Input: data frame with Gene name, EA score, and Sample genotype
// Output: a numeric matrix
// Author: Saeid Parvandeh March 2019


// [[Rcpp::export]]
NumericMatrix VCFtoDM(DataFrame x, CharacterVector genes, CharacterVector patients) {
  NumericVector   ea    = x["EA"];
  CharacterVector gene  = x["GENE"];
  NumericMatrix   dMatrix(patients.size(), genes.size());
  std::fill( dMatrix.begin(), dMatrix.end(), NumericVector::get_na() ) ;
  for (int g = 0; g < genes.size(); g++){
    Rcpp::NumericVector idx;
    for(int i = 0; i < gene.size(); i++){
      if (gene[i] == genes[g]) {
        idx.push_back(i);
      }
    }
    NumericVector   sub_ea    = ea[idx];
    CharacterVector sub_gene  = gene[idx];
    for (int p = 0; p < patients.size(); p++){
      CharacterVector id  = x[p+9];
      CharacterVector sub_id  = id[idx];
      Rcpp::NumericVector EAs;
      for (int j = 0; j < sub_gene.size(); j++){
        EAs = NumericVector();
        for (int k = 0; k < sub_id.size(); k++){
          if (sub_id[k] != "./."){
            EAs.push_back(sub_ea[k]);
          }
        }
      }
      if (Rf_length(EAs)>=1){
        // dMatrix(p, g) = Rcpp::mean(EAs);
        dMatrix(p, g) = Rcpp::max(EAs);
        // dMatrix(p, g) = Rcpp::sample(EAs, 1, false)[0];
      }
    }
  }
  return(dMatrix);
}
