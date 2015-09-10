#ifndef MODELSTRUCTURE_HPP
#define MODELSTRUCTURE_HPP

#include <fitmodel/helper_types.hpp>
#include <vector>
#include <string>

class ModelStructure {
public:
  ModelStructure();

  ModelStructure(std::vector<std::string> names1,std::vector<std::string> names2);

  void set(std::vector<std::string> names1,std::vector<std::string> names2);

  std::vector<std::string> getNames();
  
  // Why does this return a reference ?
  int_matrix &getAdjacencyMatrix();

  void setAdjacencyMatrix(const int_matrix &adj);

  friend std::ostream& operator<<(std::ostream &os, const ModelStructure &mod);
    
 private:
    std::vector<std::string> names;
    int_matrix adjacencymatrix;
} ;

#include <RcppCommon.h>

RCPP_EXPOSED_AS(ModelStructure);
  

#endif 
