#ifndef MODELSTRUCTURE_HPP
#define MODELSTRUCTURE_HPP

#include "helper_types.hpp"
#include <vector>
#include <string>

class ModelStructure {
public:
  ModelStructure();

  ModelStructure(std::vector<std::string> names1,std::vector<std::string> names2);

  void set(std::vector<std::string> names1,std::vector<std::string> names2);

  void buildLinks();

  std::vector<std::string> getNames() const;
  
  const int_matrix &getAdjacencyMatrix() const;
  const int_matrix &getRidx() const;
  const GiNaC::matrix &getR() const;
  std::vector<GiNaC::symbol> getSymbols() const;
  const std::vector<std::pair<size_t, size_t> >& getSpos() const;

  void swap_symbols(size_t ii, size_t jj);

  void setAdjacencyMatrix(const int_matrix &adj);

  friend std::ostream& operator<<(std::ostream &os, const ModelStructure &mod);
    
 private:
    std::vector<std::string> names;   // names of the nodes
    int_matrix adjacencymatrix;       // adjacency matrix with 1 and 0
    int_matrix r_idx;            // give an id to each link
    std::vector<std::pair<size_t, size_t> > s_pos; // Position of each link in the adjacency matrix
    GiNaC::matrix r;       // local response matrix r
    std::vector<GiNaC::symbol> symbols_; // Symbols of the links
} ;

#include <RcppCommon.h>

RCPP_EXPOSED_AS(ModelStructure);
  

#endif 
