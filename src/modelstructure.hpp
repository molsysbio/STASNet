///////////////////////////////////////////////////////////////////////////////
//
// Prototype of the ModelStructure class
// Copyright (C) 2013- Mathurin Dorel, Bertram Klinger, Nils Bluthgen
//
// Institute of Pathology and Institute for Theoretical Biology
// Charite - Universitätsmedizin Berlin - Chariteplatz 1, 10117 Berlin, Germany
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef MODELSTRUCTURE_HPP
#define MODELSTRUCTURE_HPP

#include "helper_types.hpp"
#include <vector>
#include <string>

class ModelStructure {
public:
  ModelStructure();

  ModelStructure(std::vector<std::string> names1,std::vector<std::string> names2, std::string title="");

  void set(std::vector<std::string> names1,std::vector<std::string> names2);

  void buildLinks();

  std::vector<std::string> getNames() const;
  std::string getTitle() const;
  
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
    std::string title;                // title of the structure
    int_matrix adjacencymatrix;       // adjacency matrix with 1 and 0
    int_matrix r_idx;            // give an id to each link
    std::vector<std::pair<size_t, size_t> > s_pos; // Position of each link in the adjacency matrix
    GiNaC::matrix r;       // local response matrix r
    std::vector<GiNaC::symbol> symbols_; // Symbols of the links
} ;

#include <RcppCommon.h>

RCPP_EXPOSED_AS(ModelStructure);
  

#endif 
