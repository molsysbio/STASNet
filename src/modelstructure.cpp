///////////////////////////////////////////////////////////////////////////////
//
// Functions of the ModelStructure class
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

#include "modelstructure.hpp"
#include "wrap_matrices.hpp"

#include <Rcpp.h>

extern bool debug;

ModelStructure::ModelStructure() : names(), adjacencymatrix() {} 

ModelStructure::ModelStructure(std::vector<std::string> names1,std::vector<std::string> names2, std::vector<std::string> namesAll, std::string title) : names(), adjacencymatrix(), title(title) {
  set(names1,names2,namesAll);
} 

std::vector<std::string> ModelStructure::getNames() const { return names ; }
std::string ModelStructure::getTitle() const { return title ; }
  
const int_matrix &ModelStructure::getAdjacencyMatrix() const { return adjacencymatrix ; }
const int_matrix &ModelStructure::getRidx() const { return r_idx ; }
const GiNaC::matrix &ModelStructure::getR() const { return r ; }
std::vector<GiNaC::symbol> ModelStructure::getSymbols() const { return symbols_ ; }
const std::vector<std::pair<size_t, size_t> >& ModelStructure::getSpos() const { return s_pos ; }

void ModelStructure::setAdjacencyMatrix(const int_matrix &adj) { 
  if ( (adj.shape()[0] == adjacencymatrix.shape()[0]) || 
       (adj.shape()[1] == adjacencymatrix.shape()[1]) )
    adjacencymatrix = adj ; 
    buildLinks();
}

void ModelStructure::buildLinks() {
 // Build symbolic "r" matrix from adjacency matrix
 symbols_.clear();
 if (debug) { std::cerr << "Building the response matrix..." << std::endl; }
  size_t size=adjacencymatrix.shape()[0];   // because we need it all the time
  r_idx.resize(boost::extents[size][size]);
  r = GiNaC::matrix(size, size);
  s_pos.clear();

  size_t c=0;  
  for (size_t i=0; i<size; i++){
    for (size_t j=0; j<size; j++){
      std::string tmpstr;                
      if (i==j) r(j,i)=-1;                         // -1 on the diagonal of the local response is required for MRA
      if (adjacencymatrix[j][i]!=0){
        if (names.size()==size) {
          tmpstr= "r_" + names[j] + "_" + names[i];      // concatenates to "r_j_i" in tmpstr which can then be converted to str
        } else {
          tmpstr = "r_" + boost::lexical_cast<std::string>(j) + "_" + boost::lexical_cast<std::string>(i); 
        }
        symbols_.push_back(GiNaC::symbol(tmpstr));   // creates another symbol r and adds one more entry to x
        s_pos.push_back(std::make_pair(j, i));
        r(j,i)=symbols_[c];                          // puts the symbol in the matrix position
        r_idx[j][i]=c++;                      // to remember the position of the symbols
      }
      else r_idx[j][i]=-1;
    }
  }
}

// Swap two symbols in the vector, and apply the changes to r_idx and s_pos to match this change
void ModelStructure::swap_symbols(size_t ii, size_t jj) {
    assert(ii >= 0 && ii < symbols_.size() && jj >= 0 && jj < symbols.size());
    std::swap(symbols_[ii], symbols_[jj]);
    std::swap(s_pos[ii], s_pos[jj]);
    for (size_t rr=0; rr<r_idx.shape()[0]; rr++) {
        for (size_t cc=0; cc<r_idx.shape()[1]; cc++) {
            if (r_idx[rr][cc] == ii) {
                r_idx[rr][cc] = jj;
            } else if (r_idx[rr][cc] == jj) {
                r_idx[rr][cc] = ii;
            }
        }
    }
}

void ModelStructure::set(std::vector<std::string> names1, std::vector<std::string> names2, std::vector<std::string> namesAll) {
    if (names1.size()!=names2.size()) 
      throw( std::invalid_argument("outgoing and incoming node names need to have the same size!") );
    
    names=names1;
    names.insert(names.end(), names2.begin(), names2.end());
    names.insert(names.end(),namesAll.begin(),namesAll.end()); 	
    std::sort (names.begin(), names.end());
    std::vector<std::string>::iterator it = std::unique(names.begin(),names.end());
    names.resize( std::distance(names.begin(),it) );
    adjacencymatrix.resize(boost::extents[names.size()][names.size()]);
    for (size_t i=0; i<names.size(); i++) {
      for (size_t j=0; j<names.size(); j++) {
        adjacencymatrix[i][j]=0;
      }
    }
      
    for (size_t i=0; i<names1.size(); i++) {
      std::vector<std::string>::iterator iter1 = std::find(names.begin(), names.end(), names1[i]);
      std::vector<std::string>::iterator iter2 = std::find(names.begin(), names.end(), names2[i]);
      size_t index1 = std::distance(names.begin(), iter1);
      size_t index2 = std::distance(names.begin(), iter2);
      adjacencymatrix[index2][index1]=1;
    }
    buildLinks();
} 

std::ostream& operator<<(std::ostream &os, const ModelStructure &mod) {
    for (size_t i=0 ; i < mod.adjacencymatrix.shape()[0] ; i++) {
        os << mod.names[i] << "\t";
    }
    os << std::endl;
    for (size_t i=0 ; i < mod.adjacencymatrix.shape()[0] ; i++) {
        for (size_t j=0 ; j < mod.adjacencymatrix.shape()[1] ; j++) {
            os << mod.adjacencymatrix[i][j] << "\t";
        }
        os << std::endl;
    }
    return os;
}
  
//
// IMPORTANT: RcppCommon.h, dann wrap, dann Rcpp.h
//

RCPP_MODULE(ModelStructureEx){
    using namespace Rcpp ;

    class_<ModelStructure>( "ModelStructure" )
      .default_constructor()
      .constructor<std::vector<std::string>,std::vector<std::string>, std::vector<std::string>, std::string >()
      .method( "set", &ModelStructure::set )
      .property( "adjacencyMatrix", &ModelStructure::getAdjacencyMatrix )
      .method("setAdjacencyMatrix", &ModelStructure::setAdjacencyMatrix )
      .property( "names", &ModelStructure::getNames, "Names" )
      .property( "title", &ModelStructure::getTitle, "Title" )
      ;
}

