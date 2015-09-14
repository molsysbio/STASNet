#include "modelstructure.hpp"
#include "wrap_matrices.hpp"

#include <Rcpp.h>

ModelStructure::ModelStructure() : names(), adjacencymatrix() {} 

ModelStructure::ModelStructure(std::vector<std::string> names1,std::vector<std::string> names2) : names(), adjacencymatrix() {
  set(names1,names2);
} 

std::vector<std::string> ModelStructure::getNames() { return names ; }
  
int_matrix &ModelStructure::getAdjacencyMatrix(){ return adjacencymatrix ; }

void ModelStructure::setAdjacencyMatrix(const int_matrix &adj) { 
  if ( (adj.shape()[0] == adjacencymatrix.shape()[0]) || 
       (adj.shape()[1] == adjacencymatrix.shape()[1]) )
    adjacencymatrix = adj ; 
}


void ModelStructure::set(std::vector<std::string> names1,std::vector<std::string> names2) {
    if (names1.size()!=names2.size()) 
      throw( std::invalid_argument("outgoing and incoming node names need to have the same size!") );
    names=names1;
    names.insert(names.end(), names2.begin(), names2.end());
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
      .constructor<std::vector<std::string>,std::vector<std::string> >()
      .method( "set", &ModelStructure::set )
      .property( "adjacencyMatrix", &ModelStructure::getAdjacencyMatrix )
      .method("setAdjacencyMatrix", &ModelStructure::setAdjacencyMatrix )
      .property( "names", &ModelStructure::getNames, "Names" )
      ;
}

