#ifndef MODELWRAPPER_HPP
#define MODELWRAPPER_HPP
#include "wrap_matrices.hpp"
#include <fitmodel/helper_types.hpp>
#include <fitmodel/model.hpp>
#include <vector>
#include <string>
#include "modelstructure.hpp"

class ModelWrapper {
public:
  ModelWrapper();
	
  ~ModelWrapper();


  void printResponse();

  size_t modelRank();

  size_t nr_of_parameters() ;

  
  bool model_design_consistent(ExperimentalDesign &exp, ModelStructure &mod);

  void setModel(ExperimentalDesign exp, ModelStructure mod);

  SEXP fitmodel(Data data, std::vector<double> parameters);

  SEXP simulate(Data data, std::vector<double> parameters);

  SEXP getLocalResponse( std::vector<double> p );

  // Check that it is correct
  SEXP profileLikelihood(Data data, std::vector<double> parameters, size_t target, const unsigned int total_steps, const double step_size);


  std::vector<double> getParameterFromLocalResponse( const double_matrix &response, const std::vector<double> inhib);

  bool linear_approximation;

private:
  GiNaC::matrix response_full_model; 
  std::vector <GiNaC::symbol> symbols_full_model;  
  int_matrix adjacency_matrix;

  Model * model;
};


#endif 
