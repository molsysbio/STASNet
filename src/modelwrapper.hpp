#ifndef MODELWRAPPER_HPP
#define MODELWRAPPER_HPP
#include "wrap_matrices.hpp"
#include "helper_types.hpp"
#include "model.hpp"
#include "rref.hpp"
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
  SEXP annealingFit(Data data, std::vector<double> parameters, int max_it, int max_depth);

  SEXP simulate(Data data, std::vector<double> parameters);

  SEXP getLocalResponse( std::vector<double> p );

  SEXP profileLikelihood(Data data, std::vector<double> parameters, size_t target, const unsigned int total_steps);

  SEXP parallelPL(Data data, std::vector<double> parameters, const unsigned int total_steps);

  std::vector<double> getParameterFromLocalResponse( const double_matrix &response, const std::vector<double> inhib);

  bool linear_approximation;
  
  SEXP getParametersLinks();
  
  void showParameterDependencyMatrix();
  void showGUnreduced();
  SEXP getUnreducedPDM();
  SEXP getParametersNames();

  void printEquation(const size_t r, const size_t c);

private:
  GiNaC::matrix response_full_model; 
  std::vector <GiNaC::symbol> symbols_full_model;  
  int_matrix adjacency_matrix;

  Model * model;
};


#endif 
