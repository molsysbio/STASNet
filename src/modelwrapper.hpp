#ifndef MODELWRAPPER_HPP
#define MODELWRAPPER_HPP
#include "wrap_matrices.hpp"
#include "helper_types.hpp"
#include "model.hpp"
#include "modelset.hpp"
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
  virtual void setModel(ExperimentalDesign exp, ModelStructure mod);

  SEXP fitmodel_wrapper(Data data, std::vector<double> parameters, std::vector<size_t> keep_constant=std::vector<size_t>());
  SEXP fitmodelWithConstants(Data data, std::vector<double> parameters, std::vector<size_t> keep_constant);
  virtual SEXP fitmodel(Data data, std::vector<double> parameters);

  SEXP annealingFit(Data data, std::vector<double> parameters, int max_it, int max_depth);

  SEXP simulate(Data *data, std::vector<double> parameters);
  SEXP getLocalResponse( std::vector<double> p );
  SEXP profileLikelihood(Data data, std::vector<double> parameters, size_t target, const unsigned int total_steps);

  std::vector<double> getParameterFromLocalResponse( const double_matrix &response, const std::vector<double> inhib);

  bool linear_approximation;
  
  SEXP getParametersLinks();
  
  void showParameterDependencyMatrix();
  void showGUnreduced();
  SEXP getUnreducedPDM();
  SEXP getPDM();
  SEXP getParametersNames();

  void printEquation(const size_t r, const size_t c);
  void printEquationMatrix();

protected:
  GiNaC::matrix response_full_model; 
  std::vector <GiNaC::symbol> symbols_full_model;  
  int_matrix adjacency_matrix;

  Model * model;
};

class ModelSetWrapper: public ModelWrapper {
public:
    ModelSetWrapper();
    ~ModelSetWrapper();
    SEXP fitmodelset(DataSet data, std::vector<double> parameters);
    void setVariableParameters(std::vector<size_t> variable_parameters);
    void setModel(ExperimentalDesign exp, ModelStructure mod);
    void setNbModels(size_t nb_submodels);
};


#endif 
