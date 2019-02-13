///////////////////////////////////////////////////////////////////////////////
//
// Prototype of the Rcpp wrapper for the Model class and the fitting function
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
  virtual void setModel(ExperimentalDesign exp, ModelStructure mod, bool log_data=false);

  SEXP fitmodel_wrapper(Data data, std::vector<double> parameters, std::vector<size_t> keep_constant=std::vector<size_t>(), std::string optimizer="levmar");
  SEXP fitmodelWithConstants(Data data, std::vector<double> parameters, std::vector<size_t> keep_constant, std::string optimizer);
  virtual SEXP fitmodel(Data data, std::vector<double> parameters, std::string optimizer);

  SEXP annealingFit(Data data, std::vector<double> parameters, std::vector<size_t> keep_constant);

  SEXP simulate(Data *data, std::vector<double> parameters);
  SEXP simulateWithOffset(Data *data, std::vector<double> parameters);
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
  void printSymbols();

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
    SEXP getSubParameterIDs();
};


#endif 
