///////////////////////////////////////////////////////////////////////////////
//
// Prototype of the Model class
// Copyright (C) 2016 Mathurin Dorel, Bertram Klinger, Nils Bluthgen
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

#ifndef MODEL_HPP
#define MODEL_HPP 


#include <ginac/ginac.h>                 // symbolic tool box for C++ 
#include <ginac/matrix.h>                // to obtain matrices wich can contain ginac expressions 
#include <vector>
#include <map>
#include <iostream>
#include <boost/multi_array.hpp>
#include "mathtree.hpp"
#include "identifiability.hpp"
#include "helper_types.hpp"
#include "rref.hpp"
#include "modelstructure.hpp"

class Model {
public:

  typedef std::map<std::string,double> parametermap;

   Model(const GiNaC::matrix &response, 
     const std::vector<GiNaC::symbol> &symbols, 
     const ExperimentalDesign &expdes,
     const ModelStructure &structure,
     bool linear_approximation,
     bool log_data_=false
     );
             
  Model();

  Model(const Model &);

  Model &operator=(const Model &m);
  void copy_Model(const Model &model);

public:

  // Evaluation of the model
  virtual void predict(const std::vector<double> &p, double_matrix &datax, const Data *data, const bool with_offset=true ) const;
  virtual void eval(const double *p, double *datax, const Data *data) const;
  double score(const double *p, const Data *data) const;
  //  evaluates the model with parameters p  and returns datax.
  // This is to be called from the optimiser
  inline void eval(const double *p,double *datax,int m, int n, const Data *data) const {
    eval(p,datax,data);
  }
  inline void useLog() { log_data_=true; }

  // returns the number of identifiable parameter combinations
  virtual size_t nr_of_parameters() const;
  virtual size_t nr_of_parameters_per_submodel() const {};
  virtual void getSubmodelsParameters(std::vector<double> &parameters) {} // Used for ModelSet

  // TODO prints a human readable report about the identifiable parameter combinations
  void print_parameter_report(std::ostream &os, const std::vector<double> &d);

  void getParametersLinks(std::vector<std::string> &description) const;
  void showParameterDependencyMatrix();
  void showGUnreduced();
  std::vector< std::vector<int> > getUnreducedParameterDependencyMatrix();
  double_matrix getParameterDependencyMatrix();

  int find_parameter(std::string name) ;
  // TODO Maps identifiable parameters to a set of possible original parameters
  void print_dot_file(std::ostream &os, const std::vector<double> &d, const int_matrix &origadj, const int_matrix &adj, const std::vector<std::string> &names);
  
  //  void return_one_set_of_orignal_parameters( std::vector<double> &p1, const std::vector<double> &p2);
  void convert_original_parameter_into_identifiable( std::vector<double> &p1, const std::vector<double> &p2 );
  void convert_parameter_map_into_identifiable(std::vector<double> &p1 , const parametermap &p2) ;
  void convert_identifiables_to_original_parameter(std::vector<double> &p_new,  const std::vector<double> &p_old) const;
  void convert_original_parameter_to_response_matrix(double_matrix &d, std::vector<double> &inh, const std::vector<double> &p);
  void convert_response_matrix_to_original_parameter(std::vector<double> &p, const double_matrix &d, const std::vector<double> &inh);
  
  void print_original_parameters(std::ostream &os, std::vector<double> &p);
  ExperimentalDesign &exp_design() { return exp_design_; }
  const ExperimentalDesign &exp_design() const { return exp_design_; }

  double getPeneltyForConstraints(const double *p) const;
  void getConstraints( parameterlist &params, std::vector<MathTree::math_item::Ptr> &equations);


  void printResponse() {
    for (std::vector<MathTree::math_item::Ptr>::iterator iter=constraints_.begin();
     iter!=constraints_.end(); iter++) {
      std::cout << **iter << std::endl;
    }
    for (parameterlist::iterator iter=param_constraints_.begin(); iter!=param_constraints_.end(); iter++) {
      std::cout << *iter->first << " " << iter->second << std::endl;
    }
      
    //    exit(0);
  }
  void printEquation(const size_t r, const size_t c);
  void printEquationMatrix();
  virtual void setNegativeInhibitions(double *p) const;
  std::vector<size_t> getInhibitorsIds() const;
  void printSymbols();

  virtual size_t modelRank() const { return rank_; }
  
  virtual void setNbModels(const int nb_submodels) {}
  virtual void setVariableParameters(std::vector<size_t> variable_parameters) {}

protected:
  // The equations of the reduced model under GiNaC or mathtree format
  equation_matrix model_eqns_;
  std::vector<MathTree::parameter::Ptr> parameters_;
  std::vector<size_t> independent_parameters_;
  std::vector<GiNaC::ex> paths_;
  GiNaC::matrix response_;
  
  // 
  ExperimentalDesign exp_design_;
  ModelStructure structure_;
  std::vector< GiNaC::symbol > symbols_;
  double_matrix parameter_dependency_matrix_;
  int_matrix parameter_dependency_matrix_unreduced_;
  size_t rank_;
  bool linear_approximation_;
  bool log_data_;

  void do_init();
  void simplify_independent_parameters_using_k(std::vector< std::pair<MathTree::math_item::Ptr, MathTree::math_item::Ptr> > &replace_vector);
void simplify_independent_parameters_using_subtraction(std::vector< std::pair<MathTree::math_item::Ptr, MathTree::math_item::Ptr> > &replace_vector);

  // ----------
  // Contains the constraints to avoid bifurcations

  parameterlist param_constraints_;
  std::vector<MathTree::math_item::Ptr> constraints_;
  
  
  // Hier muss es noch die Abhaengigkeitsmatrix zwischen den Parametern geben
  
};


#endif // MODEL_HPP
