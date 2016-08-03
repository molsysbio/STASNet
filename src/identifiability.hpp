///////////////////////////////////////////////////////////////////////////////
//
// Prototypes of the function to perform the identifiability analysis of the
// model and generate new parameters from the symbolic equations
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

#ifndef IDENTIFIABILITY_HPP
#define IDENTIFIABILITY_HPP

#include <ginac/symbol.h> 
#include <vector> 
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/multi_array.hpp>
#include "mathtree.hpp"
#include "helper_types.hpp"
#include "modelstructure.hpp"


typedef std::pair<boost::shared_ptr<MathTree::parameter>,GiNaC::ex> parametermap;

class parameterlist : public std::vector<parametermap> {
public: 
  boost::shared_ptr<MathTree::parameter> getParameterForExpression(GiNaC::ex e);
};

void identifiability_analysis(
  // Output: Entries of the matrix after having removed non-identifiable parameters flattened (row by row).
  equation_matrix &output_matrix, 
  // Output: Vector containing all new parameters
  std::vector<GiNaC::ex> &paths,
  std::vector<MathTree::parameter::Ptr> &parameters,
  // Output: Dependency of new parameters on old parameters (and vice versa...)
  double_matrix &parameter_dependency_matrix,
  // Output: Dependency of new parameters on old parameters (not in rref)
  int_matrix &parameter_dependency_matrix_unreduced,

  // Input: Matrix of expression
  const symbolic_matrix &input_matrix, 
  // Input: vector of symbols which are the parameters in the input-matrix
  std::vector<GiNaC::symbol> & vars,
  // Input: the structure of the network, which will be modified to limit the number of -1 exponents
  ModelStructure &structure);

MathTree::math_item::Ptr put_into_mathtree_format (GiNaC::ex e, parameterlist &param, bool reduce_products=true);
template<typename MatrixType> void remove_minus_one(MatrixType &A, ModelStructure& structure, std::vector<GiNaC::symbol> & vars, size_t size);

#endif // IDENTIFIABILITY_HPP
