///////////////////////////////////////////////////////////////////////////////
//
// Functions to perform the identifiability analysis of the model and generate
// new parameters from the symbolic equations
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

#include "identifiability.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/rational.hpp>
#include <iostream>  
#include <sstream>
#include <cassert>
#include <utility>
#include <fstream>

#include "mathtree.hpp"
#include "rref.hpp"

extern bool debug;
extern int verbosity;

boost::shared_ptr<MathTree::parameter> parameterlist::getParameterForExpression(GiNaC::ex e) {
  // If expression e has already been assigned to a parameter, return that.
  for (iterator iter=begin(); iter!=end(); iter++) {
    if (iter->second==e) 
      return iter->first;
  }

  // otherwise generate a new parameter, which counter is one more than the previous one
  MathTree::parameter::Ptr par(new MathTree::parameter());
  par->set_parameter(boost::shared_ptr<double>(new double(1.0)));
  push_back(std::make_pair(par,e));
  if (verbosity > 7) { std::cerr << e << " added" << std::endl;}
  return par;
}

/*
// Helper function to output a parameterlist does not work 

std::ostream &operator<<(std::ostream &os, parameterlist &a){
  for (int i=0;i<a.size();i++){
      os << a[i] << "\t";
    }
    os << std::endl;
  return os;
}
*/

// Converts an GiNaC expression into a MathTree object. 
// If reduce_products is true, products of symbols will be assigned to one parameter.
// Parameters are stored in parameterlist param
MathTree::math_item::Ptr put_into_mathtree_format (GiNaC::ex e, parameterlist &param, bool reduce_products) {
    
    MathTree::math_item::Ptr item;

    if (GiNaC::is_a<GiNaC::mul>(e)) {
        
        item = MathTree::mul::Ptr(new MathTree::mul());
        GiNaC::ex multiplier=1;
        
        for (size_t i=0; i<e.nops(); ++i) 
        {
            // This is the place where parameter reduction is happening:
            // All multiplies which are symbolic will be reduced to one...
            if (reduce_products) {
                if (!GiNaC::is_a<GiNaC::symbol>(e.op(i))) {
                    boost::dynamic_pointer_cast<MathTree::container>(item)->add_item(put_into_mathtree_format(e.op(i), param, reduce_products));
                } else {
                    multiplier*=e.op(i);
                }
            } else {
                boost::dynamic_pointer_cast<MathTree::container>(item)->add_item(put_into_mathtree_format(e.op(i), param, reduce_products));
            }
        }
        // Reduce the multiplication of symbols to one parameter 
        if (multiplier!=1) {
            MathTree::parameter::Ptr item2=param.getParameterForExpression(multiplier);
            boost::dynamic_pointer_cast<MathTree::container>(item)->add_item(item2);
        }

    } else if (GiNaC::is_a<GiNaC::add>(e)) {
        item = MathTree::add::Ptr(new MathTree::add());
        for (size_t i=0; i<e.nops(); ++i) 
            boost::dynamic_pointer_cast<MathTree::container>(item)->add_item(put_into_mathtree_format(e.op(i), param, reduce_products));

    } else if (GiNaC::is_a<GiNaC::power>(e)) {

        item = MathTree::pow::Ptr(new MathTree::pow());
        for (size_t i=0; i<e.nops(); ++i) 
            boost::dynamic_pointer_cast<MathTree::container>(item)->add_item(put_into_mathtree_format(e.op(i), param, reduce_products));

    } else if (GiNaC::is_a<GiNaC::symbol>(e)) {
        //item = param.getParameterForExpression(e);
        // Necessary, to be able to explicit single parameters
        item = MathTree::mul::Ptr(new MathTree::mul());
        boost::dynamic_pointer_cast<MathTree::mul>(item)->add_item(param.getParameterForExpression(e));

    } else if (GiNaC::is_a<GiNaC::numeric>(e)) {
        item = MathTree::numeric::Ptr(new MathTree::numeric());
        boost::dynamic_pointer_cast<MathTree::numeric>(item)->set_numeric(GiNaC::ex_to<GiNaC::numeric>(e).to_double());

    } else if (GiNaC::is_a<GiNaC::function>(e)) {
        item = MathTree::pow::Ptr(new MathTree::exp());
        for (size_t i=0; i<e.nops(); ++i) 
            boost::dynamic_pointer_cast<MathTree::container>(item)->add_item(put_into_mathtree_format(e.op(i), param, reduce_products));

    } else {
        std::cerr << "UNKNOWN TYPE: " << e;
        exit(-1);
    }

    return item;

}

// Resorts the parameter_dependency_matrix_unreduced and the parameterlist on rows
void sort ( int_matrix &old_pdmu,
        parameterlist &param) {

  //get the number of variables
  size_t i=0,j=0,vars;  
  while(old_pdmu[i][j]!=-1){j++;}
  vars=j;
  
  // std::cerr << "number of variables: " << vars << std::endl;

  //get the sum of variables for each parameter
  int nr_of_vars[param.size()];  
  int indices[param.size()];
  for (size_t i=0; i<param.size();i++){
    nr_of_vars[i]=0;
    indices[i]=i;
    for(size_t j=0; j<vars;j++){
      nr_of_vars[i]+=std::abs(old_pdmu[i][j]);
    }
  }
  
  std::sort(indices,indices+param.size(),index_cmp<int*>(nr_of_vars) );
  
  
  //change the matrix (right part) and the parameter vector accordingly
  parameterlist tmpparam = param;
  int_matrix new_pdmu = old_pdmu; 
  for (size_t i=0;i<param.size();i++){
    tmpparam[i]=param[indices[i]];
    for (size_t j=0; j<vars;j++){
      new_pdmu[i][j]=old_pdmu[indices[i]][j];
    }
  }
  old_pdmu=new_pdmu;
  param=tmpparam;
  
}


// Performs identifiability analysis on the GiNaC matrix input_matrix
void identifiability_analysis(   equation_matrix &output_matrix, 
                 std::vector<GiNaC::ex> &paths,
                 std::vector<MathTree::parameter::Ptr> &parameters,
                 double_matrix &parameter_dependency_matrix,
                 int_matrix &parameter_dependency_matrix_unreduced,
                 const symbolic_matrix &input_matrix, 
                 std::vector<GiNaC::symbol> & vars,
                 ModelStructure &structure) {


  // Generate new parameterisation
  output_matrix.resize(boost::extents[input_matrix.rows()][input_matrix.cols()]);
  parameterlist param; // List of correspondance (mathtree, GiNaC expression)
  if (debug) { std::cerr << "Converting from GiNaC to mathtree format..." << std::endl; }
  for (size_t i=0; i<input_matrix.rows(); i++) {
    for (size_t j=0; j<input_matrix.cols(); j++) {
      if(verbosity > 9) {std::cerr << i << "," << j << " : " << input_matrix(i, j) << "\t" << std::endl;}
      // Populate the param vector with (Mathtree, symbol product) pairs
      output_matrix[i][j]= put_into_mathtree_format(input_matrix(i,j).expand(),param);
    }
  }
  if(verbosity > 9) {std::cerr << std::endl;}
  // Write parameter dependencies into matrix (this is the matrix which will 
  // be put into Row Echelon form).
  // Each row represents one new parameter, the first colums are the old parameters, then the 
  // new ones.
  // It has form A|B where B is -unity matrix and A has the exponents of the products.

  // Fills the matrix A|B with the coefficients of the paths
  parameter_dependency_matrix_unreduced.resize(boost::extents[param.size()][param.size()+vars.size()]);
  if (debug) { std::cerr << "Filling the PDM..." << std::endl; }
  size_t x=0, y=0;
  for (parameterlist::iterator iter=param.begin(); iter!=param.end(); ++iter) {
        y=0;
        if(verbosity > 9) {std::cerr << "Path " << x << " : " << iter->second << std::endl;} // DEBUGGING
        assert(GiNaC::is_a<GiNaC::mul>(iter->second) || GiNaC::is_a<GiNaC::symbol>(iter->second) );
        if (GiNaC::is_a<GiNaC::mul>(iter->second)) {
          GiNaC::mul m=GiNaC::ex_to<GiNaC::mul>(iter->second); // Isn't is already a mul ?
          for (std::vector<GiNaC::symbol>::const_iterator iter2=vars.begin(); iter2!=vars.end(); ++iter2) {
            parameter_dependency_matrix_unreduced[x][y++]=(int)round(m.degree(*iter2));
          }
        } else if (GiNaC::is_a<GiNaC::symbol>(iter->second)) {
            for (std::vector<GiNaC::symbol>::const_iterator iter2=vars.begin(); iter2!=vars.end(); ++iter2) {
                if (*iter2 == iter->second)
                  parameter_dependency_matrix_unreduced[x][y++] = 1;
                else 
                  parameter_dependency_matrix_unreduced[x][y++] = 0;
            }
        }
        parameter_dependency_matrix_unreduced[x][vars.size()+x]=-1;
        x++;
  }

    // DEBUGGING prints the matrix before reduction
    if (debug) { std::cerr << "Before sorting or reduction..." << std::endl; }
    if (verbosity > 9) {
        for (size_t i=0 ; i < vars.size() ; i++) {
            std::cerr << vars[i].get_name() << "\t";
        }
        std::cerr << "\n";
        for (size_t i=0 ; i < parameter_dependency_matrix_unreduced.shape()[0] ; i++) {
            for (size_t j=0 ; j < parameter_dependency_matrix_unreduced.shape()[1] ; j++) {
                std::cerr << parameter_dependency_matrix_unreduced[i][j] << "\t";
            }
            std::cerr << "\n";
        }
        std::cerr << std::endl;
        for (int i=0 ; i < param.size() ; i++) {
            std::cerr << "Path " << i << " : " << param[i].second << std::endl;
        }
    }

  // sort the parameter dependency matrix to facilitate reduction
  sort (parameter_dependency_matrix_unreduced, param);

  parameter_dependency_matrix.resize(boost::extents[param.size()][param.size()+vars.size()]);
  //  bring into reduced row echelon form
  rational_matrix rational_matrix_for_rref;
  convert_int_to_rational_matrix(parameter_dependency_matrix_unreduced,
                 rational_matrix_for_rref);

  to_reduced_row_echelon_form(rational_matrix_for_rref);

  if (verbosity > 9) {
    for (size_t ii=0; ii < vars.size(); ii++) { std::cout << vars[ii] << "\t"; }
    std::cout << std::endl;
    printMatrix(rational_matrix_for_rref);
  }
  remove_minus_one(rational_matrix_for_rref, structure, vars, vars.size());
  if (verbosity > 5) { printMatrix(rational_matrix_for_rref); }
  to_reduced_row_echelon_form(rational_matrix_for_rref); // To be sure
  if (verbosity > 9) {
    for (size_t ii=0; ii < vars.size(); ii++) { std::cout << vars[ii] << "\t"; }
    std::cout << std::endl;
    printMatrix(rational_matrix_for_rref);
  }

  convert_rational_to_double_matrix(rational_matrix_for_rref, 
                    parameter_dependency_matrix);
  
    if (debug) { std::cerr << "After reduction" << std::endl; }
    if (verbosity > 9) {
        for (size_t i=0 ; i < vars.size() ; i++) {
            std::cerr << vars[i].get_name() << "\t";
        }
        std::cerr << "\n";
        for (size_t i=0 ; i < parameter_dependency_matrix.shape()[0] ; i++) {
            for (size_t j=0 ; j < parameter_dependency_matrix.shape()[1] ; j++) {
                std::cerr << parameter_dependency_matrix[i][j] << "\t";
            }
            std::cerr << "\n";
        }
        std::cerr << std::endl;
    }
    

  // Store parameters with mathtree format
  parameters.resize(param.size());
  paths.resize(param.size());
  for (size_t i=0; i<param.size(); i++){ 
    parameters[i]=param[i].first;
    paths[i]=param[i].second;
  }


}
 
// Remove as many -1 as possible from a row echelon matrix by swapping columns
template<typename MatrixType>
 void remove_minus_one(MatrixType &A, ModelStructure& structure, std::vector<GiNaC::symbol> & vars, size_t size) {
  matrix_traits<MatrixType> mt;
  typedef typename matrix_traits<MatrixType>::index_type index_type;

  std::vector<index_type> swaped;
  bool go_on = true;
  while (go_on) {
    std::vector<index_type> swapable;
    std::vector<index_type> leading_row;
    std::vector<int> num_minus;
    for (index_type column = mt.min_column(A); column <= size; ++column)
    {
      int one = 0;
      int minus_one = 0;
      index_type lrow = 0;
      bool invalid = false;
      for (index_type row = mt.min_row(A); row <= mt.max_row(A); ++row)
      {
        if (mt.element(A, row, column) == 1) {
          one++;
          lrow = row;
        } else if (mt.element(A, row, column) == -1) {
          minus_one++;
        } else if (mt.element(A, row, column) == 0) {
        } else {
          invalid = true;
        }
      }
      if (one == 1 && minus_one > 0 && !invalid && std::find(swaped.begin(), swaped.end(), column) == swaped.end()) {
        swapable.push_back(column);
        leading_row.push_back(lrow);
        num_minus.push_back(minus_one);
      }
    }

    if (!swapable.empty()) {
      // Select the column with the most -1
      std::vector<index_type> indices;
      for (size_t ii=0; ii < swapable.size(); ii++) { indices.push_back(ii); }
      sort(indices.begin(), indices.end(), rev_index_cmp<std::vector<int> >(num_minus));
      index_type swap_col = swapable[indices[0]];
      index_type lrow = leading_row[indices[0]];

      // Get the leading column of this row
      index_type lcol = mt.min_column(A);
      while (mt.element(A, lrow, lcol) != 1) { lcol++; }
      swap_cols(A, swap_col, lcol);
      std::swap(vars[swap_col], vars[lcol]); // Keep vars (symbols_ in Model) to avoid having to rewrite everything
      structure.swap_symbols(swap_col, lcol);
      for (index_type row = mt.min_row(A); row <= mt.max_row(A); ++row)
      {
          if (mt.element(A, row, lcol) == -1)
          {
            add_multiple_row(A, row, lrow, 1);
          }
      }
      swaped.push_back(swap_col);
      swaped.push_back(lcol);
    } else {
      go_on = false;
    }
  }
}
