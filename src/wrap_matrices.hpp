///////////////////////////////////////////////////////////////////////////////
//
// Prototypes of the Rcpp wrappers for the helper classes
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

#ifndef WRAP_MATRICES_HPP
#define WRAP_MATRICES_HPP
#include "helper_types.hpp"

#include <RcppCommon.h>

namespace Rcpp {
  template <> SEXP wrap( const int_matrix &a );
  template <> int_matrix as( SEXP xx );
  template <> SEXP wrap( const double_matrix &a );
  template <> double_matrix as( SEXP xx );
};

RCPP_EXPOSED_AS(DataSet); // Only produces a Rcpp::as for DataSet
RCPP_EXPOSED_CLASS(Data); // Produces both Rcpp::wrap and Rcpp::as for Data
RCPP_EXPOSED_AS(ExperimentalDesign);

#endif 
