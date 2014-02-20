#ifndef WRAP_MATRICES_HPP
#define WRAP_MATRICES_HPP
#include <fitmodel/helper_types.hpp>

#include <RcppCommon.h>

namespace Rcpp {
  template <> SEXP wrap( const int_matrix &a );
  template <> int_matrix as( SEXP xx );
  template <> SEXP wrap( const double_matrix &a );
  template <> double_matrix as( SEXP xx );

};



RCPP_EXPOSED_AS(Data);
RCPP_EXPOSED_AS(ExperimentalDesign);

#endif 
