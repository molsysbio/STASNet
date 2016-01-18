#include "wrap_matrices.hpp"
#include <Rcpp.h>
#include "helper_types.hpp"

//[[Rcpp:export]]
RCPP_MODULE(ExperimentalDesignEx) {
  using namespace Rcpp ;
  class_<ExperimentalDesign>("ExperimentalDesign")
    .default_constructor()
    .field( "stim_nodes", &ExperimentalDesign::stim_nodes, "Stimulated Nodes" )
    .field( "inhib_nodes", &ExperimentalDesign::inhib_nodes, "Inhibited Nodes" )
    .field( "measured_nodes", &ExperimentalDesign::measured_nodes, "Measured Nodes" )
    .field( "basal_activity", &ExperimentalDesign::basal_activity, "Basal activity" )
    .field_readonly( "stimuli", &ExperimentalDesign::stimuli )
    .field_readonly( "inhibitor", &ExperimentalDesign::inhibitor, "Inhibitor")
    .method( "set_stimuli", &ExperimentalDesign::setStimuli )
    .method( "set_inhibitor", &ExperimentalDesign::setInhibitor )
    ;

};

//[[Rcpp:export]]
RCPP_MODULE(DataEx) {
  using namespace Rcpp ;
  class_<Data>("Data")
    .default_constructor()
    .field_readonly( "unstim_data", &Data::unstim_data, "Unstimulated Data")
    .field_readonly( "stim_data", &Data::stim_data, "Stimulated Data")
    .field_readonly( "error", &Data::error, "Measurement Error")
    .field_readonly( "scale", &Data::scale, "Scaling (typically = Unstim_Data)")
    .method( "set_unstim_data", &Data::setUnstimData, "Unstimulated Data")
    .method( "set_stim_data", &Data::setStimData, "Stimulated Data")
    .method( "set_error", &Data::setError, "Measurement Error")
    .method( "set_scale", &Data::setScale, "Scaling (typically = Unstim_Data)")
    .method( "computeDataVector", &Data::computeDataVector, "Compute the data vector used in fitmodel" )
    ;

  class_<DataSet>("DataSet")
    .derives<Data>("Data")
    .default_constructor()
    .field_readonly( "datas_list", &DataSet::datas_, "List of Data composing the DataSet" )
    //.field_readonly ("datas", &DataSet::datas_, "The individual data tables of the dataset")
    .method( "addData", &DataSet::addData, "Add a data table to the dataset")
    .method( "addDataFromMatrices", &DataSet::addDataFromMatrices, "Build a data object and add it to the dataset" )
    ;

};

namespace Rcpp {
  template <> SEXP wrap( const int_matrix &a ) {
    NumericMatrix xx(a.shape()[0], a.shape()[1]);
    for (size_t i=0; i<a.shape()[0]; i++) {
      for (size_t j=0; j<a.shape()[1]; j++) {
        xx(i,j)=a[i][j];
      }
    }
    return xx;  
  }
  
  template <> int_matrix as( SEXP xxtmp ) {
    NumericMatrix xx=xxtmp;
    int_matrix a(boost::extents[xx.nrow()][xx.ncol()]);
    for (size_t i=0; i<a.shape()[0]; i++) {
      for (size_t j=0; j<a.shape()[1]; j++) {
        a[i][j]=xx(i,j);
      }
    }
    return a;  
  }
  
  template <> SEXP wrap( const double_matrix &a ) {
    NumericMatrix xx(a.shape()[0], a.shape()[1]);
    for (size_t i=0; i<a.shape()[0]; i++) {
      for (size_t j=0; j<a.shape()[1]; j++) {
        xx(i,j)=a[i][j];
      }
    }
    return xx;  
  }
  
  template <> double_matrix as( SEXP xxtmp ) {
    NumericMatrix xx=xxtmp;
    double_matrix a(boost::extents[xx.nrow()][xx.ncol()]);
    for (size_t i=0; i<a.shape()[0]; i++) {
      for (size_t j=0; j<a.shape()[1]; j++) {
        a[i][j]=xx(i,j);
      }
    }
    return a;  
  }
};


