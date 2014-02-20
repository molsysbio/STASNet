#ifndef HELPER_TYPES_HPP
#define HELPER_TYPES_HPP

#include <boost/multi_array.hpp>
#include <boost/rational.hpp>
#include <boost/lexical_cast.hpp>
#include <ginac/ginac.h>                 // symbolic tool box for C++ 
#include <ginac/matrix.h>                // to obtain matrices wich can contain ginac expressions 
#include <iostream>
#include <string>
#include <map>
#include <limits>
#include "mathtree.hpp"

// ---------------------------------
// typedefs for divers matrix type 

typedef GiNaC::matrix symbolic_matrix;
typedef boost::multi_array< double, 2 >  double_matrix;
typedef boost::multi_array< int, 2 >  int_matrix;
typedef boost::multi_array< MathTree::math_item::Ptr,2 > equation_matrix; 
typedef boost::multi_array< boost::rational<int>,2> rational_matrix;

// Convert matrix types 
// XXX Could be replaced by numerical_cast ...
void convert_int_to_rational_matrix(const int_matrix &i, rational_matrix &r);
void convert_rational_to_double_matrix(const rational_matrix &r, double_matrix &d);

// Output matrices
std::ostream &operator<<(std::ostream &os, const double_matrix &a) ;
std::ostream &operator<<(std::ostream &os, const int_matrix &a);
std::ostream &operator<<(std::ostream &os, const equation_matrix &a) ;
std::ostream &operator<<(std::ostream &os, const std::vector<int> &a);
std::ostream &operator<<(std::ostream &os, const std::vector<double> &a);
std::istream &operator>>(std::istream &is, std::map<std::string,double> &p);
 
// ---------------------------------
// Class that stores all data
 

template<typename T> 
void copy_matrix(const T &from, T &to ) {
  std::vector<std::size_t> extendlist(2);
  extendlist[0]=from.shape()[0];
  extendlist[1]=from.shape()[1];
  to.resize(extendlist);
  to=from;
}

class ExperimentalDesign {

public:
  ExperimentalDesign &operator=(const ExperimentalDesign &expdesign);
public:
  std::vector<int> stim_nodes;
  std::vector<int> inhib_nodes;
  std::vector<int> measured_nodes;
  int_matrix stimuli;
  int_matrix inhibitor;
  std::vector<int> basal_activity;
  
  void setStimuli(int_matrix new_stimuli) { copy_matrix(new_stimuli,stimuli); }
  void setInhibitor(int_matrix new_inhibitor) { copy_matrix(new_inhibitor,inhibitor); }


  bool read_from_stream(std::istream &is);
};



class Data {

public:
  Data &operator=(const Data &data);
public:
  double_matrix unstim_data;
  double_matrix stim_data;
  double_matrix error;
  double_matrix scale;

  void setUnstimData(double_matrix new_unstim) { copy_matrix(new_unstim,unstim_data); }
  void setStimData(double_matrix new_stim) { copy_matrix(new_stim,stim_data); }
  void setError(double_matrix new_error) { copy_matrix(new_error,error); }
  void setScale(double_matrix new_scale) { copy_matrix(new_scale,scale); }


  bool read_from_stream(std::istream &is);
  bool data_consistent(const ExperimentalDesign &) const;
};


std::ostream &operator<<(std::ostream &os, Data &d);

// Helper function to determine rank
// Rank will be only determined correctly for matrices in reduced row echolon form!!!
template<typename T>
size_t rank(const T &array, double eps=0.00000001) {
  size_t rank_=0;
  for ( ; rank_ < array.shape()[0]; ++rank_ )  {
    bool finished=true;
    for (size_t j=0; j < array.shape()[1]; ++j ) {
      if (std::abs(array[rank_][j])>eps)
	finished=false;
    }
    if (finished) break;
  }
  return rank_;
}

// Helper function to read in multiarray

template<typename T>
bool read_array(std::istream &is, boost::multi_array<T,2> &matrix) {
  size_t sizex, sizey;
  is >> sizex >> sizey;
  if (!is.good()) return false;
  matrix.resize(boost::extents[sizex][sizey]);
  std::string token;
  for (size_t i=0; i<sizex; i++)
    for (size_t j=0; j<sizey; j++) {
      // Here there should be some kind of casting/checking of types!
      is >> token;
      try {
	matrix[i][j]=boost::lexical_cast<T>(token);
      } catch (std::exception e) {
	return false;
      }
    }
    if (!is.good()) return false;
  return true;
}


template<typename T>
bool read_array(std::istream &is, std::vector<T> &vector) {
  size_t sizex;
  is >> sizex;
  if (!is.good()) return false;
  vector.resize(sizex);
  std::string token;
  for (size_t i=0; i<sizex; i++) {
    // Here there should be some kind of casting/checking of types!
    is >> token;
    try {
      vector[i]=boost::lexical_cast<T>(token);
    } catch (std::exception e) {
      return false;
    }
  }
  if (!is.good()) return false;
  return true;
}
  

long unsigned int getSeed();



#endif // HELPER_TYPES_HPP
