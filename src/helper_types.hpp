///////////////////////////////////////////////////////////////////////////////
//
// Prototypes of helper classe to hold the data and the experimental design
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

std::vector< std::vector<double> > LHSampling (const int nb_samples, const int sample_size, const int decimals=10);
std::vector< std::vector<double> > normalLHS (const int nb_samples, const int sample_size, const int sd=1, const double decimals=10);
double uniform_sampling(double precision=10000);

// Functors to sort indexes of an array with std::sort
template<class T> struct index_cmp {
  index_cmp(const T arr) : arr(arr) {}
  bool operator()(const size_t a, const size_t b) const
  {
      return arr[a] > arr[b]; // Normal sort
  }
  const T arr;
};
template<class T> struct rev_index_cmp {
  rev_index_cmp(const T arr) : arr(arr) {}
  bool operator()(const size_t a, const size_t b) const
  {
      return arr[a] < arr[b]; // Reverted sort
  }
  const T arr;
};

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
 
// Print a 2-dimensions boost array
template<typename Matrix>
void print_matrix(const Matrix &m) {
    for (size_t i=0 ; i < m.shape()[0] ; i++) {
        for (size_t j=0 ; j < m.shape()[1] ; j++) {
            std::cout << m[i][j];
        }
        std::cout << std::endl;
    }
}

// Converts any data type into a string
template <typename T>
std::string to_string(const T &a) {
    std::ostringstream oss;
    oss << a;
    return oss.str();
}

template<typename T> 
void copy_matrix(const T &from, T &to ) {
  std::vector<std::size_t> extendlist(2);
  extendlist[0]=from.shape()[0];
  extendlist[1]=from.shape()[1];
  to.resize(extendlist);
  to=from;
}

template<typename T, typename V>
void copy_to_view(const T &from, V &to) {
    assert(from.shape()[0] == to.shape()[0]);
    assert(from.shape()[1] == to.shape()[1]);
    to = from;
}

template<typename T, typename S>
void rbind_matrix(const T &from, S &dest) {
    std::cout << "Another round" << std::endl;
    std::cout << "from dim: " << from.shape()[0] << ", " << from.shape()[1] << std::endl; //Debug
    std::cout << "dest dim: " << dest.shape()[0] << ", " << dest.shape()[1] << std::endl; //Debug
    std::vector<std::size_t> extendlist(2);
    extendlist[0] = dest.shape()[0] + from.shape()[0];
    extendlist[1] = dest.shape()[1];
    std::cout << "Ready for extension" << std::endl; // Debug
    if (dest.shape()[1] == 0) { extendlist[1] = from.shape()[1]; }
    assert(extendlist[1] == from.shape()[1]); // <---- generates a "malloc(): memory corruption"
    std::cout << "new size: " << extendlist[0] << ", " << extendlist[1] << std::endl; // Debug
    const size_t *old_shape = dest.shape();
    std::cout << "Performing resize" << std::endl; // Debug
    dest.resize(extendlist);
    std::cout << "yolo" << std::endl;
    for (int ii=0; ii < from.shape()[0]; ii++) {
        for (int jj=0; jj < from.shape()[1]; jj++) {
            //std::cout << "Ok for " << ii << ", " << jj << std::endl;
            dest[old_shape[0]+ii][jj] = from[ii][jj];
        }
    }
}

// ---------------------------------
// Class that stores all data
 
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
  Data();
public:
  double_matrix unstim_data;
  double_matrix stim_data;
  double_matrix error;
  double_matrix scale;
//  std::vector<double> dataVector;
  size_t nb_measurements;
  double *dataVector;
  bool dataVectorComputed;

  void setUnstimData(double_matrix new_unstim) { copy_matrix(new_unstim,unstim_data); }
  void setStimData(double_matrix new_stim) { copy_matrix(new_stim,stim_data); computeDataVector(); }
  void setError(double_matrix new_error) { copy_matrix(new_error,error); computeDataVector(); }
  void setScale(double_matrix new_scale) { copy_matrix(new_scale,scale); }

  bool read_from_stream(std::istream &is);
  virtual bool data_consistent(const ExperimentalDesign &expdesign) const;

  virtual void computeDataVector();
};

// Stores several dataset, used with ModelSet
class DataSet: public Data {
public:
    DataSet();
    ~DataSet();
public:
    std::vector<Data> datas_;
public:
    void addData(Data &data, bool doDataVectorComputation=false);
    void addDataFromMatrices(double_matrix unstim_data, double_matrix stim_data, double_matrix error, double_matrix scale, bool doDataVectorComputation=false);
    virtual bool data_consistent(const ExperimentalDesign &expdesign) const;
    virtual void computeDataVector();
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
