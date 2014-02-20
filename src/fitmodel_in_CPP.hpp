#ifndef FITMODEL_IN_CPP_HPP
#define FITMODEL_IN_CPP_HPP

#include <ginac/ginac.h>                 // symbolic tool box for C++ 
#include <ginac/matrix.h>                // to obtain matrices which can contain ginac expressions
#include "helper_types.hpp"
#include "model.hpp"

void fitmodel( std::vector<double> &bestfit,
	       double  * bestresid,
	       double_matrix &prediction,
	       const Model * model,	       
	       const Data * data,
	       std::vector<size_t> keep_constant = std::vector<size_t>());


#endif // FITMODEL_IN_CPP_HPP


