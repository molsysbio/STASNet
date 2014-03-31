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

void profile_likelihood( const Data &data,
			std::vector<double> parameters,
            const std::vector<size_t> keep_constant,
	        std::vector< std::vector<double> > &residual_track,
            std::vector<double> &explored,
            const double param_value,
			const Model *model,
			bool* n_identifiability,
            std::vector<double> &thresholds);

void profile_likelihood( const Data &data,
			std::vector<double> parameters,
            const std::vector<size_t> keep_constant,
	        std::vector< std::vector<double> > &residual_track,
            std::vector<double> &explored,
            const double param_value,
			const Model *model,
			bool* n_identifiability,
            std::vector<double> &thresholds,
			const unsigned int total_steps);

double choose_step_size(const Data &data,
                        const std::vector<double> parameters,
                        const double param_value,
                        const std::vector<size_t> keep_constant,
                        const Model *model,
                        const double threshold,
                        const double init_residual,
                        double step);

#endif // FITMODEL_IN_CPP_HPP


