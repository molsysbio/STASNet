///////////////////////////////////////////////////////////////////////////////
//
// Prototype and helper classes for the fitting function
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

#ifndef FITMODEL_IN_CPP_HPP
#define FITMODEL_IN_CPP_HPP

#include <ginac/ginac.h>                 // symbolic tool box for C++ 
#include <ginac/matrix.h>                // to obtain matrices which can contain ginac expressions
#include "helper_types.hpp"
#include "model.hpp"
#include "modelset.hpp"


void fitmodel( std::vector<double> &bestfit,
           double  * bestresid,
           double_matrix &prediction,
           const Model * model,        
           const Data * data,
           std::vector<size_t> keep_constant = std::vector<size_t>());

void simulated_annealing(std::vector<double> &bestfit, double &bestresid, double_matrix &prediction, const Model *model, const Data *data, std::vector<size_t> keep_constant=std::vector<size_t>(), const double starting_temperature=1.);

// Structure for the profile likelihood
struct pl_analysis {
    // Uncertainty
    double decision;
    // Values of the thresholds
    double pointwise_threshold;
    double simultaneous_threshold;

    // Tells whether the high/low threshold has been reached for positive/negative steps
    bool ln_threshold;
    bool hn_threshold;
    bool lp_threshold;
    bool hp_threshold;

    // Indexes for which the threshold is ascendingly reached
    std::vector<size_t> negative_uncertainty;
    std::vector<size_t> positive_uncertainty;

    pl_analysis() {
        ln_threshold = false;
        hn_threshold = false;
        lp_threshold = false;
        hp_threshold = false;
    }
};
typedef struct pl_analysis pl_analysis;

void profile_likelihood( const Data &data,
            std::vector<double> parameters,
            const std::vector<size_t> keep_constant,
            std::vector< std::vector<double> > &residual_track,
            std::vector<double> &explored,
            const double param_value,
            const Model *model,
            pl_analysis &thresholds,
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


