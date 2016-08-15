///////////////////////////////////////////////////////////////////////////////
//
// Prototype of the ModelSet class
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

#ifndef MODELSET_HPP
#define MODELSET_HPP 

#include <vector>
#include "model.hpp"
#include "helper_types.hpp"

class ModelSet: public Model {
public:
    ModelSet(const GiNaC::matrix &response, const std::vector<GiNaC::symbol> &symbols, const ExperimentalDesign &expdes, const ModelStructure &structure, const int nb_submodels, const std::vector<size_t> inter_variable_parameters=std::vector<size_t>(), bool linear_approximation=false );
    ModelSet();

    void predict(const std::vector<double> &p, double_matrix &datax, const Data *data ) const;
    virtual void eval(const double *p, double *datax, const Data *data) const;

    void setNbModels(const int nb_submodels);
    void setVariableParameters(const std::vector<size_t> variable_parameters = std::vector<size_t>());
    unsigned int getNbModels() const;
    void getSubmodelsParameters(std::vector<double> &parameters);
    void setNegativeInhibitions(double *p) const;

    virtual size_t nr_of_parameters() const;
    size_t nr_of_parameters_per_submodel() const;
    size_t modelRank() const { return nb_submodels_ * rank_; }

protected:
    unsigned int nb_submodels_;
    std::vector<size_t> subparameters_ids_; // IDs of the parameters that should be fitted independently
    /*
    std::vector< std::vector<MathTree::parameter::Ptr> > subparameters_;
    std::vector<equation_matrix> submodels_eqns_;
    */
};

#endif // MODELSET_HPP
