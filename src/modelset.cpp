#include "modelset.hpp"

extern bool debug;
extern int verbosity;


ModelSet::ModelSet(const GiNaC::matrix &response, const std::vector<GiNaC::symbol> &symbols, const ExperimentalDesign &expdes, const int nb_submodels, const std::vector<size_t> inter_variable_parameters, bool linear_approximation ) : Model(response, symbols, expdes, linear_approximation), nb_submodels_(nb_submodels) {
    do_init();

    /*
    // Copy the equation for each model, the parameters are copied as pointers and thus the same in every models
    for (size_t ii=0 ; ii < nb_submodels_ ; ii++) {
        submodels_eqns_.append(model_eqns_);
    }
    */
    // Make some parameters independent accross models
    subparameters_ids_ = inter_variable_parameters;
    /*
    for (size_t id=0 ; id < subparameters_ids_.size() ; id++) {
        subparameters_.append(std::vector<MathTree::parameter::Ptr>(nb_submodels_));
        for (size_t mod=0 ; mod < nb_submodels_ ; mod++) {
            submodels_eqns_.replace(parameters_[subparameters_ids_[id]], subparameters_[id][mod]);
        }
    }
    */
}

void ModelSet::eval(const double *p, double *datax, const Data **data) const {
    // Returns the simulation matrix for the ModelSet, simulating each submodel with either the same set or different parameters

    size_t rows=data[0]->unstim_data.shape()[0], cols=data[0]->unstim_data.shape()[1];

    double dataxm[rows * cols];
    double *ptmp;
    std::copy(p, p+parameters_.size(), ptmp);
    for (size_t mod=0 ; mod < nb_submodels_ ; mod++) {
        std::copy(p, p * independent_parameters_.size(), ptmp);
        for (std::vector<size_t>::const_iterator id = subparameters_ids_.begin(); id != subparameters_ids_.end(); ++id) {
            ptmp[*id] = p[mod * independent_parameters_.size() + (*id)];
        }
        Model::eval(ptmp, dataxm, data[mod]);
        std::copy(dataxm, dataxm + rows * cols, datax + mod * rows * cols);
    }
}

unsigned int ModelSet::getNbModels() const {
    return(nb_submodels_);
}

size_t ModelSet::nr_of_parameters() const {
    return( independent_parameters_.size() + subparameters_ids_.size() * nb_submodels );
}

