#include "modelset.hpp"

extern bool debug;
extern int verbosity;


ModelSet::ModelSet(const GiNaC::matrix &response, const std::vector<GiNaC::symbol> &symbols, const ExperimentalDesign &expdes, const int nb_submodels, const std::vector<size_t> inter_variable_parameters, bool linear_approximation ) : Model(response, symbols, expdes, linear_approximation), nb_submodels_(nb_submodels) {
    //do_init();

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

//void ModelSet::predict(const std::vector<double> &p, double_matrix &datax, const DataSet *dataset ) const {
void ModelSet::predict(const std::vector<double> &p, double_matrix &datax, const Data *data ) const {

    //const DataSet* dataset = (const DataSet*)data;
    const DataSet* dataset = dynamic_cast<const DataSet*>(data);

    // Resize datax to the size of the combined dataset
    size_t rows=dataset->unstim_data.shape()[0], cols=dataset->unstim_data.shape()[1];
    datax.resize(boost::extents[rows][cols]);
    
    cols = dataset->datas_[0].unstim_data.shape()[1];
    rows = dataset->datas_[0].unstim_data.shape()[0];
    //std::cout << rows << " rows, " << cols << " columns" << std::endl;
    std::vector<double> ptmp(nr_of_parameters_per_submodel());
    double_matrix sub_datax;
    sub_datax.resize(boost::extents[rows][cols]);
    typedef boost::multi_array_types::index_range range;
    for (size_t mod=0 ; mod < nb_submodels_ ; mod++) {
        size_t shift = mod * nr_of_parameters_per_submodel();
        for (size_t ii=0; ii< nr_of_parameters_per_submodel(); ii++ ) {
            parameters_[independent_parameters_[ii]]->set_parameter(p[shift + ii]);
            ptmp[ii] = parameters_[independent_parameters_[ii]]->eval();
        }
        Model::predict(ptmp, sub_datax, &(dataset->datas_[mod]));
        boost::detail::multi_array::multi_array_view<double, 2> datax_view = datax[ boost::indices[range(mod*rows, (mod+1)*rows, 1)][range(0, cols, 1)] ];
        copy_to_view(sub_datax, datax_view);
    }
}

void ModelSet::eval(const double *p, double *datax, const Data *data) const {
    // Returns the simulation matrix for the ModelSet, simulating each submodel with either the same set or different parameters

    const DataSet* dataset = dynamic_cast<const DataSet*>(data);
    assert(dataset->datas_.size() == nb_submodels_);
    size_t rows = dataset->datas_[0].unstim_data.shape()[0], cols = dataset->datas_[0].unstim_data.shape()[1];

    double dataxm[rows * cols];
    double ptmp[independent_parameters_.size()];
    std::copy(p, p+independent_parameters_.size(), ptmp);
    for (size_t mod=0 ; mod < nb_submodels_ ; mod++) {
        // Change the parameters that vary accross models
        for (std::vector<size_t>::const_iterator id = subparameters_ids_.begin(); id != subparameters_ids_.end(); ++id) {
            ptmp[*id] = p[mod * independent_parameters_.size() + (*id)];
            if (debug) { std::cout << "Setting parameter " << *id << " to " << ptmp[*id] << " for model " << mod << std::endl; }
        }
        Model::eval(ptmp, dataxm, &(dataset->datas_[mod]));
        std::copy(dataxm, dataxm + rows * cols, datax + mod * rows * cols);
    }
}

unsigned int ModelSet::getNbModels() const {
    return(nb_submodels_);
}

void ModelSet::setNbModels(const int nb_submodels) {
    nb_submodels_ = nb_submodels;
}
void ModelSet::setVariableParameters(const std::vector<size_t> variable_parameters) {
    for (size_t ii=0; ii<variable_parameters.size(); ii++) {
        assert(variable_parameters[ii] < nr_of_parameters_per_submodel());
    }
    subparameters_ids_ = variable_parameters;
    sort(subparameters_ids_.begin(), subparameters_ids_.end());
}

size_t ModelSet::nr_of_parameters() const {
    return( independent_parameters_.size() * nb_submodels_ );
    //return( independent_parameters_.size() + subparameters_ids_.size() * nb_submodels_ );
}

size_t ModelSet::nr_of_parameters_per_submodel() const {
    return( independent_parameters_.size() );
}

// Replaces the parameters vector to reflect the parameters that were effectively fitted
void ModelSet::getSubmodelsParameters(std::vector<double> &parameters) {
    for (size_t mod=0 ; mod < nb_submodels_ ; mod++) {
        // Use the values at the beginning of the vector for parameters that are fitted simultaneously for all models, do not change those that were fitted independently
        std::vector<size_t>::const_iterator sub_id = subparameters_ids_.begin();
        for (size_t ii=0; ii < independent_parameters_.size(); ii++) {
            if (subparameters_ids_.size() <= 0)
                parameters[mod * independent_parameters_.size() + ii] = parameters[ii];
            else if (ii != *sub_id) {
                parameters[mod * independent_parameters_.size() + ii] = parameters[ii];
            } else {
                sub_id++;
            }
        }
    }
}

void ModelSet::setNegativeInhibitions(double *p) const {
  std::vector<size_t> inhibs_ids = getInhibitorsIds();

  for (std::vector<size_t>::iterator it=inhibs_ids.begin(); it!=inhibs_ids.end(); it++) {
    for (size_t jj=0; jj<nb_submodels_; jj++) {
        if (debug) { std::cout << "submodel " << jj << ", inhibitor " << *it << std::endl; }
        p[ jj*nr_of_parameters_per_submodel() + *it ] = -std::abs(p[ jj*nr_of_parameters_per_submodel() + *it ]);
    }
  }
}
