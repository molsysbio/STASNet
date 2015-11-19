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

    double *dataxm;
    double *ptmp;
    std::copy(p, p+parameters_.size(), ptmp);
    for (size_t mod=0 ; mod < nb_submodels_ ; mod++) {
        std::copy(p, p + nr_of_parameters(), ptmp);
        for (std::vector<size_t>::const_iterator id = subparameters_ids_.begin(); id != subparameters_ids_.end(); ++id) {
            ptmp[*id] = p[mod * nr_of_parameters() +(*id)];
        }
        Model::eval(ptmp, dataxm, data[mod]);
        //evalSubmodel(ptmp, dataxm, data[mod], mod);
    }
}

unsigned int ModelSet::getNbModels() const {
    return(nb_submodels_);
}

/*
ModelSet::evalSubmodel(const double *p, double *datax, const Data *data, size_t submodel_id) {
    size_t rows=data->unstim_data.shape()[0], cols=data->unstim_data.shape()[1];
    
    for (size_t i=0; i< nr_of_parameters(); i++ ) {
        parameters_[independent_parameters_[i]]->set_parameter(p[i]);
    }
        
    double penelty=getPeneltyForConstraints(p);
    for (unsigned int i=0; i<cols;i++) { 
        for (unsigned int j=0; j<rows;j++) {
            if (penelty>1) {
                // Positive feedback loops create forking, so we eliminate the parameters sets which involve such feedback
                datax[i*rows+j] = 100000*penelty*data->stim_data[j][i]/data->error[j][i];
            } else if (linear_approximation_) {
                datax[i*rows+j]=( data->unstim_data[j][i] + submodels_eqns_[submodel_id][i*rows+j][0]->eval()*data->scale[j][i])/data->error[j][i];
            } else {
                datax[i*rows+j]=( data->unstim_data[j][i] *exp( submodels_eqns_[submodel_id][i*rows+j][0]->eval()))/data->error[j][i];
            }
//}

            if (std::isnan(data->error[j][i]) || std::isnan(data->stim_data[j][i])) {
                datax[i*rows+j]=0;
            } else if ((std::isnan(datax[i*rows+j])) || 
             (std::isinf(datax[i*rows+j])) ) {
                if (verbosity > 10) {
                    std::cerr << datax[i*rows+j] << ", " << data->unstim_data[j][i] << ", " << data->error[j][i] << ", " << submodels_eqns_[submodel_id][i*rows+j][0]->eval() << std::endl;
                }
                datax[i*rows+j]=5*data->stim_data[j][i]/data->error[j][i];
            } else if ((datax[i*rows+j]<0.00001) || (datax[i*rows+j]>100000)){
    // to exclude extreme values, where the algorithm can't find a way out somehow 
                datax[i*rows+j]=log(datax[i*rows+j])*data->stim_data[j][i]/data->error[j][i];
            } 
        }
    }

}
*/

