#ifndef MODELSET_HPP
#define MODELSET_HPP 

#include <vector>
#include "model.hpp"
#include "helper_types.hpp"

class ModelSet: public Model {
public:
    ModelSet(const GiNaC::matrix &response, const std::vector<GiNaC::symbol> &symbols, const ExperimentalDesign &expdes, const int nb_submodels, const std::vector<size_t> inter_variable_parameters=std::vector<size_t>(), bool linear_approximation=false );
    ModelSet();

    virtual void eval(const double *p, double *datax, const Data *data) const;

    void setNbModels(const int nb_submodels);
    unsigned int getNbModels() const;
    void getSubmodelsParameters(std::vector<double> &parameters);

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
