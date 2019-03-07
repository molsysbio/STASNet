///////////////////////////////////////////////////////////////////////////////
//
// Definition of the functions of the Model class
// Copyright (C) 2016 Mathurin Dorel, Bertram Klinger, Nils Bluthgen
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

#include "model.hpp"
#include "fstream"

extern bool debug;
extern int verbosity;

// Constructor
// Stores all data neccessary for evaluation 
// and performs identifiability analysis

Model::Model() : rank_(0) {
};

Model::Model(const GiNaC::matrix &response, 
             const std::vector<GiNaC::symbol> &symbols, 
             const ExperimentalDesign &expdesign,
             const ModelStructure &structure,
             bool linear_approximation) : 
    response_(response), structure_(structure), exp_design_(expdesign), symbols_(symbols), rank_(0), linear_approximation_(linear_approximation)
{
    do_init();
} 

// Copy an existing model
Model::Model(const Model &model) {
    copy_Model(model);
}

Model & Model::operator=(const Model &model) {
    copy_Model(model);
    return *this;
}

void Model::copy_Model(const Model &model) {
    copy_matrix(model.model_eqns_,model_eqns_);
    parameters_=model.parameters_;
    independent_parameters_=model.independent_parameters_;
    paths_=model.paths_;
    response_ = model.response_;

    exp_design_ = model.exp_design_;
    symbols_=model.symbols_;
    copy_matrix(model.parameter_dependency_matrix_,parameter_dependency_matrix_);
    copy_matrix(model.parameter_dependency_matrix_unreduced_,parameter_dependency_matrix_unreduced_);
    rank_=model.rank_;
    linear_approximation_ = model.linear_approximation_;

    param_constraints_ = model.param_constraints_;
    constraints_ = model.constraints_;
}

void Model::do_init () {
 
    // Perform identifiability analysis
    // parameter_dependency_matrix has the following structure:
    // ( A | B ) with A is parameters.size x symbols.size and B is parameters.size * parameters_size
    // B( rank(A) .. end : 1..end ) informs us about which parameters have to be identified.

    identifiability_analysis( model_eqns_,
                    paths_,
                    parameters_,
                    parameter_dependency_matrix_,
                    parameter_dependency_matrix_unreduced_,
                    response_,
                    symbols_ ,
                    structure_);
    // What is a small number? 
    double eps=0.00001;

    // determine rank of the left submatrix for the symbols
    boost::multi_array_types::index_range 
        xrange(0,symbols_.size()), 
        yrange(0,parameter_dependency_matrix_.shape()[0]);

    rank_=rank(parameter_dependency_matrix_[ boost::indices[yrange][xrange] ], eps );
    // determine replacement rules 
    std::vector < std::pair < MathTree::math_item::Ptr, MathTree::math_item::Ptr > > replace_vector;
    
    // Use the right hand side of the matrix, rows below rank, to replace redundant parameters 
    size_t j=0;
    for ( size_t i=rank_; i< parameters_.size() ; ++i ) {
        // The first path will be replaced by the combination of the others
        while (j<parameters_.size()) {
            if (std::abs(parameter_dependency_matrix_[i][j+symbols_.size()])>eps) 
                break;
            j++;
        }
        if (j<parameters_.size()) {
            std::pair < MathTree::math_item::Ptr, MathTree::math_item::Ptr > replace_rule;
            replace_rule.first=parameters_[j]; // The power is always 1 for the first parameter on each row because the matrix is in row echelon form
            MathTree::mul::Ptr tmp2(new MathTree::mul );
            // Find the combination
            for (size_t k=j+1; k<parameters_.size() ; k++ ) {
                if (std::abs(parameter_dependency_matrix_[i][k+symbols_.size()])>eps) {
                    MathTree::pow::Ptr tmp3(new MathTree::pow );
                    tmp3->add_item(parameters_[k]);
                    tmp3->add_item(new MathTree::numeric(-parameter_dependency_matrix_[i][k+symbols_.size()]));
                    tmp2->add_item(tmp3);
                }
            }
            replace_rule.second=tmp2;
        
            replace_vector.push_back(replace_rule);
        }
    }
        
    // Collect the independent parameters
    for (size_t i=0; i<parameters_.size(); ++i) {
        bool replaced=false;
        for (size_t j=0; j<replace_vector.size(); ++j) {
            if (parameters_[i]==replace_vector[j].first) {
                replaced=true;
            }
        }
        if (!replaced) {
            independent_parameters_.push_back(i);
        }
    }

    // Removes singletons from the other parameters
    size_t k=0;
    if (verbosity > 7) {
        std::cout << "Substitution before simplification ";
        for (k=0; k<replace_vector.size(); k++) {
            replace_vector[k].first->print(std::cout);
            std::cout << " to ";
            replace_vector[k].second->print(std::cout);
            std::cout << std::endl;
        }
    }
    simplify_independent_parameters_using_k(replace_vector);
    if (verbosity > 7) {
        std::cout << "Substitution after simplification ";
        for (size_t l=k; l<replace_vector.size(); l++) {
            replace_vector[l].first->print(std::cout);
            std::cout << " to ";
            replace_vector[l].second->print(std::cout);
            std::cout << std::endl;
        }
    }
    if (verbosity > 7) {
        showParameterDependencyMatrix();
        showGUnreduced();
    }
    
    // Replacement of non identifiable links into combination of links (paths) that are identifiable in the mathtree object
    for (unsigned int i=0; i<model_eqns_.shape()[0];i++) { 
        for (unsigned int j=0; j<model_eqns_.shape()[1];j++) {
            if (verbosity > 10) {
                model_eqns_[i][j]->print(std::cout);
                std::cout << " converted to ";
            }
            if (boost::dynamic_pointer_cast<MathTree::container>(model_eqns_[i][j]).get()!=0) {
                for (size_t k=0; k<replace_vector.size(); k++) {
                    boost::dynamic_pointer_cast<MathTree::container>(model_eqns_[i][j])
                        ->replace_subitem(replace_vector[k].first,replace_vector[k].second);
                }
            }
            if (verbosity > 10) {
                model_eqns_[i][j]->print(std::cout);
                std::cout << std::endl;
            }
        }
    }

    getConstraints(param_constraints_, constraints_);
    
}

// Use the k to build reduced parameters
void Model::simplify_independent_parameters_using_k(std::vector< std::pair<MathTree::math_item::Ptr, MathTree::math_item::Ptr> > &replace_vector) {

    double eps = 0.0000001;
    size_t previous_size = parameters_.size();

    // Build the new independent parameters ("singletons") from the reduced matrix
    parameterlist singletons;
    for (size_t i=0 ; i < rank_ ; i++) {
        // Extract the link composition
        GiNaC::ex independent = 1;
        for (size_t j=0 ; j < symbols_.size() ; j++) {
            if (parameter_dependency_matrix_[i][j] == 1) {
                independent *= symbols_[j];
            } else if (parameter_dependency_matrix_[i][j] == -1) {
                independent /= symbols_[j];
            }
        }
        // Create the corresponding parameter
        MathTree::parameter::Ptr par(new MathTree::parameter);
        par->set_parameter(boost::shared_ptr<double>(new double(1.0)));

        // Link the parameter and the expression and add them to the corresponding vectors
        singletons.push_back(std::make_pair(par, independent));
        parameters_.push_back(par);
        paths_.push_back(independent);
        if (verbosity > 8) {
            std::cout << "New singleton : " << independent << std::endl;
        }

        // Add the new parameters to the matrices (can be duplicate of old parameters, not important)
        // Necessary for the conversions between links and independent paths
        //
        // Reduced matrix
        parameter_dependency_matrix_.resize(boost::extents[parameter_dependency_matrix_.shape()[0]+1][parameter_dependency_matrix_.shape()[1]+1]);
        // Fill the last column with zero
        for (size_t row=0 ; row < parameter_dependency_matrix_.shape()[0] ; row++) {
            parameter_dependency_matrix_[row][parameter_dependency_matrix_.shape()[1]-1] = 0;
        }
        // Copy the link composition of the path from the reduced matrix
        for (size_t col=0 ; col < symbols_.size() ; col++) {
            parameter_dependency_matrix_[parameter_dependency_matrix_.shape()[0]-1][col] = parameter_dependency_matrix_[i][col];
        }
        // An independent parameter only depends on itself
        for (size_t l=symbols_.size() ; l < (parameter_dependency_matrix_.shape()[1]-1) ; l++) {
            parameter_dependency_matrix_[parameter_dependency_matrix_.shape()[0]-1][l] = 0;
        }
        parameter_dependency_matrix_[parameter_dependency_matrix_.shape()[0]-1][parameter_dependency_matrix_.shape()[1]-1] = -1;
        //
        // Unreduced matrix
        parameter_dependency_matrix_unreduced_.resize(boost::extents[parameter_dependency_matrix_unreduced_.shape()[0]+1][parameter_dependency_matrix_unreduced_.shape()[1]+1]);
        // Fill the last column with zero
        for (size_t row=0 ; row < parameter_dependency_matrix_unreduced_.shape()[0] ; row++) {
            parameter_dependency_matrix_unreduced_[row][parameter_dependency_matrix_unreduced_.shape()[1]-1] = 0;
        }
        // Copy the link composition of the path from the reduced matrix
        for (size_t col=0 ; col < symbols_.size() ; col++) {
            parameter_dependency_matrix_unreduced_[parameter_dependency_matrix_unreduced_.shape()[0]-1][col] = parameter_dependency_matrix_[i][col];
        }
        // An independent parameter only depends on itself
        for (size_t l=symbols_.size() ; l < (parameter_dependency_matrix_unreduced_.shape()[1]-1) ; l++) {
            parameter_dependency_matrix_unreduced_[parameter_dependency_matrix_unreduced_.shape()[0]-1][l] = 0;
        }
        parameter_dependency_matrix_unreduced_[parameter_dependency_matrix_unreduced_.shape()[0]-1][parameter_dependency_matrix_unreduced_.shape()[1]-1] = -1;
    }

    size_t ipi;
    // Put the expression of k in term of paths in a matrix
    // Express the old parameters as combination of the new ones
    int_matrix p_from_k_unreduced;
    p_from_k_unreduced.resize(boost::extents[independent_parameters_.size()][2*independent_parameters_.size()]);
    for (size_t i=0 ; i < independent_parameters_.size() ; i++) {
        for (size_t j=0 ; j < independent_parameters_.size() ; j++) {
            p_from_k_unreduced[i][j] = parameter_dependency_matrix_[i][symbols_.size() + independent_parameters_[j]];
            p_from_k_unreduced[i][independent_parameters_.size() + j] = 0;
        }
        // Identity, will contain expression of paths in terms of k after reduction
        p_from_k_unreduced[i][independent_parameters_.size() + i] = -1;
    }

    rational_matrix rational_for_rref;
    convert_int_to_rational_matrix(p_from_k_unreduced, rational_for_rref);

    // Invert the matrix to get the expression of paths in k
    to_reduced_row_echelon_form(rational_for_rref);

    double_matrix p_from_k;
    convert_rational_to_double_matrix(rational_for_rref, p_from_k);

    // Replace the paths by the corresponding combination of k
    for (size_t i=0 ; i < independent_parameters_.size() ; i++) {
        ipi = independent_parameters_[i];

        // Find which singletons (k) are in the old parameter
        MathTree::mul::Ptr tmp (new MathTree::mul);
        for (size_t j=0 ; j < independent_parameters_.size() ; j++) {
            if (p_from_k[i][independent_parameters_.size() + j] == 1) {
                tmp->add_item(parameters_[previous_size + j]);
            } else if (p_from_k[i][independent_parameters_.size() + j] == -1) {
                MathTree::pow::Ptr tmp2(new MathTree::pow);
                tmp2->add_item(parameters_[previous_size + j]);
                tmp2->add_item(new MathTree::numeric(-1.0));
                tmp->add_item(tmp2);
            }
        }
        // Replace the pointers in the equation matrix
        replace_vector.push_back(std::make_pair(parameters_[ipi], tmp));

        // Replace the old parameters by the reduced ones in the reduced dependency matrix
        for (size_t row=0 ; row < previous_size ; row++) {
            if (std::abs(parameter_dependency_matrix_[row][symbols_.size() + ipi]) > eps) {
                for (size_t l=0 ; l < independent_parameters_.size() ; l++) {
                    parameter_dependency_matrix_[row][symbols_.size() + previous_size + l] += parameter_dependency_matrix_[row][symbols_.size() + ipi] * p_from_k[i][independent_parameters_.size() + l];
                }
                parameter_dependency_matrix_[row][symbols_.size() + ipi] = 0;
            }
        }

    }
    
    // Tell where the new independent parameters are in the matrices
    for (size_t i=0 ; i < independent_parameters_.size() ; i++) {
        independent_parameters_[i] = previous_size + i;
    }

}

// Compares pairs, first by the first term, then by the second if the first are equals
template<class T> struct pair_rev_index_cmp {
    pair_rev_index_cmp(const T arr) : arr(arr) {}
    bool operator()(const size_t a, const size_t b) const
    {
        if (arr.first[a] != arr.first[b]) {
            return arr.first[a] < arr.first[b]; // Leftmost 1 get prioririty for the row echelon form
        } else {
            return arr.second[a] > arr.second[b]; // Complex parameters are on top, important for small ones
        }
    }
    const T arr;
};

// Locate the sums to the power of -1 (denominators) from the GiNaC expression for the constraints analysis
void recurse( GiNaC::ex e, std::vector<GiNaC::ex> &set) {
    MathTree::math_item::Ptr tmp;
    if (GiNaC::is_a<GiNaC::power>(e)) {
        if (e.nops()==2) { // Sanity check
            if (e.op(1)==-1) {
                if (GiNaC::is_a<GiNaC::add>(e.op(0))) {
                    // Add each sum once to the set
                    std::vector<GiNaC::ex>::iterator iter=set.begin();
                    while ( !((iter==set.end()) || ((*iter) == e.op(0)))) {
                        iter++;
                    }
                    if (iter==set.end()) {
                        set.push_back(e.op(0));
                    }
                }

            }
        } else {
            std::cout << "Versteh ich nicht!" << std::endl;
            exit(-1);
        }
    } else {
        for (size_t i=0; i<e.nops(); ++i) {
            recurse(e.op(i),set);

        }
    }
}

double Model::getPeneltyForConstraints(const double *p) const {

    std::vector<double> p_id(Model::nr_of_parameters());
    for (size_t i=0; i< Model::nr_of_parameters(); i++ ) { p_id[i]=p[i]; }
    std::vector<double> p_new;
    convert_identifiables_to_original_parameter(p_new, p_id) ;

    for (parameterlist::const_iterator iter=param_constraints_.begin(); 
             iter!=param_constraints_.end(); iter++) {
        GiNaC::ex search=iter->second;
        size_t i=0;
        while ((i<symbols_.size() ) && (symbols_[i]!=search)) { i++; }
        if (i==symbols_.size()) { 
            std::cerr << "getPeneltyForConstraints: HIER IST WAS FAUL!!!" << std::endl; exit(-1); }
        //      std::cout <<    symbols_[i] << " " << iter->second << " " << p_new[i] << " | ";
        iter->first->set_parameter(p_new[i]);
    }
    double penelty=0.0, tmp;
    for (std::vector<MathTree::math_item::Ptr>::const_iterator iter=constraints_.begin(); 
             iter!=constraints_.end(); iter++ ) {
        tmp=(*iter)->eval() ;
        if (tmp>=0) penelty+=tmp;
    }
    return penelty;
}

// Collect the sums to the power of -1 from the GiNaC equation matrix and put them under the mathtree format in a vector (put into constraints_ by do_init)
void Model::getConstraints( parameterlist &params, std::vector<MathTree::math_item::Ptr> &equations) {
    std::vector<GiNaC::ex> set;
    for (size_t i=0; i<response_.rows(); i++) 
        for (size_t j=0; j<response_.cols(); j++) {
            GiNaC::ex e=response_(i,j).expand();
            recurse(e,set);
        }
    params.clear();
    equations.clear();
    for (std::vector<GiNaC::ex>::iterator iter=set.begin(); iter!=set.end(); iter++) {
        equations.push_back(put_into_mathtree_format (*iter, params, false));

    }
}

// Give the value of the model for each condition for the set of parameters p
void Model::predict(const std::vector<double> &p, double_matrix &datax, const Data *data, const bool with_offset ) const {
    size_t rows=data->unstim_data.shape()[0], cols=data->unstim_data.shape()[1];

    datax.resize(boost::extents[rows][cols]);

    //  assert( (unsigned int)m == nr_of_parameters() );
    //  assert( (unsigned int)n == rows*cols );
    
    for (size_t i=0; i< Model::nr_of_parameters(); i++ ) {
        parameters_[independent_parameters_[i]]->set_parameter(p[i]);
    }
    for (unsigned int i=0; i<cols;i++) { 
        for (unsigned int j=0; j<rows;j++) {
            if (linear_approximation_) {
                datax[j][i]=( data->unstim_data[j][i] + model_eqns_[i*rows+j][0]->eval()*data->scale[j][i]);
            } else {
                datax[j][i]= data->unstim_data[j][i] * exp( model_eqns_[i*rows+j][0]->eval());
                if (with_offset && datax[j][i] < data->scale[j][i]) {
                    datax[j][i] = data->scale[j][i]; // Minimum value is blank
                }
            }
        }
    }
}

// Give the value of the model for each condition for the set of parameters p normalised by the error for the fit
void Model::eval(const double *p,double *datax, const Data *data ) const {

    size_t rows=data->unstim_data.shape()[0], cols=data->unstim_data.shape()[1];
    //  assert( (unsigned int)m == nr_of_parameters() );
    //  assert( (unsigned int)n == rows*cols );
    //std::cout << "Unstim_data: " << rows << ", " << cols << std::endl;
    //std::cout << "Stim_data: " << data->stim_data.shape()[0] << ", " << data->stim_data.shape()[1] << std::endl;
    
    for (size_t i=0; i< Model::nr_of_parameters(); i++ ) {
        parameters_[independent_parameters_[i]]->set_parameter(p[i]);
    }
        
    double penelty=getPeneltyForConstraints(p);
    for (unsigned int i=0; i<cols;i++) { 
        for (unsigned int j=0; j<rows;j++) {
            //if (penelty>1) {
            //    // Positive feedback loops create forking, so we eliminate the parameters sets which involve such feedback
            //    datax[i*rows+j] = 100000*penelty*data->stim_data[j][i]/data->error[j][i];
            //} else
            if (linear_approximation_) {
                datax[i*rows+j]=( data->unstim_data[j][i] + model_eqns_[i*rows+j][0]->eval()*data->scale[j][i])/(sqrt(2) * data->error[j][i]);
            } else {
                datax[i*rows+j]= ( data->unstim_data[j][i] * exp( model_eqns_[i*rows+j][0]->eval()))/(sqrt(2) * data->error[j][i]);
                if (datax[i*rows+j] < data->scale[j][i]/(sqrt(2) * data->error[j][i])) {
                    datax[i*rows+j] = data->scale[j][i]/(sqrt(2) * data->error[j][i]);
                }
            }
            if (std::isnan(data->error[j][i]) || std::isnan(data->stim_data[j][i])) {
                datax[i*rows+j]=0;
            } else if ((std::isnan(datax[i*rows+j])) || (std::isinf(datax[i*rows+j])) ) {
                if (verbosity > 6) {
                    std::cerr << datax[i*rows+j] << ", " << data->unstim_data[j][i] << ", " << data->error[j][i] << ", " << model_eqns_[i*rows+j][0]->eval() << std::endl;
                }
                datax[i*rows+j]=5*data->stim_data[j][i]/(sqrt(2) * data->error[j][i]);
            } else if ((datax[i*rows+j]<0.00001) || (datax[i*rows+j]>100000)){
                // to exclude extreme values, where the algorithm can't find a way out somehow 
                datax[i*rows+j]=log(datax[i*rows+j])*data->stim_data[j][i]/(sqrt(2) * data->error[j][i]);
            }
        }
    }

}

// Calculates the residual for the set of parameters p
double Model::score(const double *p, const Data *data) const {
    int number_of_measurements=data->stim_data.shape()[1] * data->stim_data.shape()[0]; 

    // Data and simulation normalised by the error
    double simData[number_of_measurements];
    eval(p, simData, data);
    double datax[number_of_measurements];
    for (size_t i=0; i<data->stim_data.shape()[1]; i++ ) {
        for (size_t j=0; j<data->stim_data.shape()[0]; j++) {
            if (std::isnan(data->error[j][i]) || std::isnan(data->stim_data[j][i])) {
                datax[i*data->stim_data.shape()[0]+j]=0;
            } else {
                datax[i*data->stim_data.shape()[0]+j]=data->stim_data[j][i]/data->error[j][i];
            }
        }
    }
    assert(data->stim_data.shape()[1]*data->stim_data.shape()[0] == (size_t)number_of_measurements);

    // Sum of squares
    double residual = 0;
    for (size_t i=0 ; i < number_of_measurements ; i++) {
        residual += pow(simData[i] - datax[i], 2);
    }
    return residual;
}

size_t Model::nr_of_parameters() const { return independent_parameters_.size(); }

int Model::find_parameter(std::string name) {
    double eps = 0.000000001;
    int symbol_pos=-1;
    for (size_t j=0; j<symbols_.size(); ++j) {
        if (symbols_[j].get_name() == name) { symbol_pos=j; }
    }
    if (symbol_pos<0) return -3;    // No such symbol
    int row=-1;
    for (size_t i=0; i<rank_; ++i) {
        if (std::abs(parameter_dependency_matrix_[i][symbol_pos])>eps) {
            if (row>-1) { return -2; // Parameter is in more than one parameter combi
            } else {
    row=i;
            }
        }
    }
    for (size_t j=0; j<symbols_.size(); ++j) {
        if (std::abs(parameter_dependency_matrix_[row][j])>eps) {
            if ((int)j!=symbol_pos) return -1; // Parameter is in a combi, not alone
        }
    }

    int ind_param_nr = -1;
    for (size_t j=0; j<independent_parameters_.size(); ++j) {
        if (std::abs(parameter_dependency_matrix_[row][independent_parameters_[j]+symbols_.size()])>eps) {
            if (ind_param_nr>-1) return -1;
            ind_param_nr=j;
        }
    }
    
    return ind_param_nr;
}

// TODO prints a human readable report about the identifiable parameter combinations
void Model::print_parameter_report(std::ostream &os, const std::vector<double> &d) {

    assert(d.size()==independent_parameters_.size());
    double eps = 0.000000001;

    if (rank_>0) {
        os << "RANK: " << rank_ << std::endl;
        for (size_t i=0; i<rank_; ++i) {
            bool first=true;
            for (size_t j=0; j<symbols_.size(); ++j) {
                if (std::abs(parameter_dependency_matrix_[i][j])>eps) {
                    //  expressions_to_be_substituted.push_back(symbols_[j]);
                    if (first)
                        first=false;
                    else
                        os << "*";
                    os << symbols_[j] << "\t";
                    if (std::abs(parameter_dependency_matrix_[i][j]-1)>eps) 
                        os << "^" << -std::abs(parameter_dependency_matrix_[i][j]);
                }
            }
            double tmp=1;
            os << "(subparams:" << " ";
            for (size_t j=0; j<independent_parameters_.size(); ++j) 
            if (std::abs(parameter_dependency_matrix_[i][independent_parameters_[j]+symbols_.size()])>eps) {
                tmp*=std::pow(d[j],-parameter_dependency_matrix_[i][independent_parameters_[j]+symbols_.size()]);
                os << j << " " << 
                     -parameter_dependency_matrix_[i][independent_parameters_[j]+symbols_.size()] << " " << d[j] << ",";
            }
            os << ") = " << tmp << std::endl;
        }
    } else {
            for (size_t j=0; j<independent_parameters_.size(); ++j) 
                os << symbols_[independent_parameters_[j]] << "\t=\t" << d[independent_parameters_[j]] << " (" << independent_parameters_[j] << ")" << std::endl;
    }
}

// Returns the strings corresponding to the expression of the parameters
void Model::getParametersLinks(std::vector<std::string> &description) const {
    description = std::vector<std::string>();
    if (debug && verbosity > 4) { std::cerr << "Retrieving parameters links" << std::endl; }
    for (size_t j=0; j<independent_parameters_.size(); ++j) {
        description.push_back(to_string(paths_[independent_parameters_[j]]));
        if (verbosity > 6) { std::cout << paths_[independent_parameters_[j]] << std::endl; }
    }
}

// Displays the parameter dependency matrix G reduced
void Model::showParameterDependencyMatrix() {
    std::string output="";

    for (size_t i=0 ; i < symbols_.size() ; i++) {
        output += symbols_[i].get_name() + "\t";
    }
    output += "\n";
    for (size_t i=0 ; i < parameter_dependency_matrix_.shape()[0] ; i++) {
        for (size_t j=0 ; j < parameter_dependency_matrix_.shape()[1] ; j++) {
            output += to_string(parameter_dependency_matrix_[i][j]) + "\t";
            if (j < symbols_.size()) {
                output += "\t";
            }
        }
        output += "\n";
    }
    std::cout << output << std::endl;
}

double_matrix Model::getParameterDependencyMatrix() {
    return(parameter_dependency_matrix_);
}

std::vector< std::vector<int> > Model::getUnreducedParameterDependencyMatrix() {
    std::vector< std::vector<int> > return_matrix;
    for (size_t i=0 ; i < parameter_dependency_matrix_unreduced_.shape()[0] ; i++) {
        return_matrix.push_back(std::vector<int>());
        for (size_t j=0 ; j < parameter_dependency_matrix_unreduced_.shape()[1] ; j++) {
            return_matrix[i].push_back(parameter_dependency_matrix_unreduced_[i][j]);
        }
    }
    return return_matrix;
}

// Displays G before the reduction pb because it is modified directly by the reduction
void Model::showGUnreduced() {
    std::string output="";

    for (size_t i=0 ; i < symbols_.size() ; i++) {
        output += symbols_[i].get_name() + "\t";
    }
    output += "\n";
    for (size_t i=0 ; i < parameter_dependency_matrix_unreduced_.shape()[0] ; i++) {
        for (size_t j=0 ; j < parameter_dependency_matrix_unreduced_.shape()[1] ; j++) {
            output += to_string(parameter_dependency_matrix_unreduced_[i][j]) + "\t";
            if (j < symbols_.size()) {
                output += "\t";
            }
        }
        output += "\n";
    }
    std::cout << output << std::endl;
}

// TODO prints a human readable report about the identifiable parameter combinations
void Model::print_dot_file(std::ostream &os, const std::vector<double> &d, const int_matrix &origadj, const int_matrix &adj, const std::vector<std::string> &names)
{
    os << "digraph G {" << std::endl;

    for (size_t i=0; i<origadj.shape()[1]; i++) {
        for (size_t j=0; j<origadj.shape()[2]; j++) {
            if ((origadj[i][j]==0)&&(adj[i][j]==1)) {
    os << names[j] << " -> " << names[i] << " [color=\"blue\"]" << std::endl;
            }
            if (origadj[i][j]!=0) {
    os << names[j] << " -> " << names[i];
    if (adj[i][j]==0) {
        os << " [color=\"grey\"]" ;
    }
    os << std::endl;
            }
            
        }
    }
    os << "}" << std::endl;


    double_matrix determined_parametermatrix;
    /*
    double eps = 0.000000001;

    if (rank_>0) {
        for (size_t i=0; i<rank_; ++i) {
            bool first=true;
            for (size_t j=0; j<symbols_.size(); ++j) {
    if (std::abs(parameter_dependency_matrix_[i][j])>eps) {
        //  expressions_to_be_substituted.push_back(symbols_[j]);
        if (first)
            first=false;
        else
            os << "*";
        os << symbols_[j] << "\t";
        if (std::abs(parameter_dependency_matrix_[i][j]-1)>eps) 
            os << "^" << std::abs(parameter_dependency_matrix_[i][j]);
    }
            }
            double tmp=1;
            os << "(subparams:" << " ";
            for (size_t j=0; j<independent_parameters_.size(); ++j) 
    if (std::abs(parameter_dependency_matrix_[i][independent_parameters_[j]+symbols_.size()])>eps) {
        tmp*=std::pow(d[j],-parameter_dependency_matrix_[i][independent_parameters_[j]+symbols_.size()]);
         os << i << " " << j << " " << 
             -parameter_dependency_matrix_[i][independent_parameters_[j]+symbols_.size()] << " " << d[j] << ",";
    }
            os << ") = " << tmp << std::endl;
        }
 
    } else {
        for (size_t j=0; j<independent_parameters_.size(); ++j) 
            os << symbols_[independent_parameters_[j]] << "\t=\t" << d[independent_parameters_[j]] << " (" << independent_parameters_[j] << ")" << std::endl;

            }*/
}

// Print the equation for condition r and measurements c
void Model::printEquation(const size_t r, const size_t c) {
    assert(r >= 0 && r < exp_design_.measured_nodes.size());
    assert(c >= 0 && c < exp_design_.stimuli.shape()[0]);
    std::cout << response_[r*exp_design_.stimuli.shape()[0] + c] << std::endl;
}

void Model::printEquationMatrix() {
    for (size_t rr=0; rr < exp_design_.measured_nodes.size(); rr++) {
        for (size_t cc=0; cc < exp_design_.stimuli.shape()[0]; cc++) {
            std::cout << response_[rr*exp_design_.stimuli.shape()[0] + cc] << " ;"<< std::endl;
        }
        std::cout << std::endl;
    }
}

// Ensure inhibitors values are negative as they are treated as such in the equations.
void Model::setNegativeInhibitions(double *p) const {
  std::vector<size_t> inhibs_ids = getInhibitorsIds();
  for (std::vector<size_t>::iterator it=inhibs_ids.begin(); it!=inhibs_ids.end(); it++) {
    p[*it] = -std::abs(p[*it]);
  }
}

std::vector<size_t> Model::getInhibitorsIds() const {
  std::vector<size_t> inhibs_ids;
  std::vector<std::string> plinks;

  getParametersLinks(plinks);
  for (size_t ii=0; ii < plinks.size() ; ii++) {
    if (plinks[ii][0] == 'i') {
        inhibs_ids.push_back(ii);
    }
  }
  return(inhibs_ids);
}


// TODO Maps identifiable parameters to a set of possible original parameters
//void Model::return_one_set_of_orignal_parameters( std::vector<double> &p1, const std::vector<double> &p2) 
//{
//}

void Model::convert_parameter_map_into_identifiable(std::vector<double> &p1 , 
                                const parametermap &pm) {
    std::vector<double> p2(symbols_.size(),0);
    for (size_t i=0; i<symbols_.size(); i++) {
        if (pm.find(symbols_[i].get_name())==pm.end()) {
            p2[i]=0;
            std::cerr << "No value for " << symbols_[i].get_name() << ". assuming 0" << std::endl;
        } else {
            p2[i]=pm.find(symbols_[i].get_name())->second;
        }
    }
    convert_original_parameter_into_identifiable(p1,p2);
}

void Model::convert_original_parameter_into_identifiable( std::vector<double> &p1, const std::vector<double> &p2) 
{
 
    p1.resize(independent_parameters_.size());
    assert(p2.size()==symbols_.size());

    for (size_t i=0; i<independent_parameters_.size(); ++i) {
        assert(std::abs(parameter_dependency_matrix_unreduced_[independent_parameters_[i]][symbols_.size()+independent_parameters_[i]]+1.0)<0.000001);
        double tmp=1.0;
        for (size_t j=0; j<symbols_.size(); ++j) {
            if (std::abs(parameter_dependency_matrix_unreduced_[independent_parameters_[i]][j])>0.000001) {
                tmp*=pow(p2[j],parameter_dependency_matrix_unreduced_[independent_parameters_[i]][j]);
            }
        }
        p1[i]=tmp;
    }
}

// Converts paths to links
void Model::convert_identifiables_to_original_parameter(std::vector<double> &p_new, const std::vector<double> &p_old) const
{
    p_new.resize(symbols_.size());
    assert(p_old.size()==independent_parameters_.size() && p_new.size()==symbols_.size());
    // check for occurances of original parameters, which are initialized to 1
    for (size_t i=0; i<rank_; i++) {
        std::vector<size_t> not_zero;
        for (size_t k=0;k<symbols_.size();k++){
            if (std::abs(parameter_dependency_matrix_[i][k])>0.000001) {not_zero.push_back(k);}
        }
        //std::cout << symbols_[not_zero[0]]; //
        if (not_zero.size() >1){
            for(size_t l=1;l<not_zero.size();l++){
                p_new[not_zero[l]]=1.0;
                //std::cout << "*" << symbols_[not_zero[l]]; //
            }
        }
        //std::cout << " = "; //
        // Rebuild the value of the original parameters using the values of the independent parameters and the path informations
        double tmp=1.0;
        for (size_t j=0; j<independent_parameters_.size(); j++) {
            if(std::abs(parameter_dependency_matrix_[i][symbols_.size()+independent_parameters_[j]])>0.000001){
                tmp*=pow(p_old[j],-parameter_dependency_matrix_[i][symbols_.size()+independent_parameters_[j]]);
                //std::cout << paths_[independent_parameters_[j]] << " * "; //
            }
        }
        //std::cout << std::endl; //
        p_new[not_zero[0]]=tmp; // We give all the effect to the first link of the path
    }
}

void Model::print_original_parameters(std::ostream &os, std::vector<double> &p) {
    assert( p.size()==symbols_.size());
    for (size_t i=0; i<symbols_.size(); i++) {
        os << symbols_[i] << "\t" << p[i] << std::endl;
    }
}

void Model::printSymbols() {
    for (size_t i=0; i<symbols_.size(); i++) {
        std::cout << symbols_[i] << "\t";
    }
    std::cout << std::endl;
}

void Model::convert_original_parameter_to_response_matrix( 
    double_matrix &d, 
    std::vector<double> &inh, 
    const std::vector<double> &p) 
{
    size_t counter=0;
    const int_matrix adj = structure_.getAdjacencyMatrix();
    size_t size=adj.shape()[0];
    double_matrix::extent_gen extents;
    d.resize(extents[adj.shape()[0]][adj.shape()[1]]);
    for (size_t i=0; i<size; i++){
        for (size_t j=0; j<size; j++){
            d[j][i]=0;
        }
    }
    const std::vector<std::pair<size_t, size_t> > s_pos = structure_.getSpos();
    for (size_t ii=0; ii < s_pos.size(); ii++) {
        assert(counter<p.size());
        d[s_pos[ii].first][s_pos[ii].second] = p[counter++];
    }
    inh.clear();
    for (; counter < p.size(); counter++) {
        inh.push_back(p[counter]);
    }
}


void Model::convert_response_matrix_to_original_parameter(std::vector<double> &p, const double_matrix &d, const std::vector<double> &inh) {
    p.clear();
    const std::vector<std::pair<size_t, size_t> > s_pos = structure_.getSpos();
    for (size_t ii=0; ii < s_pos.size(); ii++) {
        p.push_back(d[s_pos[ii].first][s_pos[ii].second]);
    }
    for (size_t i=0; i<inh.size(); ++i) {
        p.push_back(inh[i]);
    }
}

