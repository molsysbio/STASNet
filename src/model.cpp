#include "model.hpp"
#include "fstream"


// Constructor
// Stores all data neccessary for evaluation 
// and performs identifiability analysis

Model::Model() : rank_(0) {
};

// COPY CONSTRUCTOR IS BUGGY!
Model::Model(const Model &model) : symbols_(model.symbols_), response_(model.response_), rank_(0), linear_approximation_(true) {
    std::cerr << "Scheisse 2!" << std::endl;
    exit(1);
    
    do_init();
}


Model::Model(const GiNaC::matrix &response, 
             const std::vector<GiNaC::symbol> &symbols, 
             const ExperimentalDesign &expdesign,
             bool linear_approximation) : 

    exp_design_(expdesign), symbols_(symbols), response_(response), rank_(0), linear_approximation_(linear_approximation)

{
     do_init();
} 



Model & Model::operator =(const Model &model) {
    exp_design_ = model.exp_design_;
    symbols_=model.symbols_;

    response_ = model.response_;
    copy_matrix(model.model_eqns_,model_eqns_);
    parameters_=model.parameters_;
    independent_parameters_=model.independent_parameters_;

    copy_matrix(model.parameter_dependency_matrix_,parameter_dependency_matrix_);
    copy_matrix(model.parameter_dependency_matrix_unreduced_,parameter_dependency_matrix_unreduced_);

    rank_=model.rank_;
    return *this;
}


void Model::do_init () {
 

    // Perform identifiability analysis
    // parameter_dependency_matrix has the following structure:
    // ( A | B ) with A is parameters.size x symbols.size and B is parameters.size * parameters_size
    // B( rank(A) .. end : 1..end ) informs us about which parameters have to be identified.

    identifiability_analysis( model_eqns_,
                    parameters_,
                    parameter_dependency_matrix_,
                    parameter_dependency_matrix_unreduced_,
                    response_,
                    symbols_ );
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
        while (j<parameters_.size()) { 
            if (std::abs(parameter_dependency_matrix_[i][j+symbols_.size()])>eps) 
                break;
            j++;
        }
        if (j<parameters_.size()) {
            std::pair < MathTree::math_item::Ptr, MathTree::math_item::Ptr > replace_rule;
            replace_rule.first=parameters_[j];
            MathTree::mul::Ptr tmp2(new MathTree::mul );
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
    

    
    for (unsigned int i=0; i<model_eqns_.shape()[0];i++) { 
        for (unsigned int j=0; j<model_eqns_.shape()[1];j++) {
            if (boost::dynamic_pointer_cast<MathTree::container>(model_eqns_[i][j]).get()!=0) {
                for (size_t k=0; k<replace_vector.size(); k++) {
                    boost::dynamic_pointer_cast<MathTree::container>(model_eqns_[i][j])
                        ->replace_subitem(replace_vector[k].first,replace_vector[k].second);
                }
            }
        }

    }

    getConstraints(param_constraints_, constraints_);
    
}

void recurse( GiNaC::ex e, std::vector<GiNaC::ex> &set) {
    MathTree::math_item::Ptr tmp;
    if (GiNaC::is_a<GiNaC::power>(e)) {
        if (e.nops()==2) {
            if (e.op(1)==-1) {
                if (GiNaC::is_a<GiNaC::add>(e.op(0))) {
                    for (size_t i=0; i<e.op(0).nops(); ++i) {
                        if (GiNaC::is_a<GiNaC::numeric>(e.op(0).op(i))) {
                            if (GiNaC::ex_to<GiNaC::numeric>(e.op(0).op(i)).to_double()!=-1) {
                    std::cout << "Was soll das?" << 
                        GiNaC::ex_to<GiNaC::numeric>(e.op(0).op(i)).to_double() << std::endl;
                    exit(-1);
                            }
                            
                        }
                    }

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
    std::vector<double> p_id(nr_of_parameters());
    for (size_t i=0; i< nr_of_parameters(); i++ ) { p_id[i]=p[i]; }
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

void Model::predict(const std::vector<double> &p, double_matrix &datax, const Data *data ) const {


    size_t rows=data->unstim_data.shape()[0], cols=data->unstim_data.shape()[1];

    datax.resize(boost::extents[rows][cols]);

    //  assert( (unsigned int)m == nr_of_parameters() );
    //  assert( (unsigned int)n == rows*cols );
    
    for (size_t i=0; i< nr_of_parameters(); i++ ) {
        parameters_[independent_parameters_[i]]->set_parameter(p[i]);
    }
    for (unsigned int i=0; i<cols;i++) { 
        for (unsigned int j=0; j<rows;j++) {
            if (linear_approximation_) {
                datax[j][i]=( data->unstim_data[j][i] + model_eqns_[i*rows+j][0]->eval()*data->scale[j][i]);
            } else {
                datax[j][i]=( data->unstim_data[j][i] *exp( model_eqns_[i*rows+j][0]->eval()));
            }
        }
    }
}

void Model::eval(const double *p,double *datax, const Data *data ) const {

    size_t rows=data->unstim_data.shape()[0], cols=data->unstim_data.shape()[1];
    //  assert( (unsigned int)m == nr_of_parameters() );
    //  assert( (unsigned int)n == rows*cols );
    
    for (size_t i=0; i< nr_of_parameters(); i++ ) {
        parameters_[independent_parameters_[i]]->set_parameter(p[i]);

    }
        
 
    for (unsigned int i=0; i<cols;i++) { 
        for (unsigned int j=0; j<rows;j++) {
            if (linear_approximation_) {
                datax[i*rows+j]=( data->unstim_data[j][i] + model_eqns_[i*rows+j][0]->eval()*data->scale[j][i])/data->error[j][i];
            } else {
                datax[i*rows+j]=( data->unstim_data[j][i] *exp( model_eqns_[i*rows+j][0]->eval()))/data->error[j][i];
            }

            if (std::isnan(data->error[j][i])) {
    datax[i*rows+j]=0;
            } else if ((std::isnan(datax[i*rows+j])) || 
         (std::isinf(datax[i*rows+j])) ) {
    datax[i*rows+j]=5*data->stim_data[j][i]/data->error[j][i];
            } else if ((datax[i*rows+j]<0.00001) || (datax[i*rows+j]>100000)){
    // to exclude extreme values, where the algorithm can't find a way out somehow 
    datax[i*rows+j]=log(datax[i*rows+j])*data->stim_data[j][i]/data->error[j][i];
            } 
        }
    }
    
 

    double penelty=getPeneltyForConstraints(p);
    if (penelty>1) {
        //      std::cerr << "P:" << penelty << " " ;
        for (unsigned int i=0; i<cols;i++) { 
            for (unsigned int j=0; j<rows;j++) {
    
    datax[i*rows+j]=datax[i*rows+j]+100000*penelty*data->stim_data[j][i]/data->error[j][i];
            }
        }
    }
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
void Model::print_parameter_report(std::ostream &os, const std::vector<double> &d)
{

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

// Put the litteral description of the identifiable parameters in a vector of strings
std::string to_string(double a); // Dirty hack
void Model::getParametersLinks(std::vector<std::string> &description) {
    
    double eps = 0.000000001;
    for (size_t i=0; i<rank_; ++i) {
        description.push_back(std::string(""));
        bool first=true;
        for (size_t j=0; j<symbols_.size(); ++j) {
            if (std::abs(parameter_dependency_matrix_[i][j])>eps) {
                if (first) {
                    first=false;
                }
                else {
                    description[i] = description[i] + "*";
                }
                description[i] = description[i] + symbols_[j].get_name();
                if (std::abs(parameter_dependency_matrix_[i][j]-1)>eps) {
                    description[i] = description[i] + "^";
                    description[i] = description[i] + to_string(-std::abs(parameter_dependency_matrix_[i][j]));
                }
            }
        }
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
        }
        output += "\n";
    }
    std::cout << output << std::endl;
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
        }
        output += "\n";
    }
    std::cout << output << std::endl;
}

// Global function, only used here, dirty hack part II
std::string to_string(double a) {
    std::ostringstream oss;
    oss << a;
    return oss.str();
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
            std::cerr << "No value for " << symbols_[i].get_name() << ". assuming 0" << std::cout;
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

void Model::convert_identifiables_to_original_parameter(std::vector<double> &p_new, const std::vector<double> &p_old) const
{
    p_new.resize(symbols_.size());
    assert(p_old.size()==independent_parameters_.size() && p_new.size()==symbols_.size());
    // check all rows i for occurances of original parameters k of which take the 1. which will be replaced by the product  
    for (size_t i=0; i<rank_; i++) {
    std::vector<size_t> not_zero;
    for (size_t k=0;k<symbols_.size();k++){
        if (std::abs(parameter_dependency_matrix_[i][k])>0.000001) {not_zero.push_back(k);}
    }
    if (not_zero.size() >1){
        for(size_t l=1;l<not_zero.size();l++){
            p_new[not_zero[l]]=1.0;
        }
    }
    double tmp=1.0;
    for (size_t j=0; j<independent_parameters_.size(); j++) {
        if(std::abs(parameter_dependency_matrix_[i][symbols_.size()+independent_parameters_[j]])>0.000001){
            tmp*=pow(p_old[j],-parameter_dependency_matrix_[i][symbols_.size()+independent_parameters_[j]]);
        }
    }
    p_new[not_zero[0]]=tmp;
 }
}

void Model::print_original_parameters(std::ostream &os, std::vector<double> &p) {
    assert( p.size()==symbols_.size());
    for (size_t i=0; i<symbols_.size(); i++) {
        os << symbols_[i] << "\t" << p[i] << std::endl;
    }
}

void Model::convert_original_parameter_to_response_matrix( 
    double_matrix &d, 
    std::vector<double> &inh, 
    const std::vector<double> &p, 
    const int_matrix &adj) 
{
    size_t counter=0;
    double_matrix::extent_gen extents;
    d.resize(extents[adj.shape()[0]][adj.shape()[1]]);
    size_t size=adj.shape()[0];
    for (size_t i=0; i<size; i++){
        for (size_t j=0; j<size; j++){
            if ((i!=j)&&(adj[j][i]!=0)) {
                assert(counter<p.size());
                d[j][i]=p[counter++];
            } else {
                d[j][i]=0;
            }
        }
    }
    inh.clear();
    for (; counter < p.size(); counter++) {
        inh.push_back(p[counter]);
    }
}


void Model::convert_response_matrix_to_original_parameter(std::vector<double> &p, const double_matrix &d, const std::vector<double> &inh, const int_matrix &adj) {
    p.clear();
    size_t size=adj.shape()[0];
    for (size_t i=0; i<size; i++){
        for (size_t j=0; j<size; j++){
            if ((i!=j)&&(adj[j][i]!=0)) {
    p.push_back(d[j][i]);
            }
        }
    }
    for (size_t i=0; i<inh.size(); ++i) {
        p.push_back(inh[i]);
    }
}
