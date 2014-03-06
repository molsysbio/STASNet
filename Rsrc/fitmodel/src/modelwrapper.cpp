#include "modelwrapper.hpp"

#include <fitmodel/model.hpp>
#include <fitmodel/generate_response.hpp>
#include <fitmodel/fitmodel_in_CPP.hpp>
#include <Rcpp.h>
//#include <boost/thread>

ModelWrapper::ModelWrapper() : model(NULL), linear_approximation(FALSE) { }

ModelWrapper::~ModelWrapper() {
  if (model!=NULL) delete model;
}

void ModelWrapper::printResponse() {

    if (model == NULL) throw(std::logic_error("Model not set"));
    model-> printResponse();
} 

size_t ModelWrapper::modelRank() {
    if (model == NULL) throw(std::logic_error("Model not set"));
    else return model->modelRank();
    return -1;
} 

size_t ModelWrapper::nr_of_parameters() {
  if (model == NULL) throw(std::logic_error("Model not set"));
  else return model->nr_of_parameters();
  return -1;
} 

  
bool ModelWrapper::model_design_consistent(ExperimentalDesign &exp, ModelStructure &mod) {
  int nodes= mod.getAdjacencyMatrix().shape()[0]; 
  if (nodes!=mod.getAdjacencyMatrix().shape()[0]) 
    throw (std::invalid_argument("adj matrix not square"));
  if (nodes!=mod.getNames().size()) 
    throw (std::invalid_argument("names size wrong"));
  if (exp.basal_activity.size() != nodes) 
    throw (std::invalid_argument("basal activity size wrong"));
  if (exp.inhib_nodes.size() != exp.inhibitor.shape()[1]) 
    throw (std::invalid_argument("inhib_nodes and inhibitor don't match"));
  if (exp.stim_nodes.size() != exp.stimuli.shape()[1]) 
    throw (std::invalid_argument("inhib_nodes and inhibitor don't match"));
  if (exp.stimuli.shape()[0] != exp.inhibitor.shape()[0]) 
    throw (std::invalid_argument("number of stimuli and inhibitors don't match"));
    return true;
}

void ModelWrapper::setModel(ExperimentalDesign exp, ModelStructure mod) {
  model_design_consistent(exp,mod);

  std::cout << mod; // DEBUGGING  
  generate_response(response_full_model,  
		      symbols_full_model,
		      mod.getAdjacencyMatrix(),
		      exp,
		      mod.getNames());
  adjacency_matrix.resize(boost::extents[mod.getAdjacencyMatrix().shape()[0]]
			  [mod.getAdjacencyMatrix().shape()[1]]);
  adjacency_matrix=mod.getAdjacencyMatrix();
  
  if (model != NULL) delete model;
  model = new Model(response_full_model,  
		    symbols_full_model,
		    exp, linear_approximation);
} 

SEXP ModelWrapper::simulate(Data data, std::vector<double> parameters) {

  if ( model == NULL ) 
    throw std::logic_error("Model not initialized yet. use setModel() before");

  if ( parameters.size() != model->nr_of_parameters() ) 
    throw std::invalid_argument("length of parameter vector invalid");

  double_matrix datax;
  model->predict(&(parameters[0]), datax, &data);
  Rcpp::List ret;
  ret["prediction"]=datax;

  return ret;
}

SEXP ModelWrapper::fitmodel(Data data, std::vector<double> parameters)  { 

  if ( parameters.size() != model->nr_of_parameters() ) 
    throw std::invalid_argument("length of parameter vector invalid");

    double residual;
    double_matrix  predictions;
    ::fitmodel(parameters, &residual, predictions, model, &data);
    Rcpp::List ret;
    Rcpp::NumericVector pars( parameters.begin(), parameters.end() );
    ret["parameter"]=pars;
    ret["residuals"]=residual;
    return ret;
}

// Computes the profile likelihood of one parameter and the functionnal relationships of the others with this parameter
SEXP ModelWrapper::profileLikelihood(Data data, std::vector<double> parameters, size_t target, const unsigned int total_steps = 10000, const double step_size = 0.01) {
    if ( parameters.size() != model->nr_of_parameters() ) 
        throw std::invalid_argument("length of parameter vector invalid");

    double param_value = parameters[target-1];
    std::vector<size_t> keep_constant(1, target-1); // -1 for R users
    std::vector< std::vector<double> > residual_track;
    std::vector<double> explored;
    bool identifiability[4] = {false};
    std::vector<double> thresholds;

    ::profile_likelihood( data, parameters, keep_constant, residual_track, explored, param_value, model, identifiability, thresholds, total_steps, step_size);
    
    Rcpp::List ret;
    Rcpp::NumericMatrix track(parameters.size(), explored.size());
    for (int i=0 ; i < parameters.size() ; i++) {
        for (int j=0 ; j < explored.size() ; j++) {
            track(i,j) = residual_track[i][j];
        }
    }
    std::vector<std::string> paths;
    model->getParametersLinks(paths);

    ret["residuals"] = track;
    ret["explored"] = explored;
    ret["lowt_upper"] = identifiability[2];
    ret["lowt_lower"] = identifiability[0];
    ret["hight_upper"] = identifiability[3];
    ret["hight_lower"] = identifiability[1];
    ret["thresholds"] = thresholds;
    ret["path"] = paths[target-1];
    ret["pathid"] = target;
    ret["value"] = param_value;

    return ret;

}

/* COMMENTED BECAUSE IT NEEDS C++11 LIBRARY thread
 * IF YOU UNCOMMENT THIS, UNCOMMENT THE THREAD HEADERS AS WELL
// Returns the profile likelihood and functionnal relationship for all parameters, does the computation in parallel
SEXP ModelWrapper::parallelPL(Data data, std::vector<double> parameters, const unsigned int total_steps = 10000, const double step_size = 0.01) {
    if ( parameters.size() != model->nr_of_parameters() ) 
        throw std::invalid_argument("length of parameter vector invalid");

    Rcpp::List *returned;

    std::vector< std::vector<size_t> > keep_constant;
    std::vector< std::vector< std::vector<double> > > residual_track;
    std::vector< std::vector<double> > explored;
    std::vector< bool* > identifiability;
    std::vector< std::vector<double> > thresholds;
    std::vector<size_t> target;

    // Creation of the threads
    std:vector<std::thread> threads;
    for (int i=0 ; i<parameters.size() ; i++) {
        keep_constant.push_back(std::vector<size_t>(1, i));
        residual_track.push_back(std::vector< std::vector<double> >);
        explored.push_back(std::vector<double>);
        identifiability.push_back(new bool[4]);
        thresholds.push_back(std::vector<double>);
        target.push_back(i);

        // Create a new thread for each parameter
        threads.push_back(std::thread(::profile_likelihood, data, parameters, keep_constant[i], residual_track[i], explored[i], parameters[i], model, identifiability[i], thresholds[i], total_steps, step_size));
        Rcpp::NumericMatrix track(parameters.size(), explored.size());
        for (int k=0 ; k < parameters.size() ; k++) {
            for (int j=0 ; j < explored.size() ; j++) {
                track[i](k, j) = residual_track[i][k][j];
            }
        }
        Rcpp::CharacterVector paths = getParametersLinks();

        returned[i]["residuals"] = track[i];
        returned[i]["explored"] = explored[i];
        returned[i]["lowt_upper"] = identifiability[i][2];
        returned[i]["lowt_lower"] = identifiability[i][0];
        returned[i]["hight_upper"] = identifiability[i][3];
        returned[i]["hight_lower"] = identifiability[i][1];
        returned[i]["thresholds"] = thresholds[i];
        ret["pathid"] = target[i];

        delete[] identifiability[i];

    }

    Rcpp::List ret(parameters.size());
    for (int i=0 ; i<parameters.size() ; i++) {
        threads[i].join();
        ret[i] = returned[i];
    }

    return ret;

}
*/


SEXP ModelWrapper::getLocalResponse( std::vector<double> p ) {
  
    if ( model == NULL ) 
        throw std::logic_error("Model not initialized yet. use setModel() before");
    if ( p.size() != model->nr_of_parameters() ) 
        throw std::invalid_argument("length of parameter vector invalid");

    double_matrix fit_response_matrix; 
    std::vector<double> fit_inh;
    std::vector<double> p_new;
    {
      model->convert_identifiables_to_original_parameter( p_new, p );
      Model::convert_original_parameter_to_response_matrix( fit_response_matrix, fit_inh, p_new, adjacency_matrix );
    }
    Rcpp::List ret;
    ret["local_response"]=fit_response_matrix;
    ret["inhibitors"]=fit_inh;
    return ret;
    
}


std::vector<double> ModelWrapper::getParameterFromLocalResponse( const double_matrix &response, const std::vector<double> inhib) {
  
  if ( model == NULL ) 
    throw std::logic_error("Model not initialized yet. use setModel() before");
  
  std::vector<double> parameter;
  std::vector<double> tmpp;
  
  Model::convert_response_matrix_to_original_parameter( tmpp, response, 
							inhib, adjacency_matrix );
  model->convert_original_parameter_into_identifiable(parameter,tmpp);

  return parameter;
  
}

// Gives the links combination for each identifiable response coefficient
SEXP ModelWrapper::getParametersLinks() {

    std::vector<std::string> buf_string;

    model->getParametersLinks(buf_string);
    Rcpp::CharacterVector ret(buf_string.size());
    for (int i=0 ; i < buf_string.size() ; i++) {
        ret(i) = buf_string[i];
    }
   
    return ret;
}

// Show the parameter dependency matrix G reduced
void ModelWrapper::showParameterDependencyMatrix() {
    model->showParameterDependencyMatrix();
}


RCPP_MODULE(ModelEx) {
  using namespace Rcpp ;
  class_<ModelWrapper>( "Model" )
    .default_constructor()
    .method( "setModel", &ModelWrapper::setModel )
    .method( "simulate", &ModelWrapper::simulate )

    .method( "printResponse", &ModelWrapper::printResponse )
    .method( "modelRank", &ModelWrapper::modelRank )
    .method( "nr_of_parameters", &ModelWrapper::nr_of_parameters )
    .method( "fitmodel", &ModelWrapper::fitmodel )
    .method( "getLocalResponseFromParameter", &ModelWrapper::getLocalResponse )
    .method( "getParameterFromLocalResponse", &ModelWrapper::getParameterFromLocalResponse )
    .method( "profileLikelihood", &ModelWrapper::profileLikelihood)
    //.method( "parallelPL", &parallelPL )
    .method( "getParametersLinks", &ModelWrapper::getParametersLinks )
    .method( "showParameterDependencyMatrix", &ModelWrapper::showParameterDependencyMatrix )
    .field("linear_approximation", &ModelWrapper::linear_approximation, "Linear Approximation" )
    ;
}

int verbosity=0;
