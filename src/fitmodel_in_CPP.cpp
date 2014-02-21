// call as make fitmodel_in_CPP && ./fitmodel_in_CPP

#include "fitmodel_in_CPP.hpp"
#include "identifiability.hpp"

//#include <nlopt.hpp>

#include <vector>                       
#include <fstream> 
#include <boost/numeric/ublas/io.hpp>   
#include <boost/bind.hpp>   
#include <boost/math/distributions/chi_squared.hpp>
#include <levmar-2.5/levmar.h>               // to include Marquardt-Levenberg fitting (lsqnonlin in matlab)(written for C)
#include <float.h>                       // needed to declare the lower and upper bounds to plus and mimnus infinity DBL_MAX and -DBL_MAX
//#include<cstring>    // include these two if you want to view some ginac expressions in matlab (mex)
//#include<string>

#include <math.h>
#include <stdlib.h>


#include "helper_types.hpp"

extern int verbosity;

// convenience function to print GiNaC::symbolic vectors
std::ostream & operator<<(std::ostream &os, std::vector<GiNaC::symbol> &v) {             
  for (std::vector<GiNaC::symbol>::iterator iter=v.begin(); iter!=v.end(); ++iter) 
    os << *iter << std::endl ;
  return os;
}

// This class wraps all stuff that the optimiser needs to know about the model:
// i.e. the model itself, which parameters to vary, and which are constant
// since the optimizer thinks that all parmeters can be optimized, we have to
// convert from a vector of variable parameters to all parameters and back.
class levmar_pass_info {
public: 
  levmar_pass_info(const Model * model, std::vector<size_t> fixed_params, std::vector<double> params, const Data * data) : 
    model_(model), fixed_params_(fixed_params), params_(params), data_(data) {
    data_->data_consistent(model_->exp_design());
    if (fixed_params_.size()==0) {
      direct_pass=true; 
    } else {
      direct_pass=false;
      // some sanity checks.
      if (params_.size()!=model_->nr_of_parameters()) {
	std::cerr << "Problem: params is not of same size as model params! aborting!!!" << std::endl; exit(-1); 
      }
      mask_=std::vector<bool>(model_->nr_of_parameters(),false);
      for (std::vector<size_t>::iterator iter=fixed_params_.begin(); iter!=fixed_params_.end(); iter++) {
	if (*iter>=model_->nr_of_parameters()) {
	  std::cerr << "Problem: some non-existing parameters have been declared fixed. aborting!!!" << std::endl; exit(-1); 
	}
	mask_[*iter]=true;
      }
    }
  }

  void eval(const double *p,double *datax, int m, int n) const {
    if (direct_pass) 
      model_->eval(p,datax,data_); 
    else
      {
	double pall[model_->nr_of_parameters()];
	size_t j=0;
	for (size_t i=0; i<model_->nr_of_parameters(); i++) {
	  if (mask_[i]) {
	    pall[i]=params_[i];
	  } else {
	    pall[i]=p[j++];
	    if (j>=(size_t)m) { std::cerr << "Problem: some problem in parameter sorting. aborting!!!" << std::endl; exit(-1); }
	  }
	}
	model_->eval(pall,datax,data_);
      }
  }

  // These two functions trust that the sizes of arrays are correct !
  void convert_to_variable_params(double *pvar, double *pall) {
    size_t j=0;
    for (size_t i=0; i<model_->nr_of_parameters(); i++) {
      if (!mask_[i]) pvar[j++]=pall[i];
    }
    
  }

  // These two functions trust that the sizes of arrays are correct !
  void convert_from_variable_params(double *pall, double *pvar) {
    size_t j=0;
    for (size_t i=0; i<model_->nr_of_parameters(); i++) {
      if (mask_[i]) {
	pall[i]=params_[i];
      } else {
	pall[i]=pvar[j++];
      }
    }
  }


protected:
  const Model * model_;
  const Data * data_;
  std::vector<size_t> fixed_params_;
  std::vector<double> params_;
  std::vector<bool> mask_;
  bool direct_pass;
};

 // Fitting routine for the dlevmar_dif function
void levmar_wrapper(double *p,double *datax,int m, int n, void *pass_data_tmp) { 

  const levmar_pass_info * lpi = (const levmar_pass_info *) pass_data_tmp;
  lpi->eval(p, datax, m,  n);  

}

double fit_using_lsqnonlin(const Model * model, double *datax, size_t number_of_measurements, double *param, std::vector<size_t> keep_constant, const Data * data) {

  double info[LM_INFO_SZ]; // returns some information about the actual fitting process (info[1])=sumsq
  double opts[LM_OPTS_SZ]; // some options to determine initial  \mu (opts[0]), stopping thresholds for some functions (opts[1-3]) and difference approxiamtion to Jacobian)
  opts[0]=LM_INIT_MU*10.0;
  opts[1]=1E-10;
  opts[2]=1E-10;
  opts[3]=1E-15; 
  opts[4]=LM_DIFF_DELTA/10.0;//relevant only if the finite differences Jacobian version is used  

 

  std::vector<double> parameters_fixed;

  if (keep_constant.size()>0) {
    parameters_fixed.resize(model->nr_of_parameters());
    for (size_t i=0; i<model->nr_of_parameters(); i++) 
      parameters_fixed[i]=param[i];
  }
  levmar_pass_info lpi(model, keep_constant, parameters_fixed, data);
  double parameters[model->nr_of_parameters()];
  if (keep_constant.size()>0) {
    lpi.convert_to_variable_params(parameters,param);
    dlevmar_dif(levmar_wrapper,parameters,datax,model->nr_of_parameters(),number_of_measurements,10000,opts,info,NULL,NULL,(void *)&lpi);
    lpi.convert_from_variable_params(param,parameters);
  } else {
     dlevmar_dif(levmar_wrapper,param,datax,model->nr_of_parameters(),number_of_measurements,10000,opts,info,NULL,NULL,(void *)&lpi);
  }   

  std::string termination_why;
  switch ((int)info[6]) {
  case 1: termination_why="stopped by small gradient J^T";
    break;
  case 2: termination_why="stopped by small Dp";
    break;
  case 3: termination_why="stopped by itmax";
    break;
  case 4: termination_why="singular matrix. Restart from current p with increased mu";
    break;
  case 5: termination_why="no further error reduction is possible. Restart with increased mu";
    break;
  case 6: termination_why="stopped by small ||e||_2";
    break;
  case 7: termination_why="stopped by invalid (i.e. NaN or Inf) func values; a user error";
    break;
  default: 
    termination_why="Dont know why terminated!";
    break;
  }
  if (verbosity>1) 
    std::cerr << "Terminated optimisation after " << info[5] << 
      " iterations and "<< info[7]<< 
      " function evaluations because : " << 
      termination_why << std::endl;
  
   return info[1];
}


void fitmodel( std::vector <double> &bestfit,
	       double * bestresid,
	       double_matrix &prediction,
	       const Model *model,
	       const Data  *data, 
	       std::vector<size_t> keep_constant
	       ) {
  
  size_t number_of_parameters=model->nr_of_parameters();
  
  data->data_consistent(model->exp_design());
  double p[number_of_parameters];
  
  // starting parameter value
  if (bestfit.size()==number_of_parameters) {
    if (verbosity>1) std::cerr << "Use existing parameter vector: " ;
    for (size_t tmp=0; tmp<number_of_parameters; tmp++) { p[tmp]=bestfit[tmp]; if (verbosity>1) std::cerr << p[tmp] << ", "; }
    if (verbosity>1) std::cerr << std::endl;
  } else {
    for (size_t tmp=0; tmp<number_of_parameters; tmp++) p[tmp]=0.39;            
  }

  //  define the measurement value to compare with measured value divided by the error
  int number_of_measurements=data->stim_data.shape()[1]*
    data->stim_data.shape()[0]; 

  double datax[number_of_measurements];
  for (size_t i=0; i<data->stim_data.shape()[1]; i++ ) {
    for (size_t j=0; j<data->stim_data.shape()[0]; j++) {
      if (std::isnan(data->error[j][i])) {
	datax[i*data->stim_data.shape()[0]+j]=0;
      } else {
	datax[i*data->stim_data.shape()[0]+j]=data->stim_data[j][i]/data->error[j][i];
      }
    }
  }
  assert(data->stim_data.shape()[1]*data->stim_data.shape()[0] == (size_t)number_of_measurements);

  if (model->modelRank()!=model->nr_of_parameters()) {
    *bestresid=1000000000000000.0;
    return;
  }


  *bestresid=fit_using_lsqnonlin(model, datax, number_of_measurements, p, keep_constant, data);

  bestfit.resize(number_of_parameters);
  std::copy(p,p+number_of_parameters,bestfit.begin());

   //   model.print_parameter_report(std::cerr, bestfitvectorfuerausgabe);
  
  //  std::cerr << "bestresid "<< *bestresid << std::endl;

   prediction.resize(boost::extents[data->error.shape()[0]][data->error.shape()[1]]);
   double datay[number_of_measurements];
   model->eval(p,datay,data);
   for (size_t i=0;i<data->error.shape()[1];i++){
     for (size_t j=0;j<data->error.shape()[0];j++){
       prediction[j][i]=datay[i*(data->error.shape()[0])+j]*data->error[j][i];
     }
   }

}

// Computes the profile likelihood and the variation of the other parameters depending on the variation of one parameter
void profile_likelihood(const Data data,
	std::vector<double> parameters,
	const std::vector<size_t> keep_constant,
	std::vector<std::vector<double>> &residual_track,
	std::vector<double> &explored,
	const double param_value,
	const Model *model,
	bool* n_identifiability,
    std::vector<double> &thresholds,
	const unsigned int total_steps = 10000,
	const double step_size = 0.01) {

    double residual;
    double_matrix prediction;

    // Initial fit
    fitmodel(parameters, &residual, prediction, model, &data);
    double old_params = parameters[keep_constant[0]];

    // Thresholds, repeated twice
    thresholds.push_back(residual + boost::math::quantile( boost::math::chi_squared(1), 0.95 ));
    thresholds.push_back(residual + boost::math::quantile( boost::math::chi_squared(parameters.size()), 0.95 ));
    
    // Upper and Lower scans are separated, to be sure to scan 
    double scanned_value = param_value;
    // Lower values scan
	std::vector<std::vector<double>> dec_residual;
    for (int i=0 ; i < parameters.size() ; i++) {
        dec_residual.push_back(std::vector<double>);
    }
	std::vector<double> dec_explored;
    for (unsigned int i=1 ; i < total_steps / 2 ; i++) {
        parameters[keep_constant[0]] = scanned_value;
        scanned_value -= step_size;

        fitmodel(parameters, &residual, prediction, model, &data, keep_constant);
        dec_explored.push_back(scanned_value);
        for (int j=0 ; j < parameters.size() ; j++) {
            dec_residual[j].push_back(parameters[j]);
        }
        // In the row corresponding to the parameter's values, we put the residual
        dec_residual[keep_constant].pop_back(residual);
        dec_residual[keep_constant].push_back(residual);

        if (residual > thresholds[0]) n_identifiability[2] = false;
        if (residual > thresholds[1]) n_identifiability[3] = false;
    }
    // We reorder the scanned values in the final vector
    int size = dec_explored.size();
    for (int it=0 ; it <  size; it++) {
        explored.push_back(dec_explored[size -1 -it]);
        residual_track.push_back(dec_residual[size -1 -it]);
    }

    // Upper values scan
    scanned_value = param_value;
    for (unsigned int i=0 ; i < total_steps / 2 ; i++) {
        parameters[keep_constant[0]] = scanned_value;
        scanned_value += step_size;

        fitmodel(parameters, &residual, prediction, model, &data, keep_constant);
        explored.push_back(scanned_value);
        residual_track.push_back(residual);

        if (residual > thresholds[0]) n_identifiability[0] = false;
        if (residual > thresholds[1]) n_identifiability[1] = false;
    }

}

