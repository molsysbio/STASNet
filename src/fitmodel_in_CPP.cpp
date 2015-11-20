// call as make fitmodel_in_CPP && ./fitmodel_in_CPP

#include "fitmodel_in_CPP.hpp"
#include "identifiability.hpp"

//#include <nlopt.hpp>

#include <vector>                       
#include <fstream> 
#include <exception>
#include <boost/numeric/ublas/io.hpp>   
#include <boost/bind.hpp>   
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>
#include "levmar-2.5/levmar.h"               // to include Marquardt-Levenberg fitting (lsqnonlin in matlab)(written for C)
#include <float.h>                       // needed to declare the lower and upper bounds to plus and mimnus infinity DBL_MAX and -DBL_MAX
//#include<cstring>    // include these two if you want to view some ginac expressions in matlab (mex)
//#include<string>

#include <math.h>
#include <stdlib.h>


#include "helper_types.hpp"

extern int verbosity;
extern bool debug;

// convenience function to print GiNaC::symbolic vectors
std::ostream & operator<<(std::ostream &os, std::vector<GiNaC::symbol> &v) {             
  for (std::vector<GiNaC::symbol>::iterator iter=v.begin(); iter!=v.end(); ++iter) 
    os << *iter << std::endl ;
  return os;
}

// This class wraps all stuff that the optimiser needs to know about the model:
// i.e. the model itself, which parameters to vary, and which are constant
// since the optimizer thinks that all parameters can be optimized, we have to
// convert from a vector of variable parameters to all parameters and back.
class levmar_pass_info {
public: 
  levmar_pass_info(const Model * model, std::vector<size_t> fixed_params, std::vector<double> params, const Data * data) : 
    model_(model), data_(data), fixed_params_(fixed_params), params_(params) {
    if (!data_->data_consistent(model_->exp_design()))
      throw std::logic_error("Error: Data and Model do not match");

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
    else {
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
  double opts[LM_OPTS_SZ]; // some options to determine initial  \mu (opts[0]), stopping thresholds for some functions (opts[1-3]) and difference approximation to Jacobian)
  opts[0]=LM_INIT_MU*100.0;
  opts[1]=1E-30; // Gradient value
  opts[2]=1E-200; // Parameters variation
  opts[3]=1E-17; // Residual variation
  opts[4]=LM_DIFF_DELTA/10.0;//relevant only if the finite differences Jacobian version is used
  unsigned int max_steps = 10000;

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
    dlevmar_dif(levmar_wrapper,parameters,datax,model->nr_of_parameters(),number_of_measurements,max_steps,opts,info,NULL,NULL,(void *)&lpi);
    lpi.convert_from_variable_params(param,parameters);
  } else {
    dlevmar_dif(levmar_wrapper,param,datax,model->nr_of_parameters(),number_of_measurements,max_steps,opts,info,NULL,NULL,(void *)&lpi);
  }   

  std::string termination_why;
  switch ((int)info[6]) {
  case 1: termination_why="stopped by small gradient J^T"; // Ideally, the reason why terminated
    break;
  case 2: termination_why="stopped by small Dp"; // opts[2] to big
    break;
  case 3: termination_why="stopped by itmax"; // On the good way but from too far
    break;
  case 4: termination_why="singular matrix. Restart from current p with increased mu";
    break;
  case 5: termination_why="no further error reduction is possible. Restart with increased mu"; //
    break;
  case 6: termination_why="stopped by small ||e||_2";
    break;
  case 7: termination_why="stopped by invalid (i.e. NaN or Inf) func values; a user error"; // Bad initialisation, it's ok to have some
    break;
  default: 
    termination_why="Dont know why terminated!";
    break;
  }
  if (verbosity>1) {
    std::cerr << "Terminated optimisation after " << info[5] << 
      " iterations and "<< info[7]<< 
      " function evaluations with a residual of " << info[1] << ", because : " << 
      termination_why << std::endl;
      }
  
   return info[1];
}

void fitmodel( std::vector <double> &bestfit,
         double * bestresid,
         double_matrix &prediction,
         const ModelSet *model,
         const Data **data, 
         std::vector<size_t> keep_constant
         ) {
         // Convert a set of Data matrices to a unique Data matrix to call the Model fitmodel function

         for (size_t mod = 1 ; mod < model->getNbModels() ; mod++) {
            assert(data[mod]->stim_data.shape()[0] == data[0]->stim_data.shape()[0]);
            assert(data[mod]->stim_data.shape()[1] == data[1]->stim_data.shape()[1]);
         }

        double_matrix stim_data(boost::extents[model->getNbModels()*data[0]->stim_data.shape()[0]][data[0]->stim_data.shape()[1]]);
        double_matrix unstim_data;
        double_matrix error;
        double_matrix scale;
        boost::multi_array_types::index_gen indices; // To generate the views
        typedef boost::multi_array_types::index_range index_range;
        for (size_t mod=0 ; mod < model->getNbModels() ; mod++) {
            stim_data[ indices[index_range(mod * data[mod]->stim_data.shape()[0], (mod+1) * data[mod]->stim_data.shape()[0])][index_range(0, stim_data.shape()[1])] ] = data[mod]->stim_data;
            unstim_data[ indices[index_range(mod * data[mod]->unstim_data.shape()[0], (mod+1) * data[mod]->unstim_data.shape()[0])][index_range(0, unstim_data.shape()[1])] ] = data[mod]->unstim_data;
            scale[ indices[index_range(mod * data[mod]->scale.shape()[0], (mod+1) * data[mod]->scale.shape()[0])][index_range(0, scale.shape()[1])] ] = data[mod]->scale;
            error[ indices[index_range(mod * data[mod]->error.shape()[0], (mod+1) * data[mod]->error.shape()[0])][index_range(0, error.shape()[1])] ] = data[mod]->error;
        }
        Data* fitdata = new Data();
        fitdata->setStimData(stim_data);
        fitdata->setScale(scale);
        fitdata->setUnstimData(unstim_data);
        fitdata->setError(error);

        fitmodel(bestfit, bestresid, prediction, model, fitdata, keep_constant);
}

void fitmodel( std::vector <double> &bestfit,
         double * bestresid,
         double_matrix &prediction,
         const Model *model,
         const Data *data, 
         std::vector<size_t> keep_constant
         ) {
  
  size_t number_of_parameters=model->nr_of_parameters();
  
  if (!data->data_consistent(model->exp_design())) 
    throw std::logic_error("Error: Data and Model do not match");
  double p[number_of_parameters];
  
  // starting parameter value
  if (bestfit.size()==number_of_parameters) {
    if (verbosity>10) std::cerr << "Use existing parameter vector: ";
    for (size_t tmp=0; tmp<number_of_parameters; tmp++) {
        p[tmp]=bestfit[tmp];
        if (verbosity>10) std::cerr << p[tmp] << ", ";
    }
    if (verbosity>10) std::cerr << std::endl;
  } else {
    for (size_t tmp=0; tmp<number_of_parameters; tmp++) p[tmp]=0.39;            
  }

  //  define the measurement value to compare with simulated values divided by the error
  int number_of_measurements=data->stim_data.shape()[1] * data->stim_data.shape()[0]; 
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

  // Can't fit if the system is non identifiable
  if (model->modelRank()!=model->nr_of_parameters()) {
    *bestresid=1000000000000000.0;
    std::cout << "The system is non identifiable, you need less links or more conditions" << std::endl;
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

// Perform simulated annealing to find the global optimum for the model
void simulated_annealing(const Model *model, const Data *data, std::vector<double> &bestfit, double &bestresid, int max_it, int max_depth) {
     // Initialise the random generator, as there is no main function to do it, problem for parallelised version
    std::srand(std::time(NULL)*getpid()); // WARNING, getpid is only provided for Unix environments

    // Parameters with a suspected value will be reseted to this value between each fit
    std::vector<double> fixed_params;
    std::vector<int> fixed_index, var_index;
    if (bestfit.size() == model->nr_of_parameters()) {
        for (size_t i=0 ; i < bestfit.size() ; i++) {
            if (bestfit[i] != 0) {
                fixed_index.push_back(i);
                fixed_params.push_back(bestfit[i]);
            } else {
                var_index.push_back(i);
            }
        }
    } else {
        bestfit.resize(model->nr_of_parameters());
    }
    std::cerr << "Annealing setup..." << std::endl;
    if (max_it == 0) {
        max_it = pow(2, model->nr_of_parameters());
    }
    if (max_depth == 0) {
        max_depth = model->nr_of_parameters();
    }

    // Find the order of magnitude of the steps and calculate the parameters from it
    size_t hypercube_samples = 100;
    double exploration[hypercube_samples];
    std::vector<size_t> index(hypercube_samples);
    std::vector< std::vector<double> > parameters = normalLHS(hypercube_samples, model->nr_of_parameters(), 1);
    double test_params[model->nr_of_parameters()];
    std::cout << "Annealing init steps :" << std::endl;
    for (size_t i=0 ; i < hypercube_samples ; i++) {
        index.push_back(i);
        std::copy(parameters[i].begin(), parameters[i].end(), test_params);
        //exploration[i] = model->score(static_cast<double*>(parameters[i].begin()), data);
        exploration[i] = model->score(test_params, data);
        //std::cout << exploration[i] << std::endl;
        //std::cout << test_params[0] << " ";
    }
    std::sort(index.begin(), index.end(), rev_index_cmp<double*>(exploration));
    std::copy(parameters[index[0]].begin(), parameters[index[0]].end(), test_params);
    std::copy(test_params, test_params + model->nr_of_parameters(), bestfit.begin());
    double temperature = exploration[index[std::floor(hypercube_samples/2)]] - exploration[index[0]];
    const double threshold = temperature/10;
    double cooling = std::exp(std::log(threshold/temperature) / max_depth);

    // Simulated annealing procedure
    unsigned int no_change = 0;
    double p[model->nr_of_parameters()];
    double new_p[model->nr_of_parameters()];
    double residual = exploration[index[0]], new_residual;
    double min_new = 0; //debug
    int depth = 0;
    bestresid = residual;
    while (temperature > threshold && no_change < max_it) {
        depth ++;
        std::vector< std::vector<double> > random_steps = normalLHS(max_it, var_index.size(), 1);
        do {
            // New step reset the fixed parameters and make a step on the other dimensions
            // Try a different step from the current best
            for (size_t i=0 ; i < fixed_index.size() ; i++) {
                new_p[fixed_index[i]] = p[index[i]];
            }
            for (size_t i=0 ; i < var_index.size() ; i++) {
                new_p[var_index[i]] = p[index[i]] + random_steps[no_change][i];
            }
            new_residual = model->score(new_p, data);
            min_new = std::min(new_residual, min_new);
            no_change++;
        } while (no_change < max_it && uniform_sampling() > std::exp((bestresid - new_residual)/temperature));
        if (debug) {
            std::cout << "best = " << bestresid << " old = " << residual << " min_new = " << min_new << ", probablility = " << std::exp((bestresid - min_new)/temperature)  << std::endl;
            std::cout << "Temperature " << temperature << ", Unity probability " << std::exp(-1/temperature) << std::endl;
        }
        if (no_change < max_it) {
            no_change = 0;
            temperature *= cooling;
            if (verbosity > 8) {
                std::cout << "Jump to " << new_residual << std::endl;
                std::cout << "Transition probability = " << std::exp((bestresid - new_residual)/temperature) << std::endl;
            }
        }
        // Collect the new position and residual
        residual = new_residual;
        for (size_t i=0 ; i < var_index.size() ; i++) {
            p[var_index[i]] = new_p[var_index[i]];
        }
        // Update the bestfit if a better solution is found
        if (residual < bestresid) {
            bestresid = residual;
            std::copy(p, p+model->nr_of_parameters(), bestfit.begin());
        }
    }
    for (size_t i=0 ; i < bestfit.size() ; i++) {
        std::cerr << bestfit[i] << ", ";
    }
    std::cout << std::endl;
    if (no_change >= max_it) {
        printf("Simulated annealing terminated at depth %i because no new valid move has been found.\n", depth);
    } else {
        printf("Simulated annealing cooling completed.\n");
        std::cout << "Best residual = " << bestresid << std::endl;
    }
}

// Computes the profile likelihood and the variation of the other parameters depending on the variation of one parameter
void profile_likelihood(const Data &data,
  const std::vector<double> bestfit,
  const std::vector<size_t> keep_constant,
  std::vector< std::vector<double> > &residual_track,
  std::vector<double> &explored,
  const double param_value,
  const Model *model,
  pl_analysis &thresholds,
  const unsigned int total_steps = 10000) {

    double residual;
    double_matrix prediction;

    // Initial fit
    std::vector<double> parameters = bestfit;
    fitmodel(parameters, &residual, prediction, model, &data);
    //    double previous_residual = residual;

    // Add the values of the thresholds
    thresholds.decision = 0.95;
    thresholds.pointwise_threshold = residual + boost::math::quantile( boost::math::chi_squared(1), thresholds.decision );
    thresholds.simultaneous_threshold = residual + boost::math::quantile( boost::math::chi_squared(parameters.size()), thresholds.decision );
    
    // Upper and Lower scans are separated, to be sure to scan near the optimum each time 
    // Lower values scan
    bool first_low = true;
    bool first_high = true;
    double scanned_value = param_value;
    // By default, we explore 1.5 times the parameter with a minimum step size of 0.01
    // ?? except if the parameter is smaller than 0.01 in which case we make sure to have 10 iterations before 0
    double step_size = std::min(-std::abs(parameters[keep_constant[0]]) * 3 / total_steps, -0.01);
    /*double step_size = -0.01;
    step_size = choose_step_size(data, parameters, param_value, keep_constant, model, boost::math::quantile( boost::math::chi_squared(1), decision), residual, step_size);*/
  std::vector< std::vector<double> > dec_residual;
    for (int i=0 ; i < parameters.size() ; i++) {
        dec_residual.push_back(std::vector<double>());
        residual_track.push_back(std::vector<double>());
    }
  std::vector<double> dec_explored;
    for (unsigned int i=0 ; i < std::floor(total_steps / 2) ; i++) {
        parameters[keep_constant[0]] = scanned_value;

        fitmodel(parameters, &residual, prediction, model, &data, keep_constant);
        dec_explored.push_back(scanned_value);
        if (scanned_value == 0 && verbosity > 7) {
            std::cout << "Simulation for value 0, residuals = " << residual << std::endl;
            double_matrix simulation;
            model->predict(parameters, simulation, &data);
            for (unsigned int row=0 ; row < data.unstim_data.shape()[0] ; row++) {
                for (unsigned int col=0 ; col < data.unstim_data.shape()[1] ; col++) {
                    std::cout << "\t" << simulation[row][col];
                }
                std::cout << std::endl;
            }
        }
        
        // We write the other parameters new values
        for (int j=0 ; j < parameters.size() ; j++) {
            dec_residual[j].push_back(parameters[j]);
        }
        // In the row corresponding to the scanned parameter's value, we put the residual
        dec_residual[keep_constant[0]].pop_back();
        dec_residual[keep_constant[0]].push_back(residual);

        // Tell if the thresholds were reached and for which index (in the final array)
        if (residual > thresholds.simultaneous_threshold && first_high) {
            thresholds.hn_threshold = true;
            first_high = false;
            // If the residual reaches both thresholds at once, there is a bad fit at this point and the parameters would have no sense so we use the previous set
            if (first_low) {
                for (int j=0 ; j < parameters.size() ; j++) {
                    parameters[j] = dec_residual[j][i-1];
                }
            }
        }
        if (first_low && residual > thresholds.pointwise_threshold) {
            thresholds.ln_threshold = true;
            first_low = false;
            // If the residual went skyrocket over the threshold, use the index of the previous point below it
            if (residual < thresholds.simultaneous_threshold) {
                thresholds.negative_uncertainty.push_back(std::floor(total_steps/2)-1-i);
            } else {
                thresholds.negative_uncertainty.push_back(std::floor(total_steps/2)-1-i +1);
            }
        } else if (residual < thresholds.pointwise_threshold){
            first_low = true;
        }

        scanned_value += step_size;
        if(scanned_value - step_size < 0 && scanned_value + step_size > 0) {scanned_value = 0;}
    }
    // We reorder the scanned values in the final vectors
    int size = dec_explored.size();
    for (int it=0 ; it < size; it++) {
        explored.push_back(dec_explored[size -1 -it]);
        for (int j=0 ; j < residual_track.size() ; j++) {
            residual_track[j].push_back(dec_residual[j][size -1 -it]);
        }
    }

    // Upper values scan
    first_low = true;
    first_high = true;
    scanned_value = param_value;
    parameters = bestfit;
    step_size = -step_size;
    // By default, we explore 1.5 times the parameter
    //step_size = std::max(std::abs(parameters[keep_constant[0]]) * 3 / total_steps, 0.01);
    //step_size = 0.01;
    //step_size = choose_step_size(data, parameters, param_value, keep_constant, model, boost::math::quantile( boost::math::chi_squared(1), decision), residual, step_size);
    for (unsigned int i=0 ; i < std::floor(total_steps / 2) ; i++) {
        parameters[keep_constant[0]] = scanned_value;

        fitmodel(parameters, &residual, prediction, model, &data, keep_constant);
        explored.push_back(scanned_value);
        if (scanned_value == 0 && verbosity > 7) {
            std::cout << "Simulation for value 0, residuals = " << residual << std::endl;
            double_matrix simulation;
            model->predict(parameters, simulation, &data);
            for (unsigned int row=0 ; row < data.unstim_data.shape()[0] ; row++) {
                for (unsigned int col=0 ; col < data.unstim_data.shape()[1] ; col++) {
                    std::cout << "\t" << simulation[row][col];
                }
                std::cout << std::endl;
            }
        }

        // We write the other parameters new values
        for (int j=0 ; j < parameters.size() ; j++) {
            residual_track[j].push_back(parameters[j]);
        }
        // In the row corresponding to the parameter's values, we put the residual
        residual_track[keep_constant[0]].pop_back();
        residual_track[keep_constant[0]].push_back(residual);

        // Update boolean and index if the thresholds are reached for ascending value
        if (first_high && residual > thresholds.simultaneous_threshold) {
            thresholds.hp_threshold = true;
            first_high = false;
            // If the residual reaches both thresholds at once, there is a bad fit at this point and the parameters would have no sense so we use the previous set
            if (first_low) {
                for (int j=0 ; j < parameters.size() ; j++) {
                    parameters[j] = residual_track[j][i-1];
                }
            }
        }
        if (first_low && residual > thresholds.pointwise_threshold) {
            thresholds.lp_threshold = true;
            first_low = false;
            // Use the previous index if the residual goes skyrocketting
            if (residual > thresholds.simultaneous_threshold) {
                thresholds.positive_uncertainty.push_back(std::floor(total_steps / 2) + i);
            } else {
                thresholds.positive_uncertainty.push_back(std::floor(total_steps / 2) + i -1);
            }
        } else if (residual < thresholds.pointwise_threshold){
            first_low = true;
        }

        scanned_value += step_size;
        if(scanned_value - step_size < 0 && scanned_value + step_size > 0) {scanned_value = 0;}
    }

}

// Determine how large the step should be for the profile likelihood
double choose_step_size(const Data &data,
                        const std::vector<double> parameters,
                        const double param_value,
                        const std::vector<size_t> keep_constant,
                        const Model *model,
                        const double threshold,
                        const double init_residual,
                        double previous_step) {

    const int VARIATION = 10; // Small variation are more precise but take more time
    // To avoid getting stuck on 0
    double step = previous_step;
    if (parameters[keep_constant[0]] == 0) {
        return step;
    }

    const double relative_increase = 0.1;
    // The max step should be bigger for big parameters and not too small for the small ones
    int MAX_STEP = std::abs(param_value) / 10;
    if (MAX_STEP < 10) { MAX_STEP = 10; }

    double residual = init_residual;
    double_matrix prediction;
    std::vector<double> params = parameters;
    params[keep_constant[0]] += step;
    fitmodel(params, &residual, prediction, model, &data, keep_constant);
    double res_dif = std::abs(residual - init_residual);

    // Decrease the step size if it changes to much the residual, but maximise it
    if (res_dif > threshold * relative_increase) {
        while (res_dif > threshold * relative_increase && std::abs(step) > 0.01) {
            step /= VARIATION;
            params = parameters;
            params[keep_constant[0]] += step;

            fitmodel(params, &residual, prediction, model, &data, keep_constant);
            res_dif = std::abs(residual - init_residual);
        }
        return step;
    }
    else {
        while (res_dif < threshold * relative_increase && std::abs(step) < MAX_STEP) {
            step *= VARIATION;
            params = parameters;
            params[keep_constant[0]] += step;

            fitmodel(params, &residual, prediction, model, &data, keep_constant);
            res_dif = std::abs(residual - init_residual);
        }
        if (std::abs(step) >= MAX_STEP) { return previous_step; } // Avoid huge steps for non identifiable parameters
        else { return step / VARIATION; }
        return step / VARIATION;
    }

    /*
    
    double residual = init_residual;
    double previous_residual = init_residual;
    double_matrix prediction;
    std::vector<double> param = parameters;
    double step_size;
    if (previous_step > 0) { step_size = 0.001; }
    else { step_size = -0.001; }

    while (residual <= previous_residual && step_size < 100) {
        step_size *= VARIATION;
        // Fit for small step
        for (int i=0 ; i < VARIATION ; i++) {
            param[keep_constant[0]] += step_size;
            fitmodel(param, &residual, prediction, model, &data, keep_constant);
        }
        previous_residual = residual;

        // Fit for big step
        param = parameters;
        param[keep_constant[0]] = param_value + VARIATION * step_size;
        fitmodel(param, &residual, prediction, model, &data, keep_constant);
    }
    // If we reach big step size, it might be that the parameter is not identifiable and we limit the step size
    if (abs(step_size) >= 100) {
        if (previous_step > 0) { step_size = 1; }
        else { step_size = -1; }
    }
    else if (abs(step_size) < 0.01) { // Minimum step_size
        if (previous_step > 0) { step_size = 0.01; }
        else { step_size = -0.01; }
    }
    if (verbosity > 4) { std::cout << "Default step_size" << std::endl; }

    return step_size;
    */

}

