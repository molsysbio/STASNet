################################ create_model.R #########################################
## Functions that simplify data parsing and model creation and fitting using the package

#' @import Rgraphviz
#' @import pheatmap
#' @import lattice
#' @import parallel
#' @import Rcpp
#' @import lhs
#' @import RhpcBLASctl
#' @useDynLib STASNet

# Global variable to have more outputs
verbose = FALSE
debug = TRUE

get_running_time <- function(init_time, text="") {
  run_time = proc.time()["elapsed"]-init_time
  run_hours = run_time %/% 3600;
  run_minutes = (run_time - 3600 * run_hours) %/% 60;
  run_seconds = round(run_time - 3600 * run_hours - 60 * run_minutes,0);
  return(paste(run_hours, "h", run_minutes, "min", run_seconds, "s", text))
}

trim_num <- function(x, non_zeros=2, behind_comma = 2){
  if (non_zeros==0){
    error("number should have a digits reconsider setting non_zeros larger 0!!")
  }
trim_it <- function(x, non_zeros, behind_comma){  
  if (is.na(x)){ return(x) }

  if (!is.numeric(x)){ oldx =x; x = as.numeric(as.character(x)) } 

  if (is.na(x)){ stop(paste("Number or NA expected '", oldx ,"' received as input!")) } 
  
  if (abs(x >= 1)){ 
    newx = round(x*10^behind_comma)/10^behind_comma
  } else{
    newx =  signif(x,non_zeros)  
  }

  if (nchar(gsub("\\.|-","",as.character(newx))) > max(5,non_zeros+behind_comma)){
    newx = format(newx, scientific = 0)  
  }
  return(newx)
  }

if (is.null(dim(x))){
  return(sapply(x,"trim_it",non_zeros,behind_comma))
}else{
 newx=sapply(1:ncol(x),function(y) sapply(x[,y],"trim_it",non_zeros,behind_comma))
 dimnames(newx) <- dimnames(x)
 return(newx)
}
}

#' Creates a parameterised model from experiment files and the network structure, and fit the parameters to the data
#'
#' @param model_links Path to the file containing the network structure, either in matrix form or in list of links form. Extension .tab expected
#' @param basal_file Path to the file indicating the nodes without basal activity. Extension .dat expected.
#' @param data.stimulation Path to the file containing the data in MRA_MIDAS format. Extension .csv expected.
## CHECK THE IMPLEMENTATION FOR THE NO BASAL
#' @param data.variation Path to the file containing the coefficient of variation for each measurement in MRA_MIDAS format. If it is not provided, the function uses the replicates of the data.stimulation file to determine a variance per probe (i.e antibody/DNA fragment/...). Extension .var expected.
#' @param nb_cores Number of cores that should be used for the computation
#' @param inits Number of initialisation steps which should be performed (see method for the exact meaning of this value)
#' @param perform_plots Whether the distribution of the residuals and the correlation plots for the parameters deduced by correlation should be plotted or not
#' @param precorrelate Whether to infer some links using a linear model instead of a random initialisation
#' @param method Method to be used for the initialisation, available methods are :
#'      random : Perform a Latin Hypercube Sampling to choose \emph{inits} starting points then perform a gradient descent to find the local optimum for each of those points.
#'      correlation : Deduce some parameters from the correlation between the measurements for the target node and all of its input nodes, then perform random to find the other parameters. Recommended, very efficient for small datasets.
#'      genetic : Genetic algorithm with mutation only. \emph{inits} is the total number of points sampled.
#'      annealing : Simulated annealing and gradient descent on the best result. \emph{inits} is the maximum number of iteration without any change before the algorithm decides it reached the best value. Use not recommended.
#' @param unused_perturbations Perturbations in the dataset that should not be used
#' @param unused_readouts Measured nodes in the datasets that should not be used
#' @param MIN_CV Minimum coefficient of variation.
#' @param DEFAULT_CV Default coefficient of variation to use when none is provided and there are no replicated in the data.
#' @param model_name The name of the model is derived from the name of the data.stimulation file name. If data.stimulation is a matrix or a data.frame, 'model_name' will be used to name the model.
#' @param rearrange Whether the rows should be rearranged. "no" to keep the order of the perturbations from the data file, "bystim" to group by stimulations, "byinhib" to group by inhibitions.
#' @return An MRAmodel object describing the model and its best fit, containing the data
#' @export
#' @seealso importModel, exportModel, rebuildModel
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
#' @examples \dontrun{
#' model = createModel("links.tab", "basal.dat", "data_MIDAS.csv") # Produces a model for the network described in links.tab using the data in data_MIDAS.csv
#' model = createModel("links.tab", "basal.dat", "data_MIDAS.csv", "variation.var") # Uses the variation from a variation file
#' model = createModel("links.tab", "basal.dat", "data_MIDAS.csv", nb_cores = detectCores()) # Uses all cores available (with the package parallel)
#' model = createModel("links.tab", "basal.dat", "data_MIDAS.csv", inits = 1000000) # Uses more initialisations for a complex network
#' }
#' @family Model initialisation
# TODO completely remove examples or add datafile so they work (or use data matrices)
createModel <- function(model_links, basal_file, data.stimulation, data.variation="", nb_cores=1, inits=1000, perform_plots=F, precorrelate=T, method="geneticlhs", unused_perturbations=c(), unused_readouts=c(), MIN_CV=0.1, DEFAULT_CV=0.3, model_name="default", rearrange="bystim") {
  # Creation of the model structure object
  model_structure = extractStructure(model_links)
  basal_activity = extractBasalActivity(basal_file)
  
  core = extractModelCore(model_structure, basal_activity, data.stimulation, data.variation, unused_perturbations, unused_readouts, MIN_CV, DEFAULT_CV, rearrange=rearrange)
  expdes = core$design
  data = core$data

  # MODEL SETUP
  model = new(STASNet:::Model)
  model$setModel(expdes, model_structure)
  
  # INITIAL FIT
  results = initModel(model, core, inits, precorrelate, method, nb_cores, perform_plots)
  
  # Choice of the best fit
  params = results$params
  residuals = results$residuals
  fit_info = paste( sum(is.na(residuals)), "NA and ", sum(is.infinite(residuals)), "infinite residuals" )
  message(fit_info)
  # Refit the best residuals to make sure it converged # TODO WITH A MORE PRECISE DESCENT (more max steps and better delta-p)
  order_resid = order(residuals)
  order_id = order_resid[1:min(20, length(residuals))]
  for ( i in 1:length(order_id) ) {
    p_res = residuals[order_id[i]]
    if (!is.na(p_res) && !is.infinite(p_res)) {
      repeat {
        result = model$fitmodel(core$data, params[order_id[i],])
        if (p_res == result$residuals) {
          break
        }
        p_res = result$residuals
      }
      params[order_id[i],] = result$parameter
      residuals[order_id[i]] = result$residuals
    }
  }
  if (debug) {
    # Print the 20 smallest residuals to check if the optimum has been found several times
    message("Best residuals :")
    message(paste0(trim_num(sort(residuals)[1:20],behind_comma = 4), collapse=" "))
  }
  if (perform_plots) {
    plot(1:length(order_resid), residuals[order_resid], main=paste0("Best residuals ", model_name), ylab="Likelihood", xlab="rank", log="y")
  }
  range_var <- function(vv) { rr=range(vv); return( (rr[2]-rr[1])/max(abs(rr)) ) }
  paths = sapply(model$getParametersLinks(), simplify_path_name)
  best_sets = order_resid[signif(residuals[order_resid], 4) == signif(residuals[order_resid[1]], 4)]
  for ( ii in 1:(ncol(params)-1) ) {
    for (jj in (ii+1):ncol(params)) {
      setii = params[best_sets,ii] 
      setjj = params[best_sets,jj]
      cij = suppressWarnings(cor(setii, setjj, use="na.or.complete")) # NA if one vector has a standard deviation of 0
      if (perform_plots) {
          if ((!is.na(cij) && cij > 0.999) || (range_var(setii) > 0.05 && range_var(setjj) > 0.05)) {
            plot(setii, setjj, xlab=paths[ii], ylab=paths[jj], main=paste0("Values for the best fits\nscore=", cij), col=residuals[best_sets])
          }
      }
    }
  }
  
  best_id = order(residuals)[1]
  init_params = params[best_id,]
  init_residual = residuals[best_id]
  if (verbose) {
    message("Model simulation with optimal parameters :")
    message(model$simulate(data, init_params)$prediction)
  }
  
  # Information required to run the model (including the model itself)
  infos = generate_infos(data.stimulation, inits, sort(residuals)[1:5], method, model_links, model_name, fit_info)
  model_description = MRAmodel(model, expdes, model_structure, basal_activity, data, core$cv, init_params, init_residual, name=infos$name, infos=infos$infos, unused_perturbations=unused_perturbations, unused_readouts=unused_readouts, min_cv=MIN_CV, default_cv = DEFAULT_CV)
  message(paste("Residual score =", trim_num(model_description$bestfitscore,behind_comma = 4)))
  
  return(model_description)
}

#' Build an MRAmodelSet
#'
#' Build and fit an MRAmodelSet, which consists of the simultaneous fitting of several MRA models
#' @inheritParams createModel
#' @export
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
createModelSet <- function(model_links, basal_file, csv_files, var_files=c(), nb_cores=1, inits=1000, perform_plots=F, method="geneticlhs", unused_perturbations=c(), unused_readouts=c(), MIN_CV=0.1, DEFAULT_CV=0.3, model_name="default", rearrange="bystim") {
  if (length(csv_files) != length(var_files)) {
    if (length(var_files) == 0) {
      var_files = rep("", length(csv_files))
    } else {
      stop("'var_files' must have the same length as 'csv_files' or be of length 0")
    }
  }
  if (!rearrange %in% c("bystim", "byinhib")) { stop("Invalid 'rearrange' for createModelSet, must be 'byinhib' or 'bystim'") }
  model_structure = extractStructure(model_links)
  basal_activity = extractBasalActivity(basal_file)
  
  nb_submodels = length(csv_files)
  core0 = extractModelCore(model_structure, basal_activity, csv_files[[1]], var_files[[1]], unused_perturbations, dont_read=unused_readouts, MIN_CV, DEFAULT_CV, rearrange=rearrange)
  stim_data = core0$data$stim_data
  unstim_data = core0$data$unstim_data
  error = core0$data$error
  offset = core0$data$scale
  cv = core0$cv
  if (verbose > 8) {
    message("Data dimensions~:")
    message(paste0(dim(core0$data$stim_data),collapse=" "))
    message(paste0(dim(core0$data$unstim_data),collapse=" "))
    message(paste0(dim(core0$data$error),collapse=" "))
  }
  # Build an extended dataset that contains the data of each model
  data_ = new(STASNet:::DataSet)
  data_$addData(core0$data, FALSE)
  for (ii in 2:nb_submodels) {
    core = extractModelCore(model_structure, basal_activity, csv_files[[ii]], var_files[[ii]], unused_perturbations, dont_read=unused_readouts, MIN_CV, DEFAULT_CV, rearrange=rearrange)
    if (!all( dim(core0$data$unstim_data)==dim(core$data$unstim_data) )) {
      stop(paste0("dimension of 'unstim_data' from model ", ii, " do not match those of model 1"))
    } else if (!all( dim(core0$data$error)==dim(core$data$error) )) {
      stop(paste0("dimension of 'error' from model ", ii, " do not match those of model 1"))
    } else if (!all( dim(core0$data$stim_data)==dim(core$data$stim_data) )) {
      stop(paste0("dimension of 'stim_data' from model ", ii, " do not match those of model 1"))
    } else if (!all( dim(core0$cv)==dim(core$cv) )) {
      stop(paste0("dimension of 'cv' from model ", ii, " do not match those of model 1"))
    }
    unstim_data = rbind(unstim_data, core$data$unstim_data)
    stim_data = rbind(stim_data, core$data$stim_data)
    error = rbind(error, core$data$error)
    offset = rbind(offset, core$data$scale)
    data_$addData(core$data, FALSE)
    cv = rbind(cv, core$cv)
  }
  data_$set_stim_data(stim_data)
  data_$set_unstim_data(unstim_data)
  data_$set_error(error)
  data_$set_scale(offset)
  
  model = new(STASNet:::ModelSet)
  model$setModel(core0$design, model_structure)
  model$setNbModels(nb_submodels)
  
  message("setting completed")
  results = initModel(model, list(design=core0$design, data=data_, structure=model_structure), inits, perform_plots=perform_plots, method=method, precorrelate=F, nb_cores=nb_cores)
  bestid = order(results$residuals)[1]
  parameters = results$params[bestid,]
  bestfit = results$residuals[bestid]
  
  infos = generate_infos(csv_files, inits, sort(results$residuals)[1:5], method, model_links, model_name)
  self = MRAmodelSet(nb_submodels, model, core0$design, model_structure, basal_activity, data_, cv, parameters, bestfit, infos$name, infos$infos, unused_perturbations=unused_perturbations, unused_readouts=unused_readouts, min_cv=MIN_CV, default_cv=DEFAULT_CV)
  # param_range, lower_values, upper_values defined using profile likelihood
  return(self)
}

#' Extend a MRAmodelSet with relevant variable parameters
#'
#' See if a MRAmodelSet requires some parameters to be variable among models to explain the variations
#'
#' @param original_modelset An MRAmodelSet object
#' @param nb_cores Number of cores to use for the refitting with the new variable parameters
#' @param max_iterations Maximum number of variable parameters to add (if 0 takes as many as possible)
#' @param nb_samples Number of samples to generate to fit the new variable parameters
#' @param accuracy Cutoff probability for the chi^2 test
#' @param method Name of the LHS method to use
#' @return An updated MRAmodelSet with the new parameter sets
#' @export
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
addVariableParameters <- function(original_modelset, nb_cores=0, max_iterations=0, nb_samples=100, accuracy=0.95, method="geneticlhs") {
  # clone MRAmodelSet object to have seperate objects at hand
  modelset = cloneModel(original_modelset)
  
  init_time = proc.time()["elapsed"];
  if (nb_cores == 0) { nb_cores = detectCores()-1 }
  model = modelset$model
  
  # Look which parameters are not already variable
  total_parameters = 1:(model$nr_of_parameters()/modelset$nb_models)
  nb_sub_params = length(total_parameters)
  max_iterations = ifelse(max_iterations<=0,nb_sub_params,max_iterations)
  
  for (it in 1:max_iterations) {
    if (length(modelset$variable_parameters) > 0) {
      extra_parameters = total_parameters[-modelset$variable_parameters]
    } else{
      extra_parameters = total_parameters
    }
    # find the parameter which fitted separately to each model improves the performance most and if significant keep variable
    psets=sapply(extra_parameters,refitWithVariableParameter,modelset,nb_sub_params,nb_cores,nb_samples)
    
    bestres = min(unlist(psets["residuals",]))
    deltares = modelset$bestfit - bestres
    df1 = modelset$nb_models-1 # Extra parameters
    data_count = sum(!is.na(modelset$data$stim_data) & !is.nan(modelset$data$stim_data))
    df2 = data_count - (length(modelset$parameters) + df1) # Degrees of freedom of the augmented model
    f_score = (deltares/df1) / (bestres/df2)
    if (f_score > qf(accuracy, df1, df2)) {
      res_id = which.min(unlist(psets["residuals",]))
      par_id = psets["added_var",ceiling(res_id/nb_samples)][[1]]
      new_parameters=unlist(psets["params",ceiling(res_id/nb_samples)][[1]][ifelse(res_id %% nb_samples==0,nb_samples,res_id %% nb_samples),])
      
      message(paste0("variable parameter found: ",model$getParametersLinks()[par_id], "; p-value: ", signif(1-pf(f_score, df1, df2),2) ))
      message(paste0("fitting improvement: ", round(modelset$bestfit,2), "(old) - ", round(bestres,2), "(new) = ", round(deltares,2))) 
      message(paste0("old parameter:", signif(modelset$parameters[par_id],4), " new parameters: ", paste0(signif(new_parameters[seq(from=par_id, to=model$nr_of_parameters(), by=nb_sub_params)],4),collapse=" " ) ))
       
      var_pars = c( modelset$variable_parameters, par_id )
      modelset = setVariableParameters(modelset, var_pars)
      modelset$parameters = new_parameters
      modelset$bestfit = bestres
      get_running_time(init_time, "elapsed");
    } else {
      break
    }
  }
  get_running_time(init_time, "in total");
  return(modelset)
}

#' Refit modelset by iterative relaxation of single parameters
#'
#' @param var_par ID of the variable parameter to test
#' @param modelset A MRAmodelSet object
#' @param nb_sub_params Number of parameters per submodel
#' @param nb_cores Number of cores to use for the fitting
#' @param nb_samples Number of samples to generate for the fitting
#' @param method Method to use for the fitting
#' @return A list with the fields 'residuals' (fitted residuals), added_var (=var_par) and 'params' (the fitted parameter sets corresponding to the residuals)
refitWithVariableParameter <- function(var_par, modelset, nb_sub_params, nb_cores=0, nb_samples=5, method="geneticlhs"){
  model=modelset$model
  var_pars = unique(c( var_par, modelset$variable_parameters ))
  model$setVariableParameters(var_pars)
  new_pset = matrix( rep(modelset$parameters, nb_samples), ncol=length(modelset$parameters), byrow=T )
  
  samples = getSamples(modelset$nb_models, nb_samples, method, nb_cores)
  
  new_pset[,seq(from=var_par, to=model$nr_of_parameters(), by=nb_sub_params)] = samples
  
  refit = parallel_initialisation(model, modelset$data, new_pset, nb_cores)
  
  return(list(residuals = refit$residuals, added_var = var_par, params = refit$params))
}

#' Perform an initialisation of the model 
#'
#' Possibility to use different sampling methods
#' @param model A MRAmodel object to initialise
#' @param core A list containing the design and data
#' @param inits Number of random samples to use for the initialisation
#' @param precorrelate Whether the parameters should be precorrelated
#' @param method The LHS method to use to generate the random samples
#' @param nb_cores Maximum number of cores to use for the initialisation
#' @param perform_plots Whether the distribution of the residuals and the correlation plots for the parameters deduced by correlation should be plotted or not
#' @return A list with the initialisation results
#' @family Model initialisation
initModel <- function(model, core, inits, precorrelate=T, method="randomlhs", nb_cores=1, perform_plots=F) {
  if (inits <= 0) {
    stop("Number of initialisations must be a positive integer")
  }
  expdes = core$design
  data = core$data
  # Parallelized version uses all cores but one to keep control
  if (nb_cores == 0) { nb_cores = detectCores()-1 }
  message (paste("Initializing the model parameters... (", inits, " random samplings) with ", nb_cores, " cores", sep=""))
  # Correlate directly measured and connected nodes -> they will not be sampled
  nr_known_par=0
  if (precorrelate){
    correlated = correlate_parameters(model, core, perform_plot=perform_plots)
    nr_known_par=length(correlated$list)
  }
  samples = getSamples(model$nr_of_parameters()-nr_known_par, inits, method, nb_cores)
  
  # assign sampled and precorrelated samples to the respective parameter positions
  if(precorrelate){
    for (ii in correlated$list){
      inset=rep(correlated$values[correlated$list==ii], times=nrow(samples))
      if (ii==1){
        samples=cbind(inset,samples)
      } else {
        samples=cbind(samples[,1:(ii-1)],inset,samples[,(ii:ncol(samples))])
      }
    }
  }
  message("Sampling terminated.")
  
  #  fit all samples to the model
  results = parallel_initialisation(model, data, samples, nb_cores)
  message("Fitting completed")
  return(results)
}

#' Get LHS samples
#'
#' @param sample_size Size of each sample
#' @param nb_samples Number of samples
#' @param method Name of the LHS method to use
#' @param nb_cores Number of cores to use to generate the samples
#' @return A matrix with the random samples
#' @author Bertram Klinger \email{klinger@@charite.de}
getSamples <- function(sample_size, nb_samples, method="randomlhs", nb_cores=1) {
  if (nb_samples > 10^8){
    warning("Number of samples is too high, restricting sample size to 10^8 !")
    nb_samples = 10^8
  }
  valid_methods=data.frame(methods=c("randomlhs",
                                     "geneticlhs",
                                     "improvedlhs",
                                     "maximinlhs",
                                     "optimumlhs"),
                           calls=c("randomLHS(sample_stack_size, sample_size)",
                                   "geneticLHS(sample_stack_size, sample_size, pop=50)",
                                   "improvedLHS(sample_stack_size, sample_size, dup=5)",
                                   "maximinLHS(sample_stack_size, sample_size, dup=5)",
                                   "optimumLHS(sample_stack_size, sample_size, maxSweeps = 2, eps = 0.1)"
                           ),stringsAsFactors = F)
  idx = grep(method,valid_methods$methods,ignore.case = T)
  if (length(idx==1)){
    max_sample_stack_size=ifelse(valid_methods$methods[idx]=="optimumlhs",100,1000)
    sample_stack_size = min(max_sample_stack_size, nb_samples)
    max_proc=ceiling(nb_samples/max_sample_stack_size)
    if (nb_cores>1){
      samples=qnorm(do.call(rbind,mclapply(1:max_proc,function(x) eval(parse(text=valid_methods$calls[idx])),mc.cores = min(nb_cores,max_proc))), sd=2)
    } else {
      samples=qnorm(do.call(rbind,lapply(1:max_proc,function(x) eval(parse(text=valid_methods$calls[idx])))), sd=2)  
    }
  } else {
    stop(paste0("The selected initialisation method '", method, "' does not exist, valid methods are : ",paste(valid_methods$methods,collapse=", ")))
  }
}

# Gives a parameter vector with a value for correlated parameters and 0 otherwise
parametersFromCorrelation <- function(model, core, perform_plot=F) {
  correlated = correlate_parameters(model, core, perform_plot)
  param_vector = c()
  j = 1
  for (i in 1:model$nr_of_parameters()) {
    if (i %in% correlated$list) {
      param_vector = c(param_vector, correlated$values[j])
      j = j + 1
    } else {
      param_vector = c(param_vector, 0)
    }
  }
  if (class(model) == "Rcpp_ModelSet") {
    message("Dealing with ModelSet")
    param_vector = rep(param_vector, model$nb_submodels)
  }
  return(param_vector)
}

# Indentify the links that can be deduced by a simple correlation, calculate them, and return their index and the corresponding parameter vector with the other values set to 0
correlate_parameters <- function(model, core, perform_plot=F) {
  message("Looking for identifiable correlations...")
  expdes = core$design
  model_structure = core$structure
  data = core$data
  # Collect the nodes that can be identified by correlation
  measured_nodes = expdes$measured_nodes+1;
  valid_nodes = c()
  upstreams = list()
  for (node in measured_nodes) {
    valid = TRUE
    senders = c()
    for (sender in 1:length(model_structure$names)) {
      if (model_structure$adjacencyMatrix[node, sender] == 1) {
        # All sender nodes must be measured for the input links to be indentifiable
        if (sender %in% measured_nodes) {
          senders = c(senders, sender)
        } else {
          valid = FALSE
          if (verbose) { message(paste(model_structure$names[node], "excluded because", model_structure$names[sender], "is not measured")) }
          break
        }
      }
    }
    if (valid && length(senders) > 0) {
      upstreams[[node]] = senders
      valid_nodes = c(valid_nodes, node)
    }
  }
  if (verbose){
  message(paste0(length(valid_nodes), " target node",ifelse(length(valid_nodes)>1,"s","") ,"found correlatable with inputs")) 
  }
  # Collect the values and perform the correlation
  params_matrix = matrix(0, ncol=length(model_structure$names), nrow=length(model_structure$names)) # To store the values
  for (node in valid_nodes) {
    use = rep(TRUE, times=nrow(expdes$inhibitor))
    # We do not use the conditions where one of the sender nodes is inhibited (since it decorelates the measurements)
    for ( sender in which((expdes$inhib_nodes+1) %in% upstreams[[node]]) ) {
      for ( i in 1:nrow(expdes$inhibitor) ) {
        if (expdes$inhibitor[i, sender] == 1) {
          use[i] = FALSE
        }
      }
    }
    
    # Collect the measured data for the regression
    node_mes = which( measured_nodes == node)
    measurements = log(data$stim_data[use, node_mes,drop=F] / mean(data$unstim_data[1, node_mes],na.rm=T) )
    regression = paste0('lm(measurements[,1] ~ ')
    first = TRUE
    node_name = model_structure$names[node]
    condition = paste(node_name, "and")
    for (sender in upstreams[[node]] ) {
      mes_index = which(measured_nodes==sender)
      measurements = cbind(measurements, log(data$stim_data[use, mes_index,drop=F] / mean(data$unstim_data[1,mes_index],na.rm=T)))
      if (first) {
        first = FALSE
      } else {
        regression = paste0(regression, "+")
        condition = paste(condition, "+")
      }
      regression = paste0(regression, "measurements[,", 1+which(upstreams[[node]] == sender), "]")
      sender_name = model_structure$names[sender]
      condition = paste(condition, sender_name)
    }
    # Perform the regression and put the values in the adjacency matrix
    regression = paste0(regression, ")")
    result = eval(parse(text=regression))
    condition = paste(condition, ".\n R^2 =", signif(summary(result)$r.squared, 2) )
    for (sender in 1:length(upstreams[[node]])) {
      params_matrix[ node, upstreams[[node]][sender] ] = result$coefficients[sender+1]
      # TODO add the information on the quality of the fit
    }
    if (perform_plot) {
      # Plot the data and the fit, directly plot the correlation curve if there is only one upstream node, or plot the distance between each measurements and the hyperplane
      if (length(upstreams[[node]]) == 1) {
        plot(measurements[,2], measurements[,1], main=condition, xlab=sender_name, ylab=node_name)
        lines(measurements[,2], result$coefficients[1] + result$coefficients[2] * measurements[,2])
      } else {
        plot(1:nrow(measurements), measurements[,1], pch=4, main=condition, ylab=node_name, xlab="Conditions")
        for (measure in 1:nrow(measurements)) {
          fitted = result$coefficients[1]
          for (sender in 2:ncol(measurements)) {
            fitted = fitted + result$coefficients[sender] * measurements[measure, sender]
          }
          points(measure, fitted, col="blue", pch=20)
          lines(rep(measure, 2), c(fitted, measurements[measure,1]), col="red")
        }
      }
    }
  }
  if (verbose) {
    message(model_structure$names)
    message(params_matrix)
  }
  
  # Each link identified by correlation will correspond to one parameter
  # We identify which parameter it is and return it with its value
  params_vector = model$getParameterFromLocalResponse( params_matrix, rep(0, times=ncol(expdes$inhibitor)) )
  correlated = list()
  correlated$list = c()
  correlated$values = c()
  correlated$infos = c()
  links = model$getParametersLinks()
  for ( param in 1:length(params_vector) ) {
    if (!is.infinite(params_vector[param]) && !is.na(params_vector[param]) && !is.nan(params_vector[param]) && params_vector[param] != 0) {
      correlated$list = c(correlated$list, param)
      correlated$values = c(correlated$values, params_vector[param])
      correlated$infos = c(correlated$infos, paste0("Correlated link ", simplify_path_name(links[param]), " = ", params_vector[param]))
      message(correlated$infos[length(correlated$infos)])
    }
  }
  return(correlated)
}

# Parallel initialisation of the parameters
# @family Model initialisation
parallel_initialisation <- function(model, data, samples, NB_CORES, keep_constant=c()) {
  if (NB_CORES == 0) {
    NB_CORES = detectCores()-1
  }
  # Put the samples under list format, as mclapply only take "one dimension" objects
  nb_samples = nrow(samples)
  parallel_sampling = lapply(1:nb_samples,function(x) samples[x,])
  # Parallel initialisations
  # The number of cores used depends on the ability of the detectCores function to detect them
  fitmodel_wrapper <- function(params, data, model, keep_constant=c()) {
    init = proc.time()[3]
    if (length(keep_constant) > 0) {
      result = model$fitmodelWithConstants(data, params, keep_constant)
    } else {
      result = model$fitmodel(data, params)
    }
    if (verbose) {
      message(paste(signif(proc.time()[3]-init, 3), "s for the descent, residual =", result$residuals))
    }
    return(result)
  }
  
  fitmodelset_wrapper <- function(params, data, model) {
    # TODO add keep_constant control
    if ( class(data) != "Rcpp_DataSet" ) { stop("MRAmodelSet require a STASNet::DataSet object") }
    init = proc.time()[3]
    result = model$fitmodelset(data, params)
    if (verbose) {
      message(paste(signif(proc.time()[3]-init, 3), "s for the descent, residual =", result$residuals))
    }
    return(result)
  }
  
  get_parallel_results <- function(model, data, samplings, NB_CORES, keep_constant=c()) {
    if (class(model) == "Rcpp_ModelSet") {
      parallel_results = mclapply(samplings, fitmodelset_wrapper, data, model, mc.cores=NB_CORES)
    } else {
      parallel_results = mclapply(samplings, fitmodel_wrapper, data, model, mc.cores=NB_CORES, keep_constant)
    }
    
    # Reorder the results to get the same output as the non parallel function
    results = list()
    results$residuals = c()
    results$params = c()
    if (!is.list(parallel_results[[1]])) {
      stop(parallel_results[[1]])
    }
    for (entry in parallel_results) {
      results$residuals = c(results$residuals, entry$residuals)
      results$params = rbind(results$params, entry$parameter)
    }
    return(results)
  }
  # Since the aggregation in the end of the parallel calculations takes a lot of time for a big list, we calculate by block and take the best result of each block
  if (nb_samples > 10000) {
    
    # chop the data into optimal number of pieces
    nb_blocks = nb_samples %/% 10000 + 1
    block_size = NB_CORES * ((nb_samples %/% nb_blocks) %/% NB_CORES)
    best_keep = 10000 %/% nb_blocks
    best_results = list( residuals=c(), params=c() )
    message(paste0("Parallel initialization using block sizes of ", block_size, " samples"))
    
    for (i in 1:nb_blocks) {
      seq = (block_size*(i-1)+1) : ifelse(block_size*(i+1) <= nb_samples, block_size*i, nb_samples)
      results = get_parallel_results(model, data, parallel_sampling[seq], NB_CORES)
      
      # Only keep the best fits
      best = order(results$residuals)[1:best_keep]
      if (verbose) { message(paste(sum(is.na(results$residuals)), "NAs"), stderr()) }
      best_results$residuals = c(best_results$residuals, results$residuals[best])
      best_results$params = rbind(best_results$params, results$params[best,])
    }
    
  } else {
    best_results = get_parallel_results(model, data, parallel_sampling, NB_CORES, keep_constant)
    if (verbose) { message(paste(sum(is.na(best_results$residuals)), "NAs"), stderr()) }
  }
  
  return(best_results)
}

# Initialise the parameters with a one core processing
# needs corrections, not used anyway
classic_initialisation <- function(model, data, nb_samples) {
  residuals = c()
  for (i in 1:nb_samples) {
    result = model$fitmodel( data, samples[i,] )
    residuals = c(residuals,result$residuals)
    params = cbind(params,result$parameter)
    if (i %% (nb_samples/20) == 0) {
      message(paste(i %/% (nb_samples/20), "/ 20 initialization done."))
    }
  }
  if (debug) {
    message(sort(residuals))
  }
  results = list()
  results$residuals = resisuals
  results$params = params
  return(results)
}

#' Extraction of object from file or object
#' @param to_detect An object to check or a string. A string is interpreted as a path to a file containing a matrix compatible with the target object. 
#' @name extraction
NULL

#' Extract a ModelStructure
#'
#' 'extractStructure' detects the format of the structure file and extract the structure of the network. It expects a matrix with 2 columns (adjacency list) or an square matrix (adjacency matrix)
#' @param name Name for the structure. Extracted automatically from the file name if used.
#' @return 'extractStructure' returns a C++ object of class 'ModelStructure'
#' @rdname extraction
#' @export
extractStructure <- function(to_detect, names="") {
  model_links = to_detect
  struct_name = paste0(names, collapse="_")
  if (is.string(model_links)) {
      # Read the file, and extract the matrix corresponding to it
      structure_file = readLines(model_links)
      splitlist = strsplit(structure_file, ",|->|;|\\ |\t")
      if (length(splitlist) < 3) {
        stop("This is a trivial network, you don't need a simulation !")
      } else if ( length(splitlist[[1]]) == 1 ) {
        splitlist = splitlist[-1] # If the number of links is indicated at the beginning, remove it
      }
      struct_matrix = c()
      for (ii in 1:length(splitlist)) {
        struct_matrix = suppressWarnings( rbind(struct_matrix, unlist(splitlist[[ii]])) )
      }
      struct_name = basename(model_links)
      model_links = struct_matrix
  }
  # Control that the format is valid
  if (is.matrix(model_links) || is.data.frame(model_links)) {
      struct_matrix = model_links
  } else if (class(model_links) == "Rcpp_ModelStructure") {
      return(model_links)
  } else {
      stop("Format of the object not compatible with a ModelStructure ")
  }
  
  # Detect if it is a list of links or an adjacency matrix
  if (ncol(struct_matrix) == 2) {
    links_list = struct_matrix
  } else if (ncol(struct_matrix) == 3) {
    links_list = struct_matrix[,1:2] # Values are given, do not use them
  } else {
    # Remove the number of nodes (used for C inputs for example)
    if ( length(struct_matrix[1,]) == 1 && struct_matrix[1, 1] != 1 && struct_matrix[1, 1] != 0) {
      struct_matrix = struct_matrix[,-1]
    }
    # If absent, set the column names: from the first line if it was read from a file, otherwise use the names provided, and in last resort use numbers
    if (!is.null(colnames(struct_matrix))) {
    } else if (suppressWarnings( is.na(as.numeric(struct_matrix[1,1])) )) {
      colnames(struct_matrix) = struct_matrix[1,]
      struct_matrix = struct_matrix[-1,]
      if (suppressWarnings( is.na(as.numeric(struct_matrix[1,1])) )) {
        struct_matrix = struct_matrix[,-1]
      }
    } else if (length(names) == ncol(struct_matrix)) {
      colnames(struct_matrix) = names
    } else if (is.null(colnames(struct_matrix))) {
      warning("No names were provided for the nodes, using numbers instead")
      colnames(struct_matrix) = paste0("node", 1:ncol(struct_matrix))
    }
    rownames(struct_matrix) = colnames(struct_matrix)
    # Check that the dimensions are correct
    if ( ncol(struct_matrix) != nrow(struct_matrix) ) {
      stop("The adjacency matrix is not square. Abort...")
    }
    
    links_list = c();
    for (i in 1:nrow(struct_matrix)) {
      for (j in 1:ncol(struct_matrix)) {
        if (i != j && as.numeric(struct_matrix[i,j]) != 0) {
          links_list = rbind(links_list, c(colnames(struct_matrix)[j], rownames(struct_matrix)[i]))
        }
      }
    }
    
  }

  model_structure=getModelStructure(links_list, struct_name)
  
  return(model_structure)
}

#' Extract the basal activity
#'
#' 'extractBasalActivity' detects if the a list of names or a file perform the extraction if necessary. It expects a vector or an object convertible to a vector.
#' @return 'extractBasalActivity' returns vector of character
#' @rdname extraction
#' @export
# TODO extraction of the basal activity with different format
extractBasalActivity <- function(to_detect) {
  if (is.string(to_detect) && to_detect != "") {
    return(as.character(read.delim(to_detect,header=FALSE)[,1]))
    #unlist(read.delim(to_detect,header = F,colClasses = "character"))
  } else {
    return( as.character(as.vector(to_detect)) )
  }
}

#' Extract a MIDAS dataset
#'
#' 'extractMIDAS' detects if the argument is a MIDAS dataset or the name of a file containing such dataset. It expects a matrix in MIDAS format (See Details).
#' @details A MIDAS format matrix is a matrix with specific column names:
#' \itemize{
#'  \item { One column 'ID:type' containing the type of experiment of each row:\itemize{
#'    \item {'c' or 'control' for control experiment with no treatment}
#'    \item {'blank' for blank experiment with no treatment}
#'    \item {'t' for experiment with a treatment}
#'  }}
#'  \item { At least one column 'TR:STIMULATION' or 'TR:INHIBITIONi' containing 0 or 1 whether the perturbation is applied for the experiment in each row }
#'  \item { One column 'DA:all' or the relevant number of columns 'DA:PERTURBATION' containing the time of the perturbation (Not used in STASNet) }
#'  \item { At least one column with column name 'DV:READOUT' containing the measured values for each readout }
#' }
#' @return 'extractMIDAS' returns a matrix or a data.frame under MIDAS format
#' @rdname extraction
#' @export
#' @seealso \code{link{readMIDAS}}
extractMIDAS <- function(to_detect) {
  if (is.string(to_detect)) {
    if (nchar(to_detect)==0) {
      stop("Can't extractMIDAS from an empty string")
    }
    detected = read.delim(to_detect, sep=",")
    if (ncol(detected)<=1) { # Accept tab-delimited MIDAS files
      detected = read.delim(to_detect, sep="\t")
    }
    checkMIDAS(detected)
    return(detected)
  } else if (is.matrix(to_detect) || is.data.frame(to_detect)) {
    checkMIDAS(to_detect)
    return(to_detect)
  }
  stop("Format of the object not compatible with a MIDAS data")
}

#' Plot a graph from an adjacency list
#'
#' @param structure A 2-columns matrix or a ModelStructure object. The network as an adjacency list, the first column is the upstream nodes, the second column the downstream nodes. Or a ModelStructure object as returned by getModelStructure.
#' @param expdes An ExperimentalDesign object. The measured, stimulated and inhibited nodes are highlighted if present. Signalling strengths are indicated in leftshifted edges and inhibitor strengths are denoted in red below the inhibited node.
#' @param local_values A list with entries 'local_response' (A weighted adjacency matrix representing the values of the links) and 'inhibitors' (A list of inhibition values) both compatible with the 'structure' input
# @export
#' @family Network graph
plotNetworkGraph <- function(structure, expdes="", local_values="") {
    if (class(structure) == "matrix") {
        names = unique(as.vector(structure))
        adm=matrix(0,length(names),length(names),dimnames = list(names,names))
        for (ii in 1:nrow(structure)) {
            adm[match(structure[ii,2],rownames(adm)), structure[ii,1]] = 1
        }
    } else if (class(structure) == "Rcpp_ModelStructure") {
        adm = structure$adjacencyMatrix
        colnames(adm) = rownames(adm) = structure$names
    } else {
        stop("Invalid 'structure' in plotNetworkGraph, must be an edge list or a ModelStructure")
    }
  
  len=length(rownames(adm))
  g1 <- graphAM(adjMat=t(adm),edgemode="directed")
  
# add inhibitors as pseudo nodes downstream of inhibited nodes in order to depict their strength  
  if (class(expdes) == "Rcpp_ExperimentalDesign" && any(local_values != "")){
    if (length(expdes$inhib_nodes)>0){
      for (nn in rownames(adm)[1+expdes$inhib_nodes]){
        g1 <- addNode(paste0(nn,"i"),g1)
        g1 <- addEdge(nn,paste0(nn,"i"),g1)
      }
    }
  }
  
  # setting of general and creation of changed properties
  nodeRenderInfo(g1) <- list(shape="ellipse")
  nodeRenderInfo(g1) <- list(textCol="black")
  nodeRenderInfo(g1) <- list(lwd=1)
  edgeRenderInfo(g1) <- list(fontsize=10)
  edgeRenderInfo(g1) <- list(textCol="black")
  edgeRenderInfo(g1) <- list(col="black")
  
  g1 <- layoutGraph(g1)

  # Add the experimental setup if provided
  if (class(expdes) == "Rcpp_ExperimentalDesign") {
    nodeRenderInfo(g1)$fill[1+expdes$measured_nodes] = "#ffff66"
    
    if (length(expdes$inhib_nodes)>0){
      nodeRenderInfo(g1)$lwd[1+expdes$inhib_nodes]=4 # Populate for perturbations
      nodeRenderInfo(g1)$col[1+expdes$inhib_nodes] = "red"
      nodeRenderInfo(g1)$col[(len+1):(len+length(expdes$inhib_nodes))]="white" # mask inhibitor pseudo nodes
      nodeRenderInfo(g1)$textCol[(len+1):(len+length(expdes$inhib_nodes))]="white" 
    }
    
    if (length(expdes$stim_nodes)>0){
      nodeRenderInfo(g1)$lwd[1+expdes$stim_nodes] = 4
      nodeRenderInfo(g1)$col[1+expdes$stim_nodes] = "blue"
    }
  }
  if (any(local_values != "")) {
    # Add Edge Weights left justified
    efrom = edgeRenderInfo(g1)$enamesFrom
    eto = edgeRenderInfo(g1)$enamesTo
    edge_spline=edgeRenderInfo(g1)$splines
    
    for (idx in which(adm!=0)) {
      vv = local_values$local_response[idx]
      afrom = colnames(adm)[ceiling(idx/len)]
      ato = rownames(adm)[ifelse(idx %% len==0,len,idx %% len)]
      cc = which(afrom==efrom & ato==eto)
      cc = ifelse(length(cc)!=0, cc, which(afrom==eto & ato==efrom)) # Link in both directions
      edgeRenderInfo(g1)$lwd[cc] = ifelse(abs(vv)<=1,1,ifelse(abs(vv)<=5,2,3))
      edgeRenderInfo(g1)$label[cc] = ifelse(vv>=10, round(vv,0),signif(vv, 2))
      
      coordMat=bezierPoints(edge_spline[[cc]][[1]]) # 11 x 2 matrix with x and y coordinates
      edgeRenderInfo(g1)$labelX[cc] = coordMat[5,"x"]-ceiling(nchar(edgeRenderInfo(g1)$label[cc])*10/2)
      edgeRenderInfo(g1)$labelY[cc] = coordMat[5,"y"]
      if (vv < 0) { edgeRenderInfo(g1)$col[cc] = "orange" }
    }
    
    # Add Inhibitor estimates
    if (length(expdes$inhib_nodes)>0){
      for (idx in 1:length(expdes$inhib_nodes)) {
        vv = local_values$inhibitors[idx]
        iname = paste0(colnames(adm)[expdes$inhib_nodes[idx]+1], "i")
        nname = colnames(adm)[expdes$inhib_nodes[idx]+1]
        cc = which(nname==efrom & iname==eto)
        edgeRenderInfo(g1)$col[cc]="white" # mask inhibitor pseudo edges
        edgeRenderInfo(g1)$label[cc] = ifelse(vv>=10, round(vv,0),signif(vv, 2))
        edgeRenderInfo(g1)$textCol[cc]="red"
        
        coordMat=bezierPoints(edge_spline[[cc]][[1]]) # 11 x 2 matrix with x and y coordinates 
        edgeRenderInfo(g1)$labelX[cc] = coordMat[2,"x"]
        edgeRenderInfo(g1)$labelY[cc] = coordMat[2,"y"]
      }
    }
  }
  renderGraph(g1)
  invisible(g1)
  # TODO MARK REMOVED LINKS
}

#' Extracts the data, the experimental design and the structure from the input files
#' @param model_structure Matrix of links [node1, node2]
#' @param basal_activity The node of the structure with a basal activity
#' @param datas Experimental data under the MIDAS format. See extractMIDAS.
#' @param var_file Variation file name under MIDAS format or "" to compute an error from the data. See extractMIDAS.
#' @param dont_perturb Perturbations to be removed for the fit, will not be used nor simulated. (vector of names)
#' @param dont_read Readouts to be removed for the fit, will not be used nor simulated. (vector of names)
#' @param MIN_CV Minimum coefficient of variation.
#' @param DEFAULT_CV Default coefficient of variation to use when none is provided and there are no replicates in the data.
#' @param rearrange Whether the rows should be rearranged. "no" to keep the order of the perturbations from the data file, "bystim" to group by stimulations, "byinhib" to group by inhibitions.
#' @seealso \code{\link{extractMIDAS}}
extractModelCore <- function(model_structure, basal_activity, data_filename, var_filename="", dont_perturb=c(), dont_read=c(), MIN_CV=0.1, DEFAULT_CV=0.3, rearrange="no") {

  model_structure = extractStructure(model_structure)
  basal_activity = extractBasalActivity(basal_activity)
  vpert = c(model_structure$names, paste0(model_structure$names, "i")) # Perturbations that can be simulated
  data_file = extractMIDAS(data_filename)
  data_values = data_file[,grepl("^DV.", colnames(data_file)),drop=FALSE]
  colnames(data_values) = gsub("^[A-Z]{2}.", "", colnames(data_values))
  not_included = setdiff(colnames(data_values), model_structure$names)
  perturbations = data_file[,grepl("^TR.", colnames(data_file)),drop=FALSE]
  if (ncol(perturbations) == 0) {
    stop("Perturbation informations are required")
  }
  colnames(perturbations) = gsub("^[A-Z]{2}.", "", colnames(perturbations))
  not_perturbable = setdiff(colnames(perturbations), vpert)
  id_colums = grepl("^ID.type", colnames(data_file))
  blanks = which(tolower(data_file[,id_colums]) %in% c("b","blank"))
  controls = which(tolower(data_file[,id_colums]) %in% c("c","ctl","ctrl","control"))
  if (length(controls) == 0) { 
    stop("Control experiments are required (indicated as 'control' in the column 'ID:type', see extractMIDAS)")
    }
  # Warn for the measured nodes that have not been found in the network, and don't use them
  if (length(not_included) > 0) {
    message(paste(not_included , "measurement is not in the network structure (could be a mispelling or a case error)\n"))
    data_values = data_values[,-which(colnames(data_values) %in% not_included), drop=FALSE]
  }
  # Remove extra readouts that should not be used
  for (ro in dont_read) {
    if (ro %in% colnames(data_values)) {
      message(paste(ro, "readout will not be used for the fit\n"))
      data_values = data_values[,-which(colnames(data_values)==ro), drop=FALSE]
    }
  }
  # Remove the perturbations that cannot be simulated to get a correct fit
  for (erm in dont_perturb) {
    if (!erm %in% vpert) {
      message(paste0(erm, " is not a valid perturbation name for this dataset"))
    }
  }
  dont_perturb = intersect(dont_perturb, vpert)
  rm_rows = c()
  if (length(not_perturbable) > 0 || length(dont_perturb) > 0) {
    if (length(not_perturbable) > 0) {
      message(paste(not_perturbable , "perturbation not compatible with the network structure, it will not be used\n" ))
    }
    if (length(dont_perturb) > 0) {
      message(paste(dont_perturb, "perturbation will not be used for the fit\n"))
    }
    not_perturbable = c(not_perturbable, dont_perturb)
    rm_rows = unique(unlist(lapply(not_perturbable, function(pp) { which(perturbations[,pp]==1) })))
    perturbations = perturbations[,-which(colnames(perturbations) %in% not_perturbable), drop=F]
    if (length(perturbations)==0){
      stop("All perturbations have been removed can not continue modelling")
    }
    if (!any(perturbations[-rm_rows,]==1)){
      stop("Remaining perturbations only occur in combination with removed perturbations, please reconsider perturbation scheme")
    }
  }
  # Means of the blank fixation of the antibodies
  if (length(blanks) == 0) {
    blank_values = matrix( rep(0, ncol(data_values)), nrow=1, dimnames=list(NULL, colnames(data_values)) )
  } else {
    blank_values = colMeans(data_values[blanks,,drop=F], na.rm=T)
  }
  blank_values[is.nan(blank_values)|is.na(blank_values)] = 0 # For perturbations without blank values 
  # Means of basal activity of antibodies
  unstim_values = colMeans(data_values[controls,,drop=F], na.rm=T)
  if (any(is.nan(unstim_values)|is.na(unstim_values))) {
    stop("Unstimulated data are required to simulate the network")
  }
  # Calculate statistics for variation data
  if (length(c(rm_rows,blanks))>0){
    mean_stat = aggregate(data_values[-c(rm_rows,blanks),,drop=FALSE], by=perturbations[-c(rm_rows,blanks),,drop=FALSE], mean, na.rm=T)[,-(1:ncol(perturbations)),drop=FALSE]
    sd_stat = aggregate(data_values[-c(rm_rows,blanks),,drop=FALSE], by=perturbations[-c(rm_rows,blanks),,drop=FALSE], sd, na.rm=T)[,-(1:ncol(perturbations)),drop=FALSE]
  } else {
    mean_stat = aggregate(data_values, by=perturbations, mean, na.rm=TRUE)[,-(1:ncol(perturbations)),drop=FALSE]
    sd_stat = aggregate(data_values, by=perturbations, sd, na.rm=TRUE)[,-(1:ncol(perturbations)),drop=FALSE]
  }
  # Delete the perturbations that cannot be used, blank and controls from the dataset
  rm_rows = c(rm_rows, controls, blanks) 
  data_values = data_values[-rm_rows,,drop=F]
  perturbations = perturbations[-rm_rows,,drop=F]
  # Compute the mean and standard deviation of the data
  mean_values = aggregate(data_values, by=perturbations, mean, na.rm=T)[,-(1:ncol(perturbations)),drop=F]
  blank_values = matrix( rep(blank_values, each=nrow(mean_values)), nrow=nrow(mean_values), dimnames=list(NULL, colnames(data_values)) )
  unstim_values = matrix(rep(unstim_values, each=nrow(mean_values)), nrow=nrow(mean_values))
  colnames(unstim_values) <- colnames(mean_values)
  if (verbose > 7) {
    message("Data used:")
    message(mean_values)
  } else if (verbose > 5) {
    message("Data used:")
    message(head(mean_values))
  }

  # Error model
  if (any(var_filename != "")) {
    # We use the CV file if there is one
    # The format and the order of the conditions are assumed to be the same as the data file
    message(paste0("Using var file ", var_filename))
    variation_file = extractMIDAS(var_filename)
    # Check that the number of samples is the same for the measurements and the variation and that the names in the measurements file and in the variation file match
    if (nrow(variation_file) != nrow(data_file)) {
      stop("Different number of experiments for the variation and the measurement files") 
    }
    notInVar = setdiff(colnames(data_file),colnames(variation_file))
    if (length(notInVar)>0){
      stop(paste0("Names of the variation and measurement files do not match","\n",
                  "Columns ",paste0(notInVar,collapse=" "),
                  " from data file are not found in var file!"))
    }
    # Reorder and filter variation columns to columns present in data file
    variation_file=variation_file[,colnames(data_file)]
    # Check whether the order and type of the perturbations is preserved
    if (!all(data_file[,grepl("^TR.", colnames(data_file))]==variation_file[,grepl("^TR.", colnames(variation_file))])){
      stop("Order or type of experiments in the variation file is different from the measurement file, please adapt them to be the same!") 
    }
    # Gather the cv values corresponding to the experimental design
    pre_cv = variation_file[, grepl("^DV", colnames(variation_file))]
    colnames(pre_cv) = gsub("^[A-Z]{2}.", "", colnames(pre_cv))
    cv_values = pre_cv[,colnames(pre_cv)%in%colnames(mean_values)]
    cv_values = cv_values[-rm_rows,,drop=FALSE]
    cv_values = aggregate(cv_values, by=perturbations, mean, na.rm=T)[,-(1:ncol(perturbations))]
  } else {
    # remove measurements not different from blank
    mean_stat[mean_stat <= 2 * blank_values] = NA
    median_cv = apply(sd_stat / mean_stat, 2, median, na.rm=T)
    cv_values = matrix(rep(median_cv,each=nrow(mean_values)),nrow=nrow(mean_values))
  }
  
  colnames(cv_values)=colnames(mean_values)
  cv_values[is.nan(as.matrix(cv_values)) | is.na(cv_values)] = DEFAULT_CV
  cv_values[cv_values < MIN_CV] = MIN_CV
  error = blank_values + cv_values * mean_values
  error[error<1] = 1 # The error cannot be 0 as it is used for the fit. If we get 0 (which means blank=0 and stim_data=0), we set it to 1 (which mean the score will simply be (fit-data)^2 for those measurements). We also ensure that is is not too small (which would lead to a disproportionate fit attempt

  # Extract experimental design
  perturbations = aggregate(perturbations, by=perturbations, max, na.rm=T)[,-(1:ncol(perturbations)),drop=F]
  names = colnames(perturbations)
  stim_names = names[grepl("[^i]$", names)]
  stim_nodes = as.character( stim_names[match(model_structure$names, stim_names, nomatch = 0)] )
  names = gsub("i$", "", names[grepl("i$", names)])
  inhib_nodes = as.character( names[match(model_structure$names, names, nomatch = 0 )] )
  
  no_stim=F
  no_inh=F
  if (length(stim_nodes)>0){
    stimuli = as.matrix(perturbations[,stim_nodes,drop=F])
  } else{
    stimuli = perturbations[,0]
    no_stim=T
  }
  
  if (length(inhib_nodes)>0){  
    inhibitor = as.matrix(perturbations[,paste(inhib_nodes, "i", sep=""),drop=F])
  }else if (no_stim){
    stop("Neither a valid stimulation nor inhibition was given, at least one type is needed!")
  }else{
    inhibitor = perturbations[,0]
    no_inh=T
  }
  
  measured_nodes = colnames(mean_values)
  measured_nodes = measured_nodes[measured_nodes %in% model_structure$names] # Preserve the order in the file, allowing to specify the order of the readouts in the input
  
  if (rearrange %in% c("bystim", "byinhib")) {
    # Arrange data according to experimental design
    stim_motifs = aggregate(stimuli,by=as.data.frame(stimuli),max,na.rm=T)[,-(1:ncol(stimuli)),drop=F] # sorts automatically
    inh_motifs = aggregate(inhibitor,by=as.data.frame(inhibitor),max,na.rm=T)[,-(1:ncol(inhibitor)),drop=F] # sorts automatically
    
    error_sort = error[0, measured_nodes,drop=F]
    stim_data_sort = mean_values[0, measured_nodes,drop=F]
    cv_sort = cv_values[0,measured_nodes,drop=F]
    stim_sort = if (no_stim) {as.matrix(stimuli)} else{stimuli[0, ,drop=F]}
    inh_sort = if (no_inh) {as.matrix(inhibitor)} else{inhibitor[0, ,drop=F]}
    if (rearrange == "byinhib") { # Invert variables
      tmp = stim_sort
      stim_sort = inh_sort
      inh_sort = tmp
      tmp = no_stim
      no_stim = no_inh
      no_inh = tmp
      tmp = stim_motifs
      stim_motifs = inh_motifs
      inh_motifs = tmp
      tmp = stimuli
      stimuli = inhibitor
      inhibitor = tmp
    }
    
    if (!no_stim){
      for (ii in 1:nrow(stim_motifs)){
        stim_pos = apply(stimuli==matrix(rep(stim_motifs[ii,], each=nrow(stimuli)), nrow=nrow(stimuli)), 1, all)
        if (!no_inh){
          for (jj in 1:nrow(inh_motifs)){
            inh_pos = apply(inhibitor==matrix(rep(inh_motifs[jj,], each=nrow(inhibitor)), nrow=nrow(inhibitor)), 1, all)
            if (any(stim_pos & inh_pos)){
              error_sort = rbind(error_sort, error[stim_pos & inh_pos, measured_nodes,drop=F])
              stim_data_sort = rbind(stim_data_sort, mean_values[stim_pos & inh_pos, measured_nodes, drop=F])
              cv_sort = rbind(cv_sort, cv_values[stim_pos & inh_pos, measured_nodes, drop=F])
              stim_sort = rbind(stim_sort, stimuli[stim_pos & inh_pos, ,drop=F])
              inh_sort = rbind(inh_sort, inhibitor[stim_pos & inh_pos, ,drop=F])
            }
          }
        }else if (any(stim_pos)){
          error_sort = rbind(error_sort, error[stim_pos, measured_nodes,drop=F])
          stim_data_sort = rbind(stim_data_sort, mean_values[stim_pos, measured_nodes,drop=F])
          cv_sort = rbind(cv_sort, cv_values[stim_pos, measured_nodes,drop=F])
          stim_sort = rbind(stim_sort, stimuli[stim_pos, ,drop=F])
        }
      }
    }else{
      for (jj in 1:nrow(inh_motifs)){
        inh_pos = apply(inhibitor==matrix(rep(inh_motifs[jj,], each=nrow(inhibitor)), nrow=nrow(inhibitor)), 1, all)
        if (any(inh_pos)){
          error_sort = rbind(error_sort, error[inh_pos, measured_nodes,drop=F])
          stim_data_sort = rbind(stim_data_sort, mean_values[inh_pos, measured_nodes,drop=F])
          cv_sort = rbind(cv_sort, cv_values[inh_pos, measured_nodes,drop=F])
          inh_sort = rbind(inh_sort, inhibitor[inh_pos, ,drop=F])
        }
      }
    }

    if (rearrange == "byinhib") { # Re-invert the variables that are used later
      tmp = stim_sort
      stim_sort = inh_sort
      inh_sort = tmp
    }
  } else {
    if (!rearrange %in% c("n", "no", "")) {
        warning(paste0("Unknown option '", rearrange, "' for the arrangement of the perturbations, interpreted as 'no'"))
    }
    stim_sort = as.matrix(stimuli)
    inh_sort = as.matrix(inhibitor)
    stim_data_sort = mean_values[,measured_nodes, drop=FALSE]
    cv_sort = cv_values[,measured_nodes, drop=FALSE]
    error_sort = error[,measured_nodes, drop=FALSE]
  }
  
  if (verbose > 3) {
    message("Stimulated nodes")
    message(stim_sort)
    message("Inhibited nodes")
    message(inh_sort)
  }
  
  expdes=getExperimentalDesign(model_structure, stim_nodes, inhib_nodes, measured_nodes, stim_sort, inh_sort, basal_activity)
  
  # Create data object
  data=new(STASNet:::Data)
  data$set_unstim_data ( unstim_values[, measured_nodes, drop=FALSE] )
  data$set_scale( blank_values[, measured_nodes, drop=FALSE] )
  data$set_stim_data( as.matrix(stim_data_sort[, measured_nodes, drop=FALSE]) )
  data$set_error( as.matrix( error_sort[, measured_nodes, drop=FALSE] ))
  
  # Finalize object 
  core = list()
  core$design = expdes
  core$data = data
  core$structure = model_structure
  core$basal = expdes$basal_activity
  core$cv = as.matrix(cv_sort)

  return(core)
}

#' Import a save model with data files
#'
#' Build a fitted model from a .mra file, and import data for this model
#' Does NOT perform any initialisation
#' @param model_file A .mra file containing the information on the model
#' @param data_file A .csv file with the data for the model
#' @param var_file A .var file with the variation of the data
#' @return An MRAmodel object describing the model and its best fit, containing the data
#' @export
#' @seealso importModel, exportModel, createModel
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
#' @examples \dontrun{
#' rebuildModel("model.mra", "data.csv", "data.var")
#' }
rebuildModel <- function(model_file, data_file, var_file="", rearrange="no") {
  if (!grepl(".mra$", model_file)) {
    stop("The model file does not have the mra extension")
  }
  model = importModel(model_file)
  #links = matrix(rep(model$structure$names, 2), ncol=2)
  core = extractModelCore(model$structure, model$basal, data_file, var_file, model$unused_perturbations, model$unused_readouts, model$min_cv, model$default_cv, rearrange=rearrange)
  
  model = MRAmodel(model$model, model$design, model$structure, model$basal, core$data, core$cv, model$parameters, model$model$fitmodel(core$data, model$parameters)$residuals, model$name, model$infos, model$param_range, model$lower_values, model$upper_values, model$unused_perturbations, model$unused_readouts, model$min_cv, model$default_cv)
  
  return(model)
}

# Perform several simulated annealing, one per core
parallelAnnealing <- function(model, core, max_it, nb_cores, perform_plot=F) {
  correlation = parametersFromCorrelation(model, core, perform_plot=perform_plot)
  correlation_list = list()
  for (i in 1:nb_cores) {
    correlation_list[[i]] = correlation
  }
  # Perform several annealing to select the best
  annealings = mclapply(correlation_list, simulated_annealing_wrapper, model, core$data, max_it, 0, mc.cores=nb_cores)
  results = list()
  results$residuals = c()
  results$params = c()
  for (annealing in annealings) {
    results$residuals = c(results$residuals, annealing$residuals)
    results$params = rbind(results$params, annealing$parameter)
  }
  
  return(results)
}

# Wrapper for the C++ function because Rcpp does not take into account optionnal arguments (and R crashes if there is not exactly the correct number of arguments)
# and for the parallelisation
simulated_annealing_wrapper <- function(correlated_parameters, model, data, max_it=0, max_depth=0) {
  return(model$annealingFit(data, correlated_parameters, max_it, max_depth))
}

# TODO integrate it into the creation function
# Explore the parameter space with series of random initialisations and selection after each step of the best fits (mutation only genetic algorithm)
MIN_SAMPLE_SIZE = 100 # TODO Find a heuristic value depending on the number of parameters to estimate
deep_initialisation <- function (model, core, NB_CORES, depth=3, totalSamples=100, perform_plot=F) { # TODO Chose to give the total number or the per sampling
  
  if (depth < 1) { depth = 1 } # We do at least one initialisation
  #if (nSamples < MIN_SAMPLE_SIZE) { nSamples = MIN_SAMPLE_SIZE }
  nSamples = ceiling(sqrt(totalSamples / depth))
  message(paste("Mutations per point =", nSamples))
  
  # The first iteration does not use any shift
  kept = matrix(0, ncol=model$nr_of_parameters())
  correlation = correlate_parameters(model, core, perform_plot)
  for (i in 1:depth) {
    message(paste("Depth :", i))
    # Create the new samples, with a random initialisation shifted by the previous steps local minimum and recombinations of those previous best steps
    # We reduce the exploration range at each step
    samples = c()
    for (j in 1:nrow(kept)) {
      # Generate half of the new samples by simple mutation
      new_sample = sampleWithCorrelation(model, core, ceiling(nSamples * 0.5), kept[j,], correlation, sd=3/i)
      samples = rbind( samples, new_sample$samples) 
      # Recombination of the individual j with some of the others
      for (k in 1:ceiling(nSamples * 0.5)) {
        new_sample$samples[k,] = (kept[j,] + kept[sample(1:nrow(kept), size=1),]) / 2
      }
      samples = rbind( samples, new_sample$samples) 
    }
    results = parallel_initialisation(model, core$data, samples, NB_CORES)
    if (perform_plot) { hist(log(results$residuals, base=10), breaks="fd") }
    # Keep the number of samples, and select the best parameters sets for the next iteration
    kept = results$params[order(results$residuals)[1:nSamples],]
  }
  
  # Finish by checking that sampling around the optimum with a high variation still finds several times the optimum
  new_sample = sampleWithCorrelation(model, core, nSamples, kept[1,], correlation, sd=3)$samples
  if (perform_plot) { hist(log(parallel_initialisation(model, core$data, new_sample, NB_CORES)$residuals, base=10), breaks="fd") }
  
  return(results) # Return the last set of fit
}
