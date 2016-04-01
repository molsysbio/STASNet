################################ create_model.R #########################################
## Functions that simplify data parsing and model creation and fitting using the package

#' @import Rgraphviz
#' @import pheatmap
#' @import lattice
#' @import parallel
#' @import Rcpp
#' @import lhs
#' @import RhpcBLASctl
#' @useDynLib fitmodel

# Global variable to have more outputs
verbose = FALSE
debug = TRUE

# Generates a random number between 0 and 1
rand <- function(decimals=4) {
  return(sample(1:10^decimals, size=1) / 10^decimals)
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
#' @param method Method to be used for the initialisation, available methods are :
#'      random : Perform a Latin Hypercube Sampling to choose \emph{inits} starting points then perform a gradient descent to find the local optimum for each of those points.
#'      correlation : Deduce some parameters from the correlation between the measurements for the target node and all of its input nodes, then perform random to find the other parameters. Recommended, very efficient for small datasets.
#'      genetic : Genetic algorithm with mutation only. \emph{inits} is the total number of points sampled.
#'      annealing : Simulated annealing and gradient descent on the best result. \emph{inits} is the maximum number of iteration without any change before the algorithm decides it reached the best value. Use not recommended.
#' @return An MRAmodel object describing the model and its best fit, containing the data
#' @export
#' @seealso importModel, exportModel, rebuildModel
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
#' @examples
#' model = createModel("links.tab", "basal.dat", "data_MIDAS.csv") # Produces a model for the network described in links.tab using the data in data_MIDES.csv
#' model = createModel("links.tab", "basal.dat", "data_MIDAS.csv", "variation.var") # Uses the variation from a variation file
#' model = createModel("links.tab", "basal.dat", "data_MIDAS.csv", nb_cores = detectCores()) # Uses all cores available (with the package parallel)
#' model = createModel("links.tab", "basal.dat", "data_MIDAS.csv", inits = 1000000) # Uses more initialisations for a complex network
createModel <- function(model_links, basal_file, data.stimulation, data.variation="", nb_cores=1, inits=1000, perform_plots=F, precorrelate=T, method="geneticlhs") {
  # Creation of the model structure object
  model_structure = extractStructure(model_links)
  
  # TODO extraction of the basal activity with different format
  #basal_activity = as.character(read.delim(basal_file,header=FALSE)[,1])
  basal_activity =unlist(read.delim(basal_file,header = F,colClasses = "character"))
  # Create a valid basal activity list from a list of the nodes that do not have basal activity
  # TODO Find a day to switch
  # IDEA Consider it is basal if more than 2/3 of the nodes and not basal if less than 1/3, ask otherwise
  #   basal = c()
  #   for (node in model_structure$names) {
  #       if (!(node %in% no_basal)) {
  #           basal = c(basal, node)
  #       }
  #   }
  core = extractModelCore(model_structure, basal_activity, data.stimulation, data.variation)
  expdes = core$design
  data = core$data

  # MODEL SETUP
  model = new(fitmodel:::Model)
  model$setModel(expdes, model_structure)
  
  # INITIAL FIT
  results <- initModel(model, core, inits, precorrelate, method, nb_cores, perform_plots)
  
  # Choice of the best fit
  params = results$params
  residuals = results$residuals
  #write("All residuals : ", stderr())
  #write(residuals, stderr())
  print(paste( sum(is.na(residuals)), "NA and ", sum(is.infinite(residuals)), "infinite residuals" ))
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
    print("Best residuals :")
    print(sort(residuals)[1:20])
  }
  if (perform_plots) { hist(log(residuals, base=10), breaks="fd", main="Distribution of the residuals") }
  range_var <- function(vv) { rr=range(vv); return( (rr[2]-rr[1])/max(abs(rr)) ) }
  paths = sapply(model$getParametersLinks(), simplify_path_name)
  best_sets = order_resid[signif(residuals[order_resid], 4) == signif(residuals[order_resid[1]], 4)]
  for ( ii in 1:(ncol(params)-1) ) {
    for (jj in (ii+1):ncol(params)) {
      setii = params[best_sets,ii] 
      setjj = params[best_sets,jj]
      cij = suppressWarnings(cor(setii, setjj, use="na.or.complete")) # NA if one vector has a standard deviation of 0
      if ((!is.na(cij) && cij > 0.999) || (range_var(setii) > 0.05 && range_var(setjj) > 0.05)) {
        plot(setii, setjj, xlab=paths[ii], ylab=paths[jj], main=paste0("Values for the best fits\ncor=", cij), col=residuals[best_sets])
      }
    }
  }
  
  best_id = order(residuals)[1]
  init_params = params[best_id,]
  init_residual = residuals[best_id]
  if (verbose) {
    print("Model simulation with optimal parameters :")
    print(model$simulate(data, init_params)$prediction)
  }
  
  # Information required to run the model (including the model itself)
  infos = c(paste0(inits, " samplings"), paste0( sort("Best residuals : "), paste0(sort(residuals)[1:5], collapse=" ") ), paste0("Method : ", method), paste0("Network : ", model_links))
  model_name = gsub("\\.csv", "", gsub("_MIDAS", "", basename(data.stimulation)))
  model_description = MRAmodel(model, expdes, model_structure, basal_activity, data, core$cv, init_params, init_residual, name=model_name, infos=infos)
  print(paste("Residual score =", model_description$bestfitscore))
  
  return(model_description)
}

#' Build an MRAmodelSet
#'
#' Build and fit an MRAmodelSet, which consists of the simultaneous fitting of several MRA models
#' @export
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
createModelSet <- function(model_links, basal_nodes, csv_files, var_files=c(), nb_cores=1, inits=1000, perform_plots=F, method="geneticlhs") {
  if (length(csv_files) != length(var_files)) {
    if (length(var_files) == 0) {
      var_files = rep("", length(csv_files))
    } else {
      stop("'var_files' must have the same length as 'csv_files' or be of length 0")
    }
  }
  model_structure = extractStructure(model_links)
  basal_activity = as.character(read.delim(basal_nodes,header=FALSE)[,1])
  
  nb_submodels = length(csv_files)
  core0 = extractModelCore(model_structure, basal_activity, csv_files[1], var_files[1])
  stim_data = core0$data$stim_data
  unstim_data = core0$data$unstim_data
  error = core0$data$error
  cv = core0$cv
  if (verbose > 8) {
    print("Data dimensions~:")
    print(dim(core0$data$stim_data))
    print(dim(core0$data$unstim_data))
    print(dim(core0$data$error))
  }
  # Build an extended dataset that contains the data of each model
  data_ = new(fitmodel:::DataSet)
  data_$addData(core0$data, FALSE)
  for (ii in 2:nb_submodels) {
    core = extractModelCore(model_structure, basal_activity, csv_files[ii], var_files[ii])
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
    data_$addData(core$data, FALSE)
    cv = rbind(cv, core$cv)
  }
  data_$set_stim_data(stim_data)
  data_$set_unstim_data(unstim_data)
  data_$set_error(error)
  data_$set_scale(error)
  
  model = new(fitmodel:::ModelSet)
  model$setModel(core0$design, model_structure)
  model$setNbModels(nb_submodels)
  
  print("setting completed")
  results = initModel(model, list(design=core0$design, data=data_, structure=model_structure), inits, perform_plots=perform_plots, method=method, precorrelate=F, nb_cores=nb_cores)
  bestid = order(results$residuals)[1]
  parameters = results$params[bestid,]
  bestfit = results$residuals[bestid]
  
  infos = c(paste0(inits, " samplings"), paste0( sort("Best residuals : "), paste0(sort(results$residuals)[1:5], collapse=" ") ), paste0("Method : ", method), paste0("Network : ", model_links))
  names = gsub("\\.csv", "", gsub("_MIDAS", "", basename(csv_files)))
  self = MRAmodelSet(nb_submodels, model, core0$design, model_structure, basal_activity, data_, cv, parameters, bestfit, names, infos)
  # param_range, lower_values, upper_values)
  return(self)
}

#' See if a MRAmodelSet requires some parameters to be variable among models to explain the variations
#'
#' @param modelset An MRAmodelSet object
#' @param nb_cores Number of cores to use for the refitting with the new variable parameters
#' @param max_iterations Maximum number of variable parameters to add
#' @param accuracy Cutoff probability for the chi^2 test
#' @param method Name of the LHS method to use
#' @return An updated MRAmodelSet with the new parameter sets
#' @export
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
addVariableParameters <- function(modelset, nb_cores=0, max_iterations=1, accuracy=0.95, method="randomlhs") {
  if (nb_cores == 0) { nb_cores = detectCores()-1 }
  model = modelset$model
  
  # Look which parameters are not already variables
  extra_parameters = 1:(model$nr_of_parameters()/modelset$nb_models)
  max_iterations = max(max_iterations, length(extra_parameters))
  if (max_iterations < 0) { stop("Incorrect number of iterations") }
  for (it in 1:max_iterations) {
    nb_sub_params = length(extra_parameters)
    extra_parameters = 1:(model$nr_of_parameters()/modelset$nb_models)
    if (length(modelset$variable_parameters) > 0) {
      extra_parameters = extra_parameters[-c(modelset$variable_parameters)]
    }
    psets = list()
    for (ii in extra_parameters) {
      var_pars = unique(c( ii, modelset$variable_parameters ))
      model$setVariableParameters(var_pars)
      new_pset = matrix( rep(modelset$parameters, 5), ncol=length(modelset$parameters), byrow=T )
      samples = geneticLHS(nrow(new_pset), modelset$nb_models)
      samples = getSamples(modelset$nb_models, nrow(new_pset), method=method,nb_cores)
      for (jj in 1:modelset$nb_models) {
        new_idx = nb_sub_params*(jj-1)+ii
        new_pset[,new_idx] = new_pset[,new_idx] + samples[,jj]
      }
      
      refit = parallel_initialisation(model, modelset$data, new_pset, NB_CORES=nb_cores)
      
      psets$residuals = c(psets$residuals, refit$residuals)
      psets$added_var = c( psets$added_var, rep(ii, length(refit$residuals)) )
      psets$params = rbind(psets$params, refit$params)
    }
    print("so far...")
    ord_res = order(psets$residuals)
    if (modelset$bestfit - psets$residuals[ord_res[1]] > qchisq(accuracy, modelset$nb_models) ) {
      bid = ord_res[1]
      var_pars = c( modelset$variable_parameters, modelset$variable_parameters[bid] )
      print("so good")
      modelset = setVariableParameters(modelset, var_pars)
      modelset$parameters = psets$params[bid,]
      modelset = computeFitScore(modelset)
    } else {
      break
    }
  }
  
  return(modelset)
}

#' Perform an initialisation of the model 
#' Possibility to use different sampling methods
initModel <- function(model, core, inits, precorrelate=T, method="randomlhs", nb_cores=1, perform_plots=F) {
  expdes = core$design
  data = core$data
  # Parallelized version uses all cores but one to keep control
  if (nb_cores == 0) { nb_cores = detectCores()-1 }
  print (paste("Initializing the model parameters… (", inits, " random samplings) with ", nb_cores, " cores", sep=""))
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
  print("Sampling terminated.")
  
  #  fit all samples to the model
  results = parallel_initialisation(model, data, samples, nb_cores)
  print("Fitting completed")
  return(results)
}

#' Get LHS samples
#'
#' @param sample_size Size of each sample
#' @param nb_samples Number of samples
#' @param method Name of the LHS method to use
getSamples <- function(sample_size, nb_samples, method="randomlhs",nb_cores=1) {
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
    }else {
      samples=qnorm(do.call(rbind,lapply(1:max_proc,function(x) eval(parse(text=valid_methods$calls[idx])))), sd=2)  
    }
  }else {
    stop(paste0("The selected initialisation method '", method, "' does not exist, valid methods are : ",paste(valid_methods$methods,collapse=", ")))
  }
}

# TODO put it in deep_initialisation
#' Produces shifted random samples with some parameters inferred from correlations in the measurements
#' @param model An MRAmodel object
#' @param shift A vector giving the shift to be applied to the normal distribution. Must have the same size as the number of parameters
sampleWithCorrelation <- function(model, core, nb_samples, shift=0, correlated="", sd=2, perform_plot=F) {
  if (!is.list(correlated)) {
    correlated = correlate_parameters(model, core, perform_plot)
  }
  random_samples = qnorm(randomLHS(nb_samples, model$nr_of_parameters()-length(correlated$list)), sd=sd)
  samples = c()
  j=1
  for (i in 1:model$nr_of_parameters()) {
    if (i %in% correlated$list) {
      samples = cbind(samples, rep(correlated$values[which(correlated$list==i)], times=nb_samples))
    } else {
      if ( length(shift) == model$nr_of_parameters() ) {
        samples = cbind(samples, random_samples[,j] + shift[i])
      } else {
        samples = cbind(samples, random_samples[,j])
      }
      j = j+1
    }
  }
  
  cor_samples = list()
  cor_samples$samples = samples
  cor_samples$cor = correlated
  print("Sampling terminated.")
  return(cor_samples)
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
    print("Dealing with ModelSet")
    param_vector = rep(param_vector, model$nb_submodels)
  }
  return(param_vector)
}

# Indentify the links that can be deduced by a simple correlation, calculate them, and return their index and the corresponding parameter vector with the other values set to 0
correlate_parameters <- function(model, core, perform_plot=F) {
  print("Looking for identifiable correlations...")
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
          if (verbose) { print(paste(model_structure$names[node], "excluded because", model_structure$names[sender], "is not measured")) }
          break
        }
      }
    }
    if (valid && length(senders) > 0) {
      upstreams[[node]] = senders
      valid_nodes = c(valid_nodes, node)
    }
  }
  print(paste(length(valid_nodes), "target nodes found correlatable with their inputs"))
  
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
    measurements = log(data$stim_data[use, node_mes] / data$unstim_data[1, node_mes] )
    regression = paste0('lm(measurements[,1] ~ ')
    first = TRUE
    node_name = model_structure$names[node]
    condition = paste(node_name, "and")
    for (sender in upstreams[[node]] ) {
      mes_index = which(measured_nodes==sender)
      measurements = cbind(measurements, log(data$stim_data[use, mes_index] / data$unstim_data[1,mes_index]) )
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
    print(model_structure$names)
    print(params_matrix)
  }
  
  # Each link identified by correlation will correspond to one parameter
  # We identify which parameter it is and return it with its value
  params_vector = model$getParameterFromLocalResponse( params_matrix, rep(0, times=ncol(expdes$inhibitor)) )
  correlated = list()
  correlated$list = c()
  correlated$values = c()
  links = model$getParametersLinks()
  for ( param in 1:length(params_vector) ) {
    if (!is.nan(params_vector[param]) && params_vector[param] != 0) {
      correlated$list = c(correlated$list, param)
      correlated$values = c(correlated$values, params_vector[param])
      print(paste0("Correlated link ", simplify_path_name(links[param]), " = ", params_vector[param]))
    }
  }
  return(correlated)
}

# Parallel initialisation of the parameters
parallel_initialisation <- function(model, data, samples, NB_CORES) {
  # Put the samples under list format, as mclapply only take "one dimension" objects
  parallel_sampling = list()
  nb_samples = dim(samples)[1]
  length(parallel_sampling) = nb_samples # Prevents the resizing of the list
  for (i in 1:nb_samples) {
    parallel_sampling[[i]] = samples[i,]
  }
  # Parallel initialisations
  # The number of cores used depends on the ability of the detectCores function to detect them
  fitmodel_wrapper <- function(params, data, model) {
    init = proc.time()[3]
    result = model$fitmodel(data, params)
    if (verbose) {
      write(paste(signif(proc.time()[3]-init, 3), "s for the descent, residual =", result$residuals), stderr())
    }
    return(result)
  }
  
  fitmodelset_wrapper <- function(params, data, model) {
    # TODO add keep_constant control
    if ( class(data) != "Rcpp_DataSet" ) { stop("MRAmodelSet require a fitmodel::DataSet object") }
    init = proc.time()[3]
    result = model$fitmodelset(data, params)
    if (verbose) {
      write(paste(signif(proc.time()[3]-init, 3), "s for the descent, residual =", result$residuals), stderr())
    }
    return(result)
  }
  
  get_parallel_results <- function(model, data, samplings, NB_CORES) {
    if (class(model) == "Rcpp_ModelSet") {
      parallel_results = mclapply(samplings, fitmodelset_wrapper, data, model, mc.cores=NB_CORES)
    } else {
      parallel_results = mclapply(samplings, fitmodel_wrapper, data, model, mc.cores=NB_CORES)
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
    print("Using the block version for the parallel initialisation")
    best_results = list()
    best_results$residuals = c()
    best_results$params = c()
    for (i in 1:(nb_samples %/% 10000)) {
      results = get_parallel_results(model, data, parallel_sampling[ (10000*(i-1)+1):(10000*i) ], NB_CORES)
      
      # Only keep the best fits
      best = order(results$residuals)[1:20]
      write(paste(sum(is.na(results$residuals)), "NAs"), stderr())
      best_results$residuals = c(best_results$residuals, results$residuals[best])
      best_results$params = rbind(best_results$params, results$params[best,])
      # We make it so that the size of best_results never execeeds 10000
      if (i %% 500 == 0) {
        best = order(best_results$residuals)[1:20]
        best_results$residuals = best_results$residuals[best]
        best_results$params = best_results$params[best,]
      }
    }
    # The last block is smaller
    if (nb_samples %% 10000 != 0) {
      results = get_parallel_results(model, data, parallel_sampling[(10000 * nb_samples %/% 10000):nb_samples], NB_CORES)
      
      best = order(results$residuals)[1:20]
      best_results$residuals = c(best_results$residuals, results$residuals[best])
      best_results$params = rbind(best_results$params, results$params[best,])
    }
  } else {
    best_results = get_parallel_results(model, data, parallel_sampling, NB_CORES)
  }
  
  return(best_results)
}

# Initialise the parameters with a one core processing
# needs corrections, not used anyway
classic_initialisation <- function(model, data, samples) {
  for (i in 1:nb_samples) {
    result = model$fitmodel( data, samples[i,] )
    residuals = c(residuals,result$residuals)
    params = cbind(params,result$parameter)
    if (i %% (nb_samples/20) == 0) {
      print(paste(i %/% (nb_samples/20), "/ 20 initialization done."))
    }
  }
  if (debug) {
    print(sort(residuals))
  }
  results = list()
  results$residuals = resisuals
  results$params = params
  return(results)
}

# Detect the format of the structure file and extract the structure of the network
extractStructure = function(model_links, names="") {
  if (class(model_links)=="character") {
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
      name = unlist(strsplit(model_links, "/"))
      name = name[length(name)]
  } else {
      struct_matrix = model_links
      if (names == "") { name = "unknow" }
  }
  
  # Detect if it is a list of links or an adjacency matrix
  if (ncol(struct_matrix) == 2) {
    links_list = struct_matrix
  } else if (ncol(struct_matrix) == 3) {
    links_list = struct_matrix[,1:2] # Values are given, do not use them
  } else {
    # Remove the number of nodes
    if ( length(unique(struct_matrix[1,])) == 1 && struct_matrix[1, 1] != 1) {
      struct_matrix = struct_matrix[,-1]
    }
    # Set the column names if they are present, otherwise use the names from another file, and in last resort use numbers
    if (suppressWarnings( is.na(as.numeric(struct_matrix[1,1])) )) {
      colnames(struct_matrix) = struct_matrix[1,]
      struct_matrix = struct_matrix[-1,]
      if (suppressWarnings( is.na(as.numeric(struct_matrix[1,1])) )) {
        struct_matrix = struct_matrix[,-1]
      }
    } else if (length(names) == ncol(struct_matrix)) {
      colnames(struct_matrix) = names
    } else if (is.null(colnames(struct_matrix))) {
      print("No names were provided for the nodes, using numbers instead")
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
  
  model_structure=getModelStructure(links_list)
  
  return(model_structure)
}

#' Plot a graph from an adjacency list
#'
#' @param links_list A 2-columns matrix or a ModelStructure object. The network as an adjacency list, the first column is the upstream nodes, the second column the downstream nodes. Or a ModelStructure object as returned by getModelStructure.
#' @param expdes An ExperimentalDesign object. The measured, stimulated and inhibited nodes are highlighted if present.
# @export
plotNetworkGraph <- function(links_list, expdes="", local_values="", edge_lim = c(1, 10, 1)) {
    if (class(links_list) == "matrix") {
        names = unique(as.vector(links_list))
        adm=matrix(0,length(names),length(names),dimnames = list(names,names))
        for (ii in 1:nrow(links_list)) {
            adm[match(links_list[ii,2],rownames(adm)), links_list[ii,1]] = 1
        }
    } else if (class(links_list) == "Rcpp_ModelStructure") {
        adm = links_list$adjacencyMatrix
        colnames(adm) = rownames(adm) = links_list$names
    } else {
        stop("Invalid 'links_list' in plotNetworkGraph, must be an edge list or a ModelStructure")
    }
    g1 <- graphAM(adjMat=t(adm),edgemode="directed")
    nodeRenderInfo(g1) <- list(shape="ellipse")
    g1 <- layoutGraph(g1)
    edgeRenderInfo(g1) <- list(fontsize=10)
    # Add the experimental setup if provided
    if (class(expdes) == "Rcpp_ExperimentalDesign") {
      nodeRenderInfo(g1)$fill[1+expdes$measured_nodes] = "#ffff66"
      nodeRenderInfo(g1)$lwd = 1 # Create the field
      nodeRenderInfo(g1)$lwd[1+c(expdes$inhib_nodes, expdes$stim_nodes)] = 4 # Populate for perturbations
      nodeRenderInfo(g1)$col[1+expdes$inhib_nodes] = "red"
      nodeRenderInfo(g1)$col[1+expdes$stim_nodes] = "blue"
    }
    if (suppressWarnings(local_values != "")) {
        edgeRenderInfo(g1)$lwd = 1
        edgeRenderInfo(g1)$label = ""
        edgeRenderInfo(g1)$labelJust = 3
        from = edgeRenderInfo(g1)$enamesFrom
        to = edgeRenderInfo(g1)$enamesTo
        nodeX = nodeRenderInfo(g1)$nodeX
        nodeY = nodeRenderInfo(g1)$nodeY
        ii = 1
        for (vv in local_values$local_response) {
            if (vv != 0) {
                edgeRenderInfo(g1)$lwd[ii] = min(edge_lim[2], max(edge_lim[3] * abs(vv), edge_lim[1])) # Width between limits
                if (vv < 0) { edgeRenderInfo(g1)$col[ii] = "red" }
                edgeRenderInfo(g1)$label[ii] = signif(vv, 2)
                edgeRenderInfo(g1)$labelX[ii] = (nodeX[from[ii]] + nodeX[to[ii]]) / 2
                edgeRenderInfo(g1)$labelY[ii] = (nodeY[from[ii]] + nodeY[to[ii]]) / 2
                ii = ii + 1
            }
        }
    }
    renderGraph(g1)
    invisible(g1)
}

#' Extracts the data, the experimental design and the structure from the input files
#' @param model_structure Matrix of links [node1, node2]
#' @param data_filename Experimental data file name. Should be as follows, with one line per replicate:
#'
#'          stimulator                |          inhibitor                |                         type                       | [one column per measured nodes]
#'
#' -------------------------------------------------------------------------------------------------------------------------------------------------------------
#' stimulator name or solvant if none | inhibitor name or solvant if none | c for control, b for blank and t for experiment |    measure for the condition
extractModelCore <- function(model_structure, basal_activity, data_filename, var_file="") {
  ### READ DATA
  print(paste0("Reading data from ", data_filename))
  # Read the experiment design and extract the values
  use_midas = FALSE
  if (grepl(".data$", data_filename)) {
    data_file = read.delim(data_filename)
    data_file[data_file=="Medium"] = "DMSO"
    begin_measure = 4
    # Indicate where the conditions are
    conditions = c(1, 2)
    data.values = data_file[, colnames(data_file) %in% model_structure$names]
    not_included = (colnames(data_file)[!(colnames(data_file) %in% model_structure$names)])[-(1:(begin_measure-1))]
  } else if (grepl(".csv$", data_filename)) {
    use_midas = TRUE
    data_file = read.delim(data_filename, sep=",")
    if (ncol(data_file) <= 1) {
      data_file = read.delim(data_filename, sep="\t") # Accept tab-delimited MIDAS files
    }
    begin_measure = which(grepl("^DA.", colnames(data_file)))
    data_file = data_file[-begin_measure] # Delete the DA field which is not used
    begin_measure = begin_measure[1] # If there were several DA fields
    # Indicate where the conditions are
    conditions = 2:(begin_measure-1)
    
    # Extract the measurements of nodes in the network
    data.values = data_file[,grepl("^DV", colnames(data_file))]
    colnames(data_file) = gsub("^[A-Z]{2}.", "", colnames(data_file))
    colnames(data.values) = gsub("^[A-Z]{2}.", "", colnames(data.values))
    not_included = colnames(data.values)[!(colnames(data.values) %in% model_structure$names)]
    data.values = data.values[, colnames(data.values) %in% model_structure$names]
  } else {
    stop("Incorrect format for the data file, check the extension")
  }
  # Warn for the measured nodes that have not been found in the network
  if (length(not_included) > 0) {
    print(paste(not_included , "measurement is not in the network structure (could be a mispelling or a case error)" ))
  }
  
  # Means of basal activity of the network and of the blank fixation of the antibodies
  unstim.values = colMeans(data.values[data_file$type=="c"|data_file$type=="control",])
  lapply(unstim.values, function(X) { if (is.nan(X)|is.na(X)) {stop("Unstimulated data are required to simulate the network")}})
  blank.values = colMeans(data.values[data_file$type=="blank",])
  blank.values[is.nan(blank.values)] = 0; # For conditions without blank values
  
  # Calculates the mean and standard deviation for each condition 
  mean.values = aggregate(as.list(data.values),by=data_file[,1:(begin_measure-1)],mean, na.rm=T)
  sd.values = aggregate(as.list(data.values),by=data_file[,1:(begin_measure-1)],sd, na.rm=T)
  if (verbose) {
    print("Data used :")
    print(mean.values)
  }
  
  # Separate values and perturbation
  data.stim = mean.values[mean.values$type=="t",begin_measure:dim(mean.values)[2]]
  data.perturb = mean.values[mean.values$type=="t", conditions]
  
  ### CALCULATE ERROR MODEL
  if (grepl("\\.cv$", var_file) || grepl("\\.var$", var_file)) {
    # We use the CV file if there is one
    # The format and the order of the conditions are assumed to be the same as the data file
    print(paste0("Using var file ", var_file))
    if (use_midas) {
      variation.file = read.delim(var_file, sep=",")
      pre_cv = variation.file[, grepl("^DV", colnames(variation.file))]
      colnames(pre_cv) = gsub("^[A-Z]{2}.", "", colnames(pre_cv))
      # Check that the number of samples is the same for the measurements and the variation
      if (nrow(pre_cv) != nrow(data_file)) { stop("Different number of experiments for the variation and the measurement files") }
      # Check that the names in the measurements file and in the variation file match
      matchError = FALSE
      for (name in colnames(pre_cv)) {
        if (!(name %in% colnames(data_file))) {
          matchError = TRUE
          write(paste0(name, " from the variation file is not in the measurement file"), stderr())
        }
      }
      if (matchError) { stop("Names of the variation and measurement files do not match") }
      cv.values = aggregate(as.list( pre_cv[colnames(pre_cv) %in% model_structure$names] ), by=data_file[,1:(begin_measure-1)], mean, na.rm=T)
    } else {
      variation.file = read.delim(var_file)
      cv.values = aggregate(as.list(variation.file[, colnames(data.values) %in% model_structure$names]),by=data_file[,1:(begin_measure-1)], mean, na.rm=T)
    }
    cv.stim = cv.values[cv.values$type=="t", begin_measure:dim(cv.values)[2]]
    error = matrix(rep(blank.values,each=dim(data.stim)[1]),nrow=dim(data.stim)[1]) + cv.stim * data.stim
  } else {
    # Define the lower and default error threshold
    min.cv=0.1;     # parameters: minimal cv=0.1
    default.cv=0.3; # parameters: default cv if there are only 2 replicates
    
    # Calculate error percentage
    cv.values = sd.values[begin_measure:dim(sd.values)[2]] / mean.values[begin_measure:dim(sd.values)[2]]
    # Values too close to the blank are removed because the error is not due to antibody specific binding
    cv.values[!mean.values[,begin_measure:dim(mean.values)[2]] > 2 * matrix(rep(blank.values,each=dim(mean.values)[1]), nrow=dim(mean.values)[1])] = NA
    cv.stim = sd.values[sd.values$type=="t", begin_measure:dim(sd.values)[2]] # NA if no replicates are in stimulated data
    
    # Generation of error percentage, one cv per antibody calculated using all the replicates available, default.cv if there is only two replicate to calculate the cv
    cv = colMeans(cv.values,na.rm=TRUE)
    cv[cv<min.cv] = min.cv
    cv[is.nan(cv)|is.na(cv)]=default.cv
    
    if (FALSE) { #"Multiline comment"
      for (i in 1:dim(cv.values)[2]) {
        count = 0
        for (j in 1:dim(cv.values)[1]) {
          if (!is.na(cv[i][j]) | !is.nan(cv[i][j])) {
            count = count + 1
          }
        }
        if (count <= 2) {
          cv[i] = default.cv
          print("Defaulted")
        }
      }
    } ##
    
    if (verbose) {
      print("Error model :")
      for (i in 1:length(cv)) {
        print(paste(colnames(data.values)[i], " : ", cv[i]))
      }
    }
    error = matrix(rep(blank.values,each=dim(data.stim)[1]),nrow=dim(data.stim)[1])+matrix(rep(cv,each=dim(data.stim)[1]),nrow=dim(data.stim)[1])*data.stim
  }
  error[error<1] = 1 # The error cannot be 0 as it is used for the fit. If we get 0 (which means blank=0 and stim_data=0), we set it to 1 (which mean the score will simply be (fit-data)^2 for those measurements). We also ensure that is is not too small (which would lead to a disproportionate fit attempt
  
  ### SET UP DATA OBJECT
  
  data=new(fitmodel:::Data)
  data$set_unstim_data (matrix(rep(unstim.values,each=dim(data.stim)[1]),nrow=dim(data.stim)[1]))
  data$set_scale( data$unstim_data )
  data$set_stim_data( as.matrix(data.stim) )
  data$set_error( as.matrix(error ))
  # TODO data$data_cv
  
  ### EXTRACT EXPERIMENTAL DESIGN
  
  # Extraction of stimulated, inhibited and measured nodes
  if (use_midas) {
    names = colnames(mean.values)[conditions]
    stim_names = names[grepl("[^i]$", names)]
    stim.nodes = as.character( stim_names[ stim_names %in% model_structure$names] )
    names = gsub("i$", "", names[grepl("i$", names)])
    inhib.nodes = as.character( names[names %in% model_structure$names] )
  } else {
    stim.nodes = as.character(unique(mean.values$stimulator[mean.values$stimulator %in% model_structure$names]))
    inhib.nodes = as.character(unique(mean.values$inhibitor[mean.values$inhibitor %in% model_structure$names]))
  }
  measured.nodes=colnames(data.stim)
  
  # Inhibition and stimulation vectors for each experiment
  if (use_midas) {
    stimuli = as.matrix(data.perturb[stim.nodes])
  } else {
    stimuli=matrix(0,ncol=length(stim.nodes),nrow=dim(data.perturb)[1])
    for (i in 1:length(stim.nodes)) {
      stimuli[grepl(stim.nodes[i], data.perturb$stimulator),i]=1
    }
  }
  if (verbose) {
    print("Stimulated nodes")
    print(stim.nodes)
    print(stimuli)
  }
  if (use_midas) {
    inhibitor = as.matrix(data.perturb[paste(inhib.nodes, "i", sep="")])
  } else {
    inhibitor=matrix(0,ncol=length(inhib.nodes),nrow=dim(data.perturb)[1])
    if (length(inhib.nodes) > 0) { # Usefull for artificial networks
      for (i in 1:length(inhib.nodes)) {
        inhibitor[grepl(inhib.nodes[i], data.perturb$inhibitor),i]=1
      }
    }
  }
  if (verbose) {
    print("Inhibited nodes")
    print(inhib.nodes)
    print(inhibitor)
  }
  
  # Experimental design
  expdes=getExperimentalDesign(model_structure,stim.nodes,inhib.nodes,measured.nodes,stimuli,inhibitor,basal_activity)
  
  core = list()
  core$design = expdes
  core$data = data
  core$structure = model_structure
  core$basal = expdes$basal_activity
  core$cv = as.matrix(cv.stim)
  
  return(core)
}

#' Build a fitted model from a .mra file, and import data for this model
#' Does NOT perform any initialisation
#' @param model_file A .mra file containing the information on the model
#' @param data_file A .csv file with the data for the model
#' @param var_file A .var file with the variation of the data
#' @return An MRAmodel object describing the model and its best fit, containing the data
#' @export
#' @seealso importModel, exportModel, createModel
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
#' @examples
#' rebuildModel("model.mra", "data.csv", "data.var")
rebuildModel <- function(model_file, data_file, var_file="") {
  if (!grepl(".mra$", model_file)) {
    stop("The model file does not have the mra extension")
  }
  model = importModel(model_file)
  #links = matrix(rep(model$structure$names, 2), ncol=2)
  core = extractModelCore(model$structure, model$basal, data_file, var_file)
  
  model$data = core$data
  model$bestfit = model$model$fitmodel(model$data, model$parameters)$residuals
  model$cv = core$cv
  
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
  print(paste("Mutations per point =", nSamples))
  
  # The first iteration does not use any shift
  kept = matrix(0, ncol=model$nr_of_parameters())
  correlation = correlate_parameters(model, core, perform_plot)
  for (i in 1:depth) {
    print(paste("Depth :", i))
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
