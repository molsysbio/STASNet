############################ model_set.R ################################
# The MRAmodelSet object, to fit multiple dataset on the same network

#' Import multiple files to create a group of models
#' @param model_links An adjacency matrix
#' @param data_list Name of the file containing the list of the files for the data without the extension
#' @param basal_file File with basal activity
#' @param cores Number of cores to use
#' @param inits Number of inits to do
#' @param init_distribution Whether the initial distribution should be plotted
#' @param method The LHs method to use
createDataSet <- function(model_links, data_list, basal_file, cores=1, inits=1000, init_distribution=F, method="default") {
    files = unique(gsub("\\..*$", "", readLines(data_list)))
    folder_files = dir()
    model_set = vector("list",length(files))
    ii = 1
    for (file in files) {
        if (paste0(file, ".var") %in% folder_files) {
            model_set[[ii]] = createModel(model_links, paste0(file, ".csv"), basal_file, paste0(file, ".var"), cores, inits, init_distribution, method)
        } else {
            model_set[[ii]] = createModel(model_links, paste0(file, ".csv"), basal_file, cores, inits, init_distribution, method)
        }
        ii = ii+1
    }
}

#' Compare networks
#'
#' Compile the parameters of a model set in one matrix and plot the results.
#' highlights the variable links by scaling to it
#' @param modelset A MRAmodelSet object
#' @return matrix of parameters
#' @export
#' @rdname compareParameters
#' @author Bertram Klinger \email{bertram.klinger@@charite.de}
compareParameters <- function(multi_model, ...) { UseMethod("compareParameters", multi_model) }

#' @export
#' @rdname compareParameters
compareParameters.MRAmodelSet <- function(modelset) {
  links=matrix(modelset$parameters,ncol=modelset$nb_models,byrow = F)
  colnames(links) <- modelset$names
  rownames(links) <- unname(sapply(modelset$model$getParametersLinks(), function(x) STASNet:::simplify_path_name(x)))

  if (length(modelset$variable_parameters)==0){
     warning("No variable links detected, please run 'addVariableParameters()' before calling this function!")
  }
  
  STASNet:::plotHeatmap(mat = links,
              main = "modelSet parameters rowwise scaled to mean",
              stripOut = 0.01,
              lim = 10,
              scale_rows = T)
    invisible(links)
}

#' @export
#' @rdname compareParameters
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
# Very naive implementation TODO compare parameter names and make them match
compareParameters.modelGroup <- function(modelgroup) {
    links = sapply(modelgroup$models, function(mm){ mm$parameters }) # Assumes same length
    if (typeof(links) == "list") { stop("Different number of parameters is not supported yet") }

    colnames(links) = modelgroup$names
    rownames(links) = getParametersNames(modelgroup[[1]])

    STASNet:::plotHeatmap(mat = links,
              main = "modelGroup parameters rowwise scaled to mean",
              stripOut = 0.01,
              lim = 10,
              scale_rows = T)
    invisible(links)
}

#' Constructor for MRAmodelSet objects
#'
#' Build an MRAmodelSet, which contains all the information required to simulate a set of MRAmodels, ensuring that the sets of parameters are similar for all models
#'
#' @param nb_models Number of models in the set
#' @inheritParams MRAmodel
#' @return An MRAmodelSet object
#' @seealso \code{\link{createModel}}
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
MRAmodelSet <- function(nb_models=1, model=NULL, design=NULL, structure=NULL, basal=matrix(), data=matrix(), cv=matrix(), parameters=vector(), bestfit=NA, name=c(), infos=c(), param_range=list(), lower_values=c(), upper_values=c(), unused_perturbations=c(), unused_readouts=c(), min_cv=0.1, default_cv=0.3, use_log=FALSE) {
    if (length(name) != nb_models) {
        name = rep(name[1], nb_models)
    }

    # An MRAmodelSet is an MRAmodel
    self = MRAmodel(model, design, structure, basal, data, cv, parameters, bestfit, paste0("Model set using: ", paste0(name, collapse=" ")),  infos, param_range, lower_values, upper_values, unused_perturbations, unused_readouts, min_cv, default_cv, use_log)
    # With some extra attributes
    class(self) = c("MRAmodelSet", class(self))
    self$nb_models = nb_models
    self$names = name
    self$variable_parameters = numeric(0)

    return(self)
}

#' Set the variable parameters of a modelset
#'
#' Set the parameters that can vary across submodels for an MRAmodelSet
#' @param modelset An MRAmodelSet object
#' @param parameters_ids The ids of the parameters to set as variable across submodels
#' @return The modified MRAmodelSet
#' @export
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
setVariableParameters <- function(modelset, parameters_ids) {
    nb_parameters = modelset$model$nr_of_parameters()/modelset$nb_models
    for (pid in parameters_ids) {
        if (pid > nb_parameters) { stop(paste("Invalid parameter id:", pid)) }
    }
    parameters_ids = sort(parameters_ids)
    modelset$variable_parameters = parameters_ids
    modelset$model$setVariableParameters(parameters_ids)
    return(modelset)
}

#' Extract individual models from an MRAmodelSet
#'
#' Extract individual models from an MRAmodelSet where they were fitted together and return them in a list where the indices are the names of each submodel.
#' @param modelset An MRAmodelSet object
#' @return A list of MRAmodel objects
#' @export
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
extractSubmodels <- function(modelset) {
    model_list = vector('list',modelset$nb_models)
    model = new(STASNet:::Model)
    model$setModel(modelset$design, modelset$structure, modelset$use_log)
    nb_parameters = length(modelset$parameters)/modelset$nb_models # All submodels have the same number of parameters
    data_size = nrow(modelset$data$unstim_data)/modelset$nb_models # The data matrix dimensions are the same for all models
    for (ii in 1:modelset$nb_models) {
        parameters = modelset$parameters[((ii-1)*nb_parameters+1):(ii*nb_parameters)]
        row_selection = ((ii-1)*data_size+1):(ii*data_size)
        cv = modelset$cv[row_selection,]
        data = modelset$data$datas_list[[ii]]
        fit_value = sum( (( model$simulate( data, parameters )$prediction - data$stim_data) / data$error)^2, na.rm=T)
        model_list[[ii]] = MRAmodel(model, modelset$design, modelset$structure, modelset$basal, data, cv, parameters, fit_value, modelset$names[ii], modelset$infos, modelset$param_range, modelset$lower_values, modelset$upper_values, modelset$unused_perturbations, modelset$unused_readouts, modelset$min_cv, modelset$default_cv, modelset$use_log)
    }
    return(modelGroup(model_list))
}
