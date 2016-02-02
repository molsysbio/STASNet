############################ model_set.R ################################
# The MRAmodelSet object, to fit multiple dataset on the same network

#' Import multiple files 
#' @param data_list Name of the file containing the list of the files for the data without the extension
createDataSet <- function(model_links, data_list, basal_file, cores=1, inits=1000, init_distribution=F, method="default") {
    files = unique(gsub("\\..*$", "", readLines(data_list)))
    folder_files = dir()
    model_set = list()
    model_set[[length(files)]] = ""
    for (file in files) {
        if (paste0(file, ".var") %in% folder_files) {
            model_set[[i]] = create_model(model_links, paste0(file, ".csv"), basal_file, paste0(file, ".var"), cores, inits, init_distribution, method)
        } else {
            model_set[[i]] = create_model(model_links, paste0(file, ".csv"), basal_file, cores, inits, init_distribution, method)
        }
    }
}

#' Compare networks
#'
#' Compare the parameters of one network for different conditions (cell line, perturbations, ...)
#' The links must be the same for all the networks
#' @param files A list of .mra files
#' @return None
# TODO improve the function so that different networks with common links can be compared
compareModels <- function(files) {
  models = list()
  for (i in 1:length(files)) { models[[i]] = importModel(files[i])

  links = c()
  for (model in models) { links = cbind(links, model$parameters) }
  rownames(links) = models[[1]]$model$getParametersLinks()
  colnames(links) = files  
  med = median(abs(links))
  m = max(abs(links))
  breaks = unique(c(seq(-m, -2*med, length.out = 10), seq(-2*med, 2*med, length.out=50), seq(2*med, m, length.out=10)))
  pheatmap(links, breaks = breaks, color=colorRampPalette("deepskyblue", "black", "red")(length(breaks)-1))
  }
}

#' Constructor for MRAmodelSet objects
#'
#' Build an MRAmodelSet, which contains all the information required to simulate a set of MRAmodels, ensuring that the sets of parameters are similar for all models
#'
#' @param nb_models Number of models in the set
#' @param model An object of class Model
#' @param design An object of class ExperimentalDesign
#' @param structure An object of class ModelStructure
#' @param basal A matrix describing the basal activity of the nodes in the network
#' @param data An object of class Data with aggregated data for all models
#' @param cv A matrix containing the coefficient of variation of the data used to build the model
#' @param parameters A vector containing the values of the parameters of the model
#' @param bestfit The residual chi-2 value associated with the best fit
#' @param basefit The chi-2 value associated with the raw data (no fit)
#' @param name Name of the model
#' @param infos Extra information on the model
#' @param param_range Alternative parameters sets for the model
#' @param lower_value Lower bound on the values of the parameters
#' @param upper_values Upper bound on the values of the parameters
#' @return An MRAmodelSet object
#' @seealso \code{\link{createModel}}
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
MRAmodelSet <- function(nb_models=1, model=NULL, design=NULL, structure=NULL, basal=matrix(), data=matrix(), cv=matrix(), parameters=vector(), bestfit=NA, names=c(), infos=c(), param_range=list(), lower_values=c(), upper_values=c()) {
    if (length(names) != nb_models) {
        names = rep(names[1], nb_models)
    }

    # An MRAmodelSet is an MRAmodel
    self = MRAmodel(model, design, structure, basal, data, cv, parameters, bestfit, paste0("Model set using: ", paste0(names, collapse=" ")),  infos, param_range, lower_values, upper_values)
    # With some extra attributes
    class(self) = c("MRAmodelSet", class(self))
    self$nb_models = nb_models
    self$names = names
    self$variable_parameters = c()

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
    model_list = list()
    model_list[[modelset$nb_models]] = NA
    model = new(fitmodel:::Model)
    model$setModel(modelset$design, modelset$structure)
    for (ii in 1:modelset$nb_models) {
        nb_parameters = length(modelset$parameters)/modelset$nb_models # All submodels have the same set of parameters
        parameters = modelset$parameters[((ii-1)*nb_parameters+1):(ii*nb_parameters)]
        data_size = nrow(modelset$data$unstim_data)/modelset$nb_models # The data matrix dimensions are the same for all models
        row_selection = ((ii-1)*data_size+1):(ii*data_size)
        cv = modelset$cv[row_selection,]
        data = modelset$data$datas_list[[ii]]
        fit_value = sum( (( model$simulate( data, parameters )$prediction - data$unstim_data) / data$error)^2, na.rm=T)
        model_list[[ii]] = MRAmodel(model, modelset$design, modelset$structure, modelset$basal, data, cv, parameters, fit_value, modelset$names[ii], modelset$infos, modelset$param_range, modelset$lower_values, modelset$upper_values)
    }
    return(modelGroup(model_list))
}
