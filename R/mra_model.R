############################ mra_model.R ################################
# Definition of the MRAmodel object, an R wrapper around the C++ Model object which also contains the data necessary to fit the model

#' Constructor for MRAmodel objects
#'
#' Construct an MRAmodel objects, which contains all the informations recquired to simulate the model, and the data used to build it.
#' The model is built using perturbation datas and an input network.
#'
#' @param model An object of class Model
#' @param design An object of class ExperimentalDesign
#' @param structure An object of class ModelStructure
#' @param basal A matrix describing the basal activity of the nodes in the network
#' @param data An object of class Data
#' @param cv A matrix containing the coefficient of variation of the data used to build the model
#' @param parameters A vector containing the values of the parameters of the model
#' @param bestfit The residual chi-2 value associated with the best fit
#' @param name Name of the model
#' @param infos Extra information on the model
#' @param param_range Alternative parameters sets for the model
#' @param lower_values Lower bound on the values of the parameters
#' @param upper_values Upper bound on the values of the parameters
#' @param unused_perturbations Perturbations from the data file that have not been used for the fitting
#' @param unused_readouts Readouts from the data file that have not been used for the fitting
#' @return An MRAmodel object
#' @seealso \code{\link{createModel}}
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
MRAmodel <- function(model, design, structure, basal=matrix(), data=matrix(), cv=matrix(), parameters=vector(), bestfit=NA, name="", infos=c(), param_range=list(), lower_values=c(), upper_values=c(), unused_perturbations=c(), unused_readouts=c(), min_cv=0.1, default_cv=0.3) {

    if (class(data) != "Rcpp_Data" && class(data) != "Rcpp_DataSet") { stop("A Data object with unstimulated measurements is required") }
    mra_model = structure(
              list(
                   # Objects to build the model
                   model=model,
                   design=design,
                   structure=structure,
                   basal=basal,
                   data=data,
                   cv=cv,
                   # Optimal parameters of the model
                   parameters=parameters,
                   bestfit=bestfit,
                   # Name of the model and extra informations
                   name=name,
                   infos=infos,
                   # Values defined by profile likelihood
                   param_range=param_range,
                   lower_values=lower_values,
                   upper_values=upper_values,
                   unused_perturbations = unused_perturbations,
                   unused_readouts = unused_readouts,
                   min_cv = min_cv,
                   default_cv = default_cv
                   ),
              class="MRAmodel")
    mra_model = computeFitScore(mra_model)
    return(mra_model)
}

#' Compute the fitting scores of a model
#'
#' Compute 2 fitting scores for the model, the fraction of the variation in the data explained by the network, and the improvement compared to a model with no links
#' Do the computation for each measured node and for the network
#' @param mra_model The MRAmodel object for which the score should be computed
#' @param refit_model Whether the model should be refitted before computing the scores (using the 'mra_model$parameters' as the initial value)
#' @return A MRAmodel object with the scores in the fields 'Rscores' and 'bestfitscore'
computeFitScore <- function(mra_model, refit_model=FALSE) {
    data = mra_model$data
# The code for ModelSet::predict in C++ generates a segfault on datax return to R for an unknown reason
# Couldn't find the bug so we do not compute the score for the MRAmodelSet objects
    if (class(data) != "Rcpp_Data" || any(dim(data$stim_data)==0)) {# && class(data) != "Rcpp_DataSet") {
        mra_model$Rscores = NA
        mra_model$bestfitscore = NA
        return(computeReducedChiScore(mra_model))
    }
    if (refit_model) {
        refit = parallel_initialisation(mra_model$model, mra_model$data, matrix(mra_model$parameters, nrow=1), NB_CORES=1)
        mra_model$bestfit = refit$residual[1]
        mra_model$parameters = refit$params[1,]
    } else {
        simulation = simulateModel(mra_model)
        mra_model$bestfit = sum( (simulation$bestfit - simulation$data)^2/simulation$error^2, na.rm=T )
    }
    prediction = getSimulation(mra_model)
    Rscores = c()
    meanScores = c()
    for ( abc in 1:ncol(prediction) ) {
        mdata = mean(data$stim_data[,abc], na.rm=T)
        Sbase = sum((data$stim_data[,abc]-mdata)^2, na.rm=T)
        Sfit = sum((data$stim_data[,abc]-prediction[,abc])^2, na.rm=T)
        Rscores[colnames(prediction)[abc]] = 1 - Sfit/Sbase
    }

    mra_model$Rscores = Rscores
    mra_model$bestfitscore = mean(Rscores)
    mra_model = computeReducedChiScore(mra_model)

    return(mra_model)
}

computeReducedChiScore <- function(mra_model) {
    real_data = mra_model$data$stim_data
    data_count = sum(!is.na(real_data) & !is.nan(real_data))
    redChi = mra_model$bestfit / (data_count - length(mra_model$parameters))

    mra_model$reducedChi = redChi
    return(mra_model)
}

getMeasuredNodesNames <- function(mra_model) {
    return(mra_model$structure$names[mra_model$design$measured_nodes+1])
}

#' Plot the network used in the model
#'
#' Plot the network used in the model with the experimental design on top of it
#' @param mra_model A MRAmodel object.
#' @export
#' @family Model plots
#' @family Network graph
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
plotModelGraph <- function(mra_model) {
    plotNetworkGraph(mra_model$structure, mra_model$design, mra_model$model$getLocalResponseFromParameter(mra_model$parameters))
}

#' Return the paths corresponding to each parameter
#'
#' @param mra_model An MRAmodel object.
#' @return A vector containing the names of the parameters
#' @export
getParametersNames <- function(mra_model) {
    return(sapply(mra_model$model$getParametersLinks(), simplify_path_name))
}

#' Build infos vector from the parameters used to build the model and 
#' @param fit_info A vector of strings with additionnal informations
generate_infos <- function(input_file, inits, best_resid, method, model_links, name, fit_info=NULL) {
    infos = list()
    infos$call = sys.call(-1)
    if (!is.string(model_links)) {
        model_links = "R data"
    }
    infos$infos = c(paste0(inits, " samplings"), paste0( "Best residuals : ", paste0(best_resid, collapse=" ") ), paste0("Method : ", method), paste0("Network : ", model_links), fit_info)
    if ( !is.matrix(input_file) && all(sapply(input_file, is.string)) ) {
      infos$name = gsub("\\.csv", "", gsub("_MIDAS", "", basename(input_file)))
    } else {
      infos$name = name
    }

    return(infos)
}
