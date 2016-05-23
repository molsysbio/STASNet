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
#' @param basefit The chi-2 value associated with the raw data (no fit)
#' @param name Name of the model
#' @param infos Extra information on the model
#' @param param_range Alternative parameters sets for the model
#' @param lower_value Lower bound on the values of the parameters
#' @param upper_values Upper bound on the values of the parameters
#' @return An MRAmodel object
#' @seealso \code{\link{createModel}}
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
MRAmodel <- function(model, design, structure, basal=matrix(), data=matrix(), cv=matrix(), parameters=vector(), bestfit=NA, name="", infos=c(), param_range=list(), lower_values=c(), upper_values=c()) {

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
                   upper_values=upper_values
                   ),
              class="MRAmodel")
    mra_model = computeFitScore(mra_model)
    return(mra_model)
}

#' Compute the fitting scores of a model
#'
#' Compute 2 fitting scores for the model, the fraction of the variation in the data explained by the network, and the improvement compared to a model with no links
#' Do the computation for each measured node and for the network
computeFitScore <- function(mra_model, refit_model=F) {
    data = mra_model$data
# The code for ModelSet::predict in C++ generates a segfault on datax return to R for an unknown reason
# Couldn't find the bug so we do not compute the score for the MRAmodelSet objects
    if (class(data) != "Rcpp_Data" || any(dim(data$stim_data)==0)) {# && class(data) != "Rcpp_DataSet") {
        mra_model$Rscores = NA
        mra_model$bestfitscore = NA
        return(mra_model)
    }
    if (refit_model) {
        refit = parallel_initialisation(mra_model$model, mra_model$data, matrix(mra_model$parameters, nrow=1), NB_CORES=1)
        mra_model$bestfit = refit$residual[1]
        mra_model$parameters = refit$params[1,]
    }
    prediction = getSimulation(mra_model)
    Rscores = c()
    meanScores = c()
    for ( abc in 1:ncol(prediction) ) {
        mdata = mean(data$stim_data[,abc])
        Sbase = sum((data$stim_data[,abc]-mdata)^2)
        Sfit = sum((data$stim_data[,abc]-prediction[,abc])^2)
        Rscores[colnames(prediction)[abc]] = 1 - Sfit/Sbase
    }

    mra_model$Rscores = Rscores
    mra_model$bestfitscore = mean(Rscores)

    return(mra_model)
}

getMeasuredNodesNames <- function(mra_model) {
    return(mra_model$structure$names[mra_model$design$measured_nodes+1])
}

#' Plot the network used in the model
#'
#' Plot the network used in the model with the experimental design on top of it
#' @param mra_model An MRAmodel object.
#' @export
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
