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

#' Compute the fitting score of a model
#'
#' Compute the fitting score of the model (fraction of the variation in the data explained by the network)
#' Do the computation for each measured node and for the network
computeFitScore <- function(mra_model) {
    if (class(data) == "Rcpp_Data" || class(data) == "Rcpp_DataSet") {
        data = mra_model$data
    } else {
        mra_model$abScores = NA
        mra_model$bestfitscore = NA
        return()
    }
    refit = parallel_initialisation(mra_model$model, mra_model$data, matrix(mra_model$parameters, nrow=1), NB_CORES=1)
    mra_model$bestfit = refit$residual[1]
    mra_model$parameters = refit$params[1,]
    prediction = getSimulation(mra_model)
    Rscores = c()
    for ( abc in 1:ncol(prediction) ) {
        mdata = mean(data$stim_data[,abc])
        Sbase = sum((data$stim_data[,abc]-mdata)^2)
        Sfit = sum((data$stim_data[,abc]-prediction[,abc])^2)
        Rscores[colnames(prediction)[abc]] = 1 - Sfit/Sbase
    }

    mra_model$abScores = Rscores
    mra_model$bestfitscore = mean(Rscores)

    return(mra_model)
}

getMeasuredNodesNames <- function(mra_model) {
    return(mra_model$structure$names[mra_model$design$measured_nodes+1])
}
