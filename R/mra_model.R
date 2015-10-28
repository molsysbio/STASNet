
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
MRAmodel <- function(model=NULL, design=NULL, structure=NULL, basal=matrix(), data=matrix(), cv=matrix(), parameters=vector(), bestfit=NA, basefit=NA, name="", infos=c(), param_range=list(), lower_values=c(), upper_values=c()) {
    structure(
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
                   basefit = basefit, # Chi-2 score without any link
                   bestfitscore = bestfit / basefit,
                   # Name of the model and extra informations
                   name=name,
                   infos=infos,
                   # Values defined by profile likelihood
                   param_range=param_range,
                   lower_values=lower_values,
                   upper_values=upper_values
                   ),
              class="MRAmodel")
}

