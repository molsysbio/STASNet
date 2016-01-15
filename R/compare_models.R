# Plot functions to compare models

#' Constructor for modelGroup objects
#'
#' Construct a modelGroup object, which contains several mraModels and provides methods to compare them
#' @param models A list of models to use in the modelGroup object
#' @return A modelGroup object
#' @export
modelGroup <- function(models=list()) {
    # Case when there is only one model in the list
    if (class(models) == "MRAmodel") { models = list(models)
    } else if (is.character(models) && !is.vector(models)) { models = list(importModel(models)) }
    if (!is.vector(models) || length(models)<1) {
        stop("You must provide a list of models to build a model group")
    }

    count = 0
    residuals = c()
    scores = c()
    links = c()
    names = c()
    models_list = list(); models_list[[length(models)]] = NA
    for (model in models) {
        count = count + 1
        if (is.character(model)) {
            model = importModel(model)
        } else if (class(model) != "MRAmodel") {
            stop(paste("Object", count, "of the list is not an MRAmodel"))
        }
        models_list[[count]] = model
        residuals = c(residuals, model$bestfit)
        scores = c(scores, model$bestfitscore)
        links = unique(c(links, model$model$getParametersLinks()))
        names = c(names, model$name)
    }
    names(scores) = names
    names(residuals) = names
    for (ii in 1:length(links)) {
        links[ii] = simplify_path_name(links[ii])
    }
    links = unique(links) # The product can have different forms in the non simplified form
    return(structure(
              list(
                   models = models_list,
                   residuals = residuals,
                   scores = scores,
                   links = links,
                   names = names
                   ),
              class="modelGroup"
              ))
}

#' Plot the residual of each model in the group
#'
#' Display the residuals of the models as a barplot
#' @param model_group Object of class modelGroup
#' @return Invisibly, the residuals (chi-2)
#' @export
plotResiduals <- function(model_group) {
    barplot(model_group$residuals, main="Models residuals")
    invisible(model_group$residuals)
}

#' Plot scores
#'
#' Generic function to plot scores
#' @export
#' @seealso plotScores.modelGroup
plotScores <- function(x, ...) { UseMethod("plotScores", x) }
#' Plot the scores of each model in the group
#'
#' Display the scores of the modelGroup object as a barplot
#' @param model_group Object of class modelGroup
#' @return Invisibly, the scores ( -log(best_fit/base_fit) )
#' @seealso \code{\link{MRAmodel}}
#' @export
plotScores.modelGroup <- function(model_group) {
    scores = -log(model_group$scores)
    barplot(scores, main="Models scores")
    invisible(scores)
}

#' Plot parameters
#'
#' Generic function to plot parameters
#' @export
#' @seealso plotParameters.modelGroup
plotParameters <- function(x, ...) { UseMethod("plotParameters", x) }
#' Plot the parameters of the models
#'
#' Plot the parameters fitted for each model as an heatmap
#' @export
#' @param model_group Object of class modelGroup
#' @export
plotParameters.modelGroup <- function(model_group) {
    all_links = model_group$links
    param_matrix = c()
    for (mraModel in model_group$models) {
        links = sapply(mraModel$model$getParametersLinks(), simplify_path_name)
        parameters = mraModel$parameters
        param_vec = c()
        for (link in all_links) {
            if (link %in% links) {
                param_vec = c(param_vec, parameters[which(links==link)])
            } else {
                param_vec = c(param_vec, NA)
            }
        }
        param_matrix = rbind(param_matrix, param_vec)
    }
    colnames(param_matrix) = all_links
    rownames(param_matrix) = model_group$names
    param_matrix = t(param_matrix)
    if (length(model_group$models) > 1) { p2 = T } else { p2=F }
# TODO Define breaks centered on 0
    return( pheatmap(param_matrix, cluster_rows=p2, cluster_cols=p2, display_numbers=T) )
}
