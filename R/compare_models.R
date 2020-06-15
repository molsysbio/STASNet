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
    if (is.list(models) && !is.null(names(models))) {
        names = names(models)
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
#' @param main Title of the plot
#' @param ... Extra parameters accepted by base::barplot
#' @return Invisibly, the residuals (chi-2)
#' @author Bertram Klinger \email{bertram.klinger@@charite.de}
#' @family Model plots
#' @export
plotResiduals <- function(model_group, main = "Model residuals", col="black", ...) {
  barx=barplot(model_group$residuals, ylab="residual", main=main, col=col, ...)
  text(barx, 1+0.1*min(model_group$residuals), labels = sapply(model_group$residuals, function(x) ifelse(abs(log10(x))>5, signif(x,1), round(x))), pos = 3, srt = 90, col = "grey80", offset = 1)
  invisible(model_group$residuals)
}

#' Plot the fitting scores of the model
#'
#' Generic function to plot scores for measured nodes fitting
#' @param x A Model object
#' @param ... Passed to base::barplot
#' @export
#' @seealso plotScores.modelGroup, plotScores.MRAmodel
#' @family Model plots
#' @rdname plotModelScores
plotModelScores <- function(x, ...) { UseMethod("plotModelScores", x) }
#' Plot the measurement scores for a group of models
#'
#' Plot the scores for each measured node fitting as a barplot with a bar for each model for each measurement
#' @param model_group A modelGroup object
#' @param ... Passed to base::barplot with a matrix of dimension nr_of_models * nr_of_measurements
#' @return See \code{\link{base::barplot}}
#' @rdname plotModelScores
#' @export
plotModelScores.modelGroup <- function(model_group, ...) {
    table = c()
    for (ii in 1:length(model_group$models)) {
        table = rbind(table, model_group$models[[ii]]$Rscores)
    }
    return(barplot(table, beside=TRUE, ...))
}

#' Plot parameters
#'
#' Generic function to plot parameters
#' @param x A Model object
#' @param ... Plot parameters
#' @export
#' @family Model plots
#' @rdname plotParameters
#' @seealso plotParameters.modelGroup
plotParameters <- function(x, ...) { UseMethod("plotParameters", x) }
#' Plot the parameters of the models
#'
#' Plot the parameters fitted for each model as an heatmap
#' @param model_group A modelGroup object
#' @export
#' @rdname plotParameters
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
