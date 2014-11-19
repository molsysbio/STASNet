########################### model_manipulation.R ###########################
# Functions to change and visualise the model

#' Print the value of each path from the model, and add the profile likelihood infos if they are provided
#' @param model_description An MRAmodel object
#' @return Nothing
#' @export
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
printParameters <- function(model_description) {
    model = model_description$model
    parameters = model_description$parameters

    print("Parameters :")
    paths = model$getParametersLinks()
    decimals = 2
    if (length(model_description$lower_values) > 0) {
        for (i in 1:length(paths)) {
            text = paste(simplify_path_name(paths[i]), "=", signif(parameters[i],decimals))
            if (is.na(model_description$lower_values[i])) {
                if (is.na(model_description$upper_values[i])) {
                    text = paste(text, "(non identifiable)")
                } else {
                    text = paste(text, "( ni - ", signif(model_description$upper_values[i], decimals), ")")
                }
            } else {
                text = paste(text, "(", signif(model_description$lower_values[i], decimals))
                if (is.na(model_description$upper_values[i])) {
                    text = paste(text, "- ni )")
                } else {
                    text = paste(text, "-", signif(model_description$upper_values[i], decimals), ")")
                }
            }
            print(text)
        }
    } else {
        for (i in 1:length(paths)) {
            print (paste( simplify_path_name(paths[i]), ":", signif(parameters[i], decimals) ))
        }
    }
}

#' Plots heatmaps of the model prediction against the data, weighted by the error, as well as the log fold change data and the prediction
#' @param model_description A list describing the model, as the one produced by createModel or importModel
#' @return Nothing
#' @export
#' @seealso plotModelPrediction, createModel
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
accuracyPlot <- function(model_description) {
    # Calculate the mismatch
    model = model_description$model
    data = model_description$data
    error = data$error
    cv = model_description$cv
    stim_data = data$stim_data
    init_params = model_description$parameter

    simulation = model$simulate(data, init_params)$prediction
    if (!is.null(nrow(cv)) && nrow(cv) == nrow(stim_data)) {
        mismatch = (stim_data - simulation) / (cv * stim_data)
    } else {
        mismatch = (stim_data - simulation) / error
    }
    simulation = log(simulation / data$unstim_data, base=2)
    stim_data = log(stim_data / data$unstim_data, base=2)

    # Rebuild the conditions from the design
    nodes = model_description$structure$names
    design = model_description$design
    treatments = c()
    for (row in 1:nrow(mismatch)) {
        stim_names = nodes[design$stim_nodes[which(design$stimuli[row,]==1)]+1]
        inhib_names = nodes[design$inhib_nodes[which(design$inhibitor[row,]==1)]+1]
        if (length(inhib_names) > 0) {
            inhib_names = paste(inhib_names, "i", sep="")
        }
        treatments = c(treatments, paste(c(stim_names, inhib_names), collapse="+", sep="") )
    }
    print("Treatments : ")
    print(treatments)
    colnames(mismatch) = colnames(stim_data) = colnames(simulation) = nodes[design$measured_nodes + 1]
    rownames(mismatch) = rownames(stim_data) = rownames(simulation) = treatments

    # Comparison of the data and the stimulation in term of error fold change and log fold change
    bk = unique( c(seq(min(mismatch, na.rm=T), 0, length=50), seq(0, max(mismatch, na.rm=T), length=50)) )
    pheatmap(mismatch, color=colorRampPalette(c("deepskyblue", "black", "red"))(length(bk)-1), breaks = bk, cluster_rows=F, cluster_col=F, display_numbers = T, main="(data - simulation) / error")
    bk = unique( c(seq(min(stim_data-simulation, na.rm=T), 0, length=50), seq(0, max(stim_data-simulation, na.rm=T), length=50)) )
    pheatmap(stim_data-simulation, color=colorRampPalette(c("deepskyblue", "black", "red"))(length(bk)-1), breaks = bk, cluster_rows=F, cluster_col=F, display_numbers = T, main="log2(data/simulation)")
    # Log fold changes for the data and the stimulation
    bk = unique( c(seq(min(stim_data, na.rm=T), 0, length=50), seq(0, max(stim_data, na.rm=T), length=50)) )
    pheatmap(stim_data, color=colorRampPalette(c("deepskyblue", "black", "red"))(length(bk)-1), breaks = bk, cluster_rows=F, cluster_col=F, display_numbers = T, main="data")
    bk = unique( c(seq(min(simulation, na.rm=T), 0, length=50), seq(0, max(simulation, na.rm=T), length=50)) )
    pheatmap(simulation, color=colorRampPalette(c("deepskyblue", "black", "red"))(length(bk)-1), breaks = bk, cluster_rows=F, cluster_col=F, display_numbers = T, main="simulation")
}

#' Selection of a minimal model by the removal of non significant links with a Chi^2 test
#' @param model_description A list describing the model, as the one produced by createModel or importModel
#' @param accuracy Probability of the confidence interval for the Chi^2 test. The link can be removed if 0 is included in this interval.
#' @return An MRAmodel object of the reduced model with the data
#' @export
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
selectMinimalModel <- function(model_description, accuracy=0.95) {
    # Extra fitting informations from the model description
    model = model_description$model
    init_params = model_description$parameters
    initial_response = model$getLocalResponseFromParameter( init_params )
    expdes = model_description$design
    model_structure = model_description$structure
    adj = model_structure$adjacencyMatrix
    data = model_description$data
    if (is.na(model_description$bestfit)) {stop("Data are recquired to reduce the model")}

    print("Performing model reduction…")
    initresidual = model_description$bestfit
    rank = model$modelRank()
    reduce=TRUE
    while (reduce) {
        links.to.test=which(adj==1)
        params=c()
        residuals=c()
        ranks=c()
        newadj=adj

        # Each link is removed and the best of those networks is compared to the previous model
        newadj=adj
        for (i in links.to.test) {
            newadj[i]=0
            model_structure$setAdjacencyMatrix( newadj )
            model$setModel ( expdes, model_structure )
            paramstmp = model$getParameterFromLocalResponse(initial_response$local_response, initial_response$inhibitors)
            result = model$fitmodel( data,paramstmp)
            response.matrix = model$getLocalResponseFromParameter( result$parameter )
            residuals = c(residuals, result$residuals);   
            params = cbind(params, c(response.matrix))
            new_rank = model$modelRank()
            ranks = c(ranks, new_rank)

            if (verbose) {
                dr = rank - new_rank
                print(paste("old :", rank, ", new : ", new_rank))
                deltares = residuals[length(residuals)]-initresidual
                print(paste(model_structure$names[(i-1) %/% dim(adj)[1]+1], "->", model_structure$names[(i-1) %% dim(adj)[1]+1], ": Delta residual = ", deltares, "; Delta rank = ", dr, ", p-value = ", pchisq(deltares, df=dr) ))
            }

            newadj[i]=1; ## Slightly accelerate the computation
        }
        
        order.res=order(residuals)
        # The loss of degree of freedom is equal to the difference in the ranks of the matrices
        new_rank = ranks[order.res[1]]
        dr = rank - new_rank
        deltares = residuals[order.res[1]]-initresidual
        # Some boundary cases might give low improvement of the fit
        if (deltares < 0) { warning(paste("Negative delta residual :", deltares)) ; deltares = -deltares  }
        if (deltares < qchisq(accuracy, df=dr)) {
            adj[links.to.test[order.res[1]]]=0
            rank = new_rank
            initial_response=params[,order.res[1]]
            initresidual = residuals[order.res[1]]
            #print(initial_response)
            print(paste0("Remove ",
                      model_structure$names[((links.to.test[order.res[1]]-1) %/% (dim(adj)[1])) +1], "->", # Line
                      model_structure$names[((links.to.test[order.res[1]]-1) %% (dim(adj)[1])) +1])); # Column (+1 because of the modulo and the R matrices starting by 1 instead of 0)

            print(paste( "New residual = ", residuals[order.res[1]], ", Delta residual = ", deltares, ",  p-value = ", pchisq(deltares, df=dr) ))
            print("------------------------------------------------------------------------------------------------------")

            model_description$bestfit = sort(residuals)[1]
        } else {
          reduce=FALSE
        }
    }
    print("Reduction complete")
    # We recover the final model
    ## Basal activity and data do not change
    model_description$structure$setAdjacencyMatrix(adj)
    model_description$model$setModel(expdes, model_description$structure)
    model_description$parameters = model_description$model$getParameterFromLocalResponse(initial_response$local_response, initial_response$inhibitors)
    model_description$infos = c(model_description$infos, "Reduced model")
    model_description$param_range = list()
    model_description$lower_values = c()
    model_description$upper_values = c()

    return(model_description)
}



