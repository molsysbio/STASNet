########################### model_manipulation.R ###########################
# Functions to change and visualise the model

#' Print the value of the parameters of the model
#'
#' Print the value of each path from the model, with the profile likelihood infos if they are provided.
#' @param model_description An MRAmodel object
#' @param precision Number of significant digits to print
#' @return Nothing
#' @details The print function is 'message' and thus produces an output in stderr
#' @export
#' @seealso \code{\link{message}}
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
printParameters <- function(model_description, precision=2) {
  model = model_description$model
  parameters = model_description$parameters
  
  message("Parameters :")
  paths = model$getParametersLinks()
  if (length(model_description$lower_values) > 0) {
    for (i in 1:length(paths)) {
      text = paste(simplify_path_name(paths[i]), "=", signif(parameters[i],precision))
      if (is.na(model_description$lower_values[i])) {
        if (is.na(model_description$upper_values[i])) {
          text = paste(text, "(non identifiable)")
        } else {
          text = paste(text, "( ni - ", signif(model_description$upper_values[i], precision), ")")
        }
      } else {
        text = paste(text, "(", signif(model_description$lower_values[i], precision))
        if (is.na(model_description$upper_values[i])) {
          text = paste(text, "- ni )")
        } else {
          text = paste(text, "-", signif(model_description$upper_values[i], precision), ")")
        }
      }
      message(text)
    }
  } else {
    for (i in 1:length(paths)) {
      message (paste( simplify_path_name(paths[i]), ":", signif(parameters[i], precision) ))
    }
  }
}

#' Plot heatmaps of the model simulation against the data weighted by the error, as well as the log fold change for the data and the prediction
#' @param model_description An MRAmodel object
#' @param limit An integer to force the limit of the heatmaps
#' @param show_values Whether the values should be printed in the heatmap boxes or not.
#' @param graphs Define which graphs should be plotted. "accuracy" for the residual as seen by the model, "qq" for the qqnorm plot of those residuals, "diff" for the delta log data-simulation, "data" for the log-fold change computed from the data, "simulation" for the log-fold change simulated by the model, "prediction" for the log-fold change that would be predicted without the blank correction.
#' @param selected_treatments A vector with the names of the subset of treatments that should be plotted
#' @param selected_readouts A vector with the names of the subset of readouts that should be plotted
#' @param name The name of the model, used as subtitle in the plot
#' @param regroup A value in c('inhib', 'stim') which indicates if the data should be regrouped by inhibition or stimulation in the heatmap
#' @return Invisibly, a list with elements 'mismatch', 'stim_data' and 'simulation' corresponding to the values plotted.
#' @export
#' @seealso createModel, importModel
#' @family Model plots
#' @rdname accuracy_plot
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
plotModelAccuracy <- function(x, ...) { UseMethod("plotModelAccuracy", x) }
#' Plot model accuracy for MRAmodel
#' @export
#' @rdname accuracy_plot
plotModelAccuracy.MRAmodel <- function(model_description, limit=Inf, show_values=TRUE, graphs=c("accuracy", "diff", "data", "simulation", "prediction"), selected_treatments = c(), selected_readouts = c(), name="", regroup="no") {
  # Calculate the mismatch
  model = model_description$model
  data = model_description$data
  error = data$error
  cv = model_description$cv
  stim_data = data$stim_data
  init_params = model_description$parameter
  if (name == "") { name = model_description$name }
  
  simulation = model$simulateWithOffset(data, init_params)$prediction
  prediction = log2(model$simulate(data, init_params)$prediction / data$unstim_data)
  if (model_description$use_log) {
      mismatch = (log(stim_data) - log(simulation)) / (log(error)*sqrt(2))
  } else {
      mismatch = (stim_data - simulation) / (error*sqrt(2))
  }
  simulation = log2(simulation / data$unstim_data)
  stim_data = log2(stim_data / data$unstim_data)
  
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

  if (verbose > 3) {
      message("Treatments : ")
      message(paste(treatments, collapse=" "))
  }
  colnames(mismatch) = colnames(stim_data) = colnames(simulation) = colnames(prediction) = nodes[design$measured_nodes + 1]
  rownames(mismatch) = rownames(stim_data) = rownames(simulation) = rownames(prediction) = treatments

  if (regroup %in% c("stim", "inhib")) {
      if (regroup == "stim") {
          reorder = order(apply(cbind(design$stimuli, design$inhibitor), 1, paste0, collapse="_"), decreasing=TRUE)
      } else {
          reorder = order(apply(cbind(design$inhibitor, design$stimuli), 1, paste0, collapse="_"), decreasing=TRUE)
      }
      mismatch = mismatch[reorder,]
      stim_data = stim_data[reorder,]
      simulation = simulation[reorder,]
      prediction = prediction[reorder,]
  }

  # Subset the matrices to only have the select readouts and treatments
  if (length(selected_readouts) > 0) {
      valid_selection = c()
      for (sr in selected_readouts) {
          if (sr %in% colnames(stim_data)) {
              valid_selection = c(valid_selection, sr)
          }
      }
      mismatch = mismatch[,valid_selection,drop=FALSE]
      stim_data = stim_data[,valid_selection,drop=FALSE]
      simulation = simulation[,valid_selection,drop=FALSE]
      prediction = prediction[,valid_selection,drop=FALSE]
  }
  if (length(selected_treatments) > 0) {
      valid_selection = c()
      for (st in selected_treatments) {
          if (st %in% rownames(stim_data)) {
              valid_selection = c(valid_selection, st)
          }
      }
      mismatch = mismatch[valid_selection,,drop=FALSE]
      stim_data = stim_data[valid_selection,,drop=FALSE]
      simulation = simulation[valid_selection,,drop=FALSE]
      prediction = prediction[valid_selection,,drop=FALSE]
  }

# Comparison of the data and the stimulation in term of error fold change and log fold change
  if (any(grepl("qq", graphs))) {
      qqnorm(mismatch)
      mis_range = max(abs(range(mismatch, na.rm=TRUE)))
      lines(c(-mis_range, mis_range), c(-mis_range, mis_range), col="red")
  }
  if (any(grepl("acc", graphs)))
      plotHeatmap(mismatch,"(data - simulation) / error", show_values=show_values, lim=2, fixedRange=TRUE, sub=name)
  if (any(grepl("diff", graphs)))
      plotHeatmap(stim_data-simulation,"log2(data/simulation)", show_values=show_values, sub=name)
# Log fold changes for the data and the stimulation with comparable color code
  lim=min(10, max(abs( range(quantile(stim_data,0.05, na.rm=TRUE),
                       quantile(simulation,0.05, na.rm=TRUE),
                       quantile(stim_data,0.95, na.rm=TRUE),
                       quantile(simulation,0.95, na.rm=TRUE)) )))
  if (!is.infinite(limit)) { lim = limit }
  if (any(grepl("data", graphs)))
      plotHeatmap(stim_data, "Log-fold change Experimental data",lim,TRUE, show_values=show_values, sub=name)
  if (any(grepl("sim", graphs)))
      plotHeatmap(simulation, "Log-fold change Simulated data",lim,TRUE, show_values=show_values, sub=name)
  if (any(grepl("pred", graphs)))
      plotHeatmap(prediction, "Log-fold change Prediction",lim,TRUE, show_values=show_values, sub=name)

  invisible(list(mismatch=mismatch, stim_data=stim_data, simulation=simulation))
}
#' Plot accuracy of all submodels of a modelset
#' @export
#' @rdname accuracy_plot
plotModelAccuracy.MRAmodelSet <- function(model_description, limit=Inf, show_values=TRUE, graphs=c("accuracy", "data", "simulation"), selected_treatments = c(), selected_readouts = c(), name="") {
    submodels = extractSubmodels(model_description)
    if (any(grepl("acc", graphs))) {
        for (subm in submodels$models) { plotModelAccuracy(subm, graphs="accuracy") }
    }
    if (any(grepl("diff", graphs))) {
        for (subm in submodels$models) { plotModelAccuracy(subm, graphs="diff") }
    }
    if (any(grepl("data", graphs))) {
        for (subm in submodels$models) { plotModelAccuracy(subm, graphs="data") }
    }
    if (any(grepl("sim", graphs))) {
        for (subm in submodels$models) { plotModelAccuracy(subm, graphs="simulation") }
    }
}

#' Compute the error of the model
#'
#' @param mra_model An MRAmodel object
#' @return A list with the simulation, the mismatch between the simulation and the data, and the residual of the fit
getModelError <- function(mra_model) {
    simulation = simulateModel(mra_model)
    mismatch = (mra_model$data$stim_data - simulation) / (mra_model$data$error*sqrt(2))
    residual = sum(mismatch^2, na.rm=T)
    return(list(simulation=simulation, mismatch=mismatch, residual=residual))
}

#' Plot the scores of each antibody
#'
#' Plot the scores of the fit for each antibody, which is how much
#' of the variation in the data is explained by the model
#' @param mra_model An MRAmodel object
#' @param ... Extra barplot parameters
#' @export
#' @family Model plots
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
plotModelScores.MRAmodel <- function(mra_model, ...) {
    low_lim=ifelse(min(mra_model$Rscores,na.rm=T)<0,min(-1,0.1*min(mra_model$Rscores,na.rm=T)),0)
    
    bb=barplot(mra_model$Rscores, xaxt="n", ylim=c(low_lim,1), las=1,ylab = "R^2", ...)
    lablist = names(mra_model$Rscores)
    text(bb, par("usr")[3] - 0.05, labels=lablist, srt=45, pos=1, xpd=TRUE)

    invisible(bb)
}

#' @rdname selectMinimalModel
reduceModel <- function(original_model, accuracy=0.95) {
    selectMinimalModel(original_model, accuracy)
}

#' Selection of a minimal model by iteratively removing insignificant links using the likelihood ratio test.
#' @param original_model An MRAmodel object, as the one produced by createModel or importModel
#' @param accuracy Probability threshold, the type I error for each link will be 1-accuracy. Multiple testing is not taken into account.
#' @return An MRAmodel object of the reduced model with the data
#' @export
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr} 
#' @author Bertram Klinger \email{bertram.klinger@@charite.de}
selectMinimalModel <- function(original_model, accuracy=0.95,verbose=F) {
  # Clone model object to not change original model specifications
  model_description = cloneModel(original_model)
  
  # Extra fitting informations from the model description
  model = model_description$model
  init_params = model_description$parameters
  expdes = model_description$design
  model_structure = model_description$structure
  adj = model_structure$adjacencyMatrix
  data = model_description$data
  
  if ("MRAmodelSet" %in% class(model_description)) {
    modelgroup=extractSubmodels(model_description)
    models = modelgroup$models
    n_par=model$nr_of_parameters()/model_description$nb_models
    initial_response = lapply(1:length(models), function(x) models[[x]]$model$getLocalResponseFromParameter( init_params[(1+(x-1)*n_par):(x*n_par)] ))
    rank = model$modelRank()/model_description$nb_models + length(model_description$variable_parameters)*(model_description$nb_models-1)
    if (length(model_description$variable_parameters)>0){
      variable_links = model$getParametersLinks()[model_description$variable_parameters]
    }else{
      variable_links=numeric()
    }
  } else {
    initial_response = model$getLocalResponseFromParameter( init_params )
    rank = model$modelRank()
  }
  
  if (is.na(model_description$bestfit)) {stop("A prior best fitting step is required to reduce the model from")}
  real_data = model_description$data$stim_data
  
  message("Performing model reduction...")
  init_residual = model_description$bestfit
  reduce=TRUE
  while (reduce) {
    links.to.test=which(adj==1)
    params=vector("list",length(links.to.test))
    residuals=c()
    ranks=c()
    newadj=adj
    c=0
  
    # Each link is removed and the best of those networks is compared to the previous model
    for (ii in links.to.test) {
      c=c+1
      newadj[ii]=0
      model_structure$setAdjacencyMatrix( newadj )
      model$setModel ( expdes, model_structure, model_description$use_log )
      
      if (class(model)== "Rcpp_ModelSet"){
        model$setNbModels(model_description$nb_models)
        paramstmp = unlist(lapply(1:length(models), function(x) model$getParameterFromLocalResponse(initial_response[[x]]$local_response, initial_response[[x]]$inhibitors)))
        # detect the parameter values that are fixed in the new network by having the exact same value in all following parameters
        leftVar = which( apply(matrix(paramstmp,ncol=model_description$nb_models,byrow=F), 1, not_duplicated) )
        model$setVariableParameters(leftVar) 
        result = model$fitmodelset(data, paramstmp, "levmar") 
        n_par=model$nr_of_parameters()/model_description$nb_models
        response.matrix = lapply(1:model_description$nb_models, function(x) model$getLocalResponseFromParameter( rep(result$parameter[(1+(x-1)*n_par):(x*n_par)], model_description$nb_models) ))
        params[[c]] = c(response.matrix)
        new_rank = model$modelRank()/model_description$nb_models + length(leftVar)*(model_description$nb_models-1)
      }else{
        paramstmp = model$getParameterFromLocalResponse(initial_response$local_response, initial_response$inhibitors)  
        result = model$fitmodel(data, paramstmp, "levmar")
        response.matrix = model$getLocalResponseFromParameter( result$parameter )
        params[[c]] = c(response.matrix)
        new_rank = model$modelRank()
      }
      
      residuals = c(residuals, result$residuals)   
      ranks = c(ranks, new_rank)
      
      if (verbose) {
        dr = rank - new_rank
        message(paste("old :", rank, ", new : ", new_rank))
        new_residual = residuals[length(residuals)]
        deltares = new_residual - init_residual
        message(paste(model_structure$names[(ii-1) %/% dim(adj)[1]+1], "->", model_structure$names[(ii-1) %% dim(adj)[1]+1], ": Delta residual = ", trim_num(deltares), "; Delta rank = ", dr, ", p-value = ", pchisq(deltares, df=dr) ))
      }
      
      newadj[ii]=1 ## Slightly accelerate the computation
    }
    
    order.res=order(residuals)
    # The loss of degree of freedom is equal to the difference in the ranks of the matrices
    new_rank = ranks[order.res[1]]
    new_residual = residuals[order.res[1]]
    dr = rank - new_rank
    if (dr==0) {
      warning(paste0("Link ",
                             model_structure$names[((links.to.test[order.res[1]]-1) %/% (dim(adj)[1])) +1], "->",
                             model_structure$names[((links.to.test[order.res[1]]-1) %% (dim(adj)[1])) +1], " belongs to a non-identifiable combination, setting df to 1."))
      dr=1
      }
    deltares = new_residual - init_residual
    chi_score =qchisq(accuracy,df=dr)
    if (deltares < 0) { warning(paste("Negative delta residual :", deltares)) ; deltares = -deltares  }
    if (deltares < chi_score) {
      adj[links.to.test[order.res[1]]]=0
      rank = new_rank
      initial_response=params[[order.res[1]]]
      #message(initial_response)
      message(paste0("Remove ",
                   model_structure$names[((links.to.test[order.res[1]]-1) %/% (dim(adj)[1])) +1], "->", # Line
                   model_structure$names[((links.to.test[order.res[1]]-1) %% (dim(adj)[1])) +1])); # Column (+1 because of the modulo and the R matrices starting by 1 instead of 0)
      
      message(paste( "New residual = ", residuals[order.res[1]], ", Delta residual = ", trim_num(deltares), ",  p-value = ", trim_num(pchisq(deltares, df=dr)) ))

      other_best = which((abs(residuals[order.res] - residuals[order.res[1]])) < 1e-4)[-1]
      if (length(other_best) > 0) {
          message("--- Other best links ---")
          for (lid in other_best) {
              deltares = residuals[order.res[lid]] - init_residual
              tmp_rank = ranks[order.res[lid]]
              tmp_dr = rank-tmp_rank
              if (tmp_dr==0) {
                warning(paste0("    Link ",
                               model_structure$names[((links.to.test[order.res[lid]]-1) %/% (dim(adj)[1])) +1], "->",
                               model_structure$names[((links.to.test[order.res[lid]]-1) %% (dim(adj)[1])) +1], " belongs to a non-identifiable combination, setting df to 1."))
                tmp_dr=1
              }
              message(paste0("    Could remove ",
                           model_structure$names[((links.to.test[order.res[lid]]-1) %/% (dim(adj)[1])) +1], "->", 
                           model_structure$names[((links.to.test[order.res[lid]]-1) %% (dim(adj)[1])) +1])); 
              message(paste( "    New residual = ", residuals[order.res[lid]], ", Delta residual = ", trim_num(deltares), ",  p-value = ", trim_num(pchisq(deltares, df=tmp_dr)) ))
          }
      }
      message("------------------------------------------------------------------------------------------------------")
      init_residual = residuals[order.res[1]]
      model_description$bestfit = sort(residuals)[1]
    } else {
      reduce=FALSE
    }
  }
  message("Reduction complete")
  # We recover the final model
  ## Basal activity and data do not change
  model_description$structure$setAdjacencyMatrix(adj)
  model_description$model$setModel(expdes, model_description$structure, original_model$use_log)
  if ("MRAmodelSet" %in% class(model_description)) {
    model_description$model$setNbModels(model_description$nb_models)
    model_description$parameters = unlist(lapply(1:length(models), function(x) model_description$model$getParameterFromLocalResponse(initial_response[[x]]$local_response, initial_response[[x]]$inhibitors)))
    # detect the parameter values that are fixed in the new network by having the exact same value in all following parameters
    leftVar = which( apply(matrix(model_description$parameters,ncol=model_description$nb_models,byrow=F), 1, not_duplicated) )
    model_description = setVariableParameters(model_description, leftVar) 
  } else {
    model_description$parameters = model_description$model$getParameterFromLocalResponse(initial_response$local_response, initial_response$inhibitors)
  }
  model_description$infos = c(model_description$infos, "Reduced model")
  model_description$param_range = list()
  model_description$lower_values = c()
  model_description$upper_values = c()

  if ("MRAmodelSet" %in% class(model_description)){
    return(model_description)
  } else{
  return(computeFitScore(model_description))
  }
}

#' Tries to locally add one link each and returns a list of links ordered by their chi-squared differences to the original model
#' 
#' A new link found to be suitable by the modeller can then added by re-running the createModel function with the altered adjacency information.
#' Note that values assigned to be exactly 1 inidicate almost always an non-identifiable link, whose combined value is assigned to another node in the combination!
#' @param original_model MRAmodel or MRAmodelSet object describing the model and its best fit, containing the data
#' @param parallel Boolean number indicating whether addition is executed in a parallel fashion
#' @param mc Number of cores that should be used for the computation
#' @param sample_range Numeric vector containing all starting values for the new link (DEFAULT: c(10^(2:-1),0,-10^(-1:2)))
#' @param padjust_method The method to use for the adjusted p-value, as defined in p.adjust. 'BY' by default which provides the FDR under general dependence assumption (conservative)
#' @param print Boolean indicating whether the result should be printed in a text file "Additional_link_suggestion.txt"
#' @param fname Name of the file where the results should be printed
#' @name suggestExtension
#' @export
#' @seealso selectMinimalModel, createModel
#' @author Bertram Klinger \email{bertram.klinger@@charite.de}
#' @examples \dontrun{
#' ext_list = suggestExtension(mramodel)
#' }
suggestExtension <- function(original_model,parallel = F, mc = 1, sample_range=c(10^(2:-1),0,-10^(-1:2)), padjust_method="bonferroni", print = F, fname="Additional_link_suggestion.txt"){
  if (mc == 0) {
      mc = detectCores() - 1
  }
  # Clone model object to not change original model specifications
  model_description = cloneModel(original_model)
  
  expdes = model_description$design
  model_structure = model_description$structure
  adj = model_structure$adjacencyMatrix
  data = model_description$data
  model = model_description$model
  
  if ("MRAmodelSet" %in% class(model_description)) {
    modelgroup=extractSubmodels(model_description)
    models = modelgroup$models
    n_par=model$nr_of_parameters()/model_description$nb_models
    initial_response = lapply(1:length(models), function(x) models[[x]]$model$getLocalResponseFromParameter( model_description$parameters[(1+(x-1)*n_par):(x*n_par)] ))
    rank = model$modelRank()/model_description$nb_models + length(model_description$variable_parameters)*(model_description$nb_models-1)
    variable_links = model_description$variable_parameters
    if (length(variable_links)>0){
      variable_links = model$getParametersLinks()[variable_links]
    }
  } else {
    initial_response = model$getLocalResponseFromParameter( model_description$parameters )
    rank = model$modelRank()
    variable_links = c()
  }

  if (is.na(model_description$bestfit)) {
    stop("A prior fit is required to reduce the model")
  }
  if (length(sample_range)==0){
    stop("'sample_range' should have at least one numeric value") 
  }
  
  message("Performing model extension...")
  init_residual = model_description$bestfit
  
  # determine the links that should be added, exclude self links and links acting on a stimulus (if not measured)
  exclude=diag(1,nrow(adj),ncol(adj))
  if (length(setdiff(expdes$stim_nodes,expdes$measured_nodes))>0){
    exclude[setdiff(expdes$stim_nodes,expdes$measured_nodes)+1,]=1
  }
  links_to_test=which( adj==0 & exclude==0)
  message(paste0(length(links_to_test)," links will be tested..."))
  
  # Each link is added and compared to the previous model
  if (parallel == T){
    extension_mat=mclapply(links_to_test,addLink,adj,rank,init_residual,model,initial_response,expdes,data,model_structure,sample_range,variable_links,mc.cores=mc)  
    extension_mat=as.data.frame(do.call("rbind",extension_mat))
  }else{
    cnames=c("adj_idx","from","to","value","residual","df","Res_delta","df_delta","pval")
    extension_mat=data.frame(matrix(NA,nrow=length(links_to_test),ncol=length(cnames),byrow=T))
    for (ii in 1:length(links_to_test))
      extension_mat[ii,]=addLink(links_to_test[ii],adj,rank,init_residual,model,initial_response,expdes,data,model_structure,sample_range,variable_links,verbose=T)
  colnames(extension_mat) <- cnames  
  } 
  extension_mat=extension_mat[order(as.numeric(as.matrix(extension_mat$Res_delta)),decreasing=T),]
  extension_mat=data.frame(extension_mat,"adj_pval"=p.adjust(as.numeric(as.matrix(extension_mat$pval)),method=padjust_method))
  
  message("Extension tests completed!")
  sig_res = sum(as.numeric(as.matrix(extension_mat$adj_pval))<=0.05, na.rm=T)
  if (sig_res > 0){
    select=match(c("from","to","value","Res_delta","adj_pval"),colnames(extension_mat))
    message(paste0(sig_res ," significant link extension",ifelse(sig_res>1,"s",""),"found"))
    sig_res=min(sig_res,10)
    message(paste0("printing the first ",sig_res," :\n"))
    message(paste(colnames(extension_mat)[select],collapse="\t"))
    for (ii in 1:sig_res){
    tmp = as.matrix(extension_mat[as.numeric(as.matrix(extension_mat$adj_pval))<=0.05,][ii,select])  
    tmp[,c("value","Res_delta","adj_pval")] = trim_num(tmp[,c("value","Res_delta","adj_pval")])  
    message(paste(tmp,collapse="\t"))
    }
  }
  if(print)
    write.table(extension_mat, fname, quote = F,row.names = F,sep="\t")
  
  return(extension_mat)
  
  #TODO so far only locally explores extension by assuming the starting values of all previously fitted parameters and a range of different starting values for the candidate link
}

#' add Link routine
#'
#' @param new_link Integer link whose addition is to be tested
#' @param adj integer Matrix original adjacency matrix excluding the new_link 
#' @param rank Rank of the input model
#' @param init_residual Numeric sum-squared error of original network
#' @param model MRAmodel or MRAmodelSet object of original network
#' @param initial_response List containing the local_response matrix and inhibitor strength of original network
#' @param expdes Design object of MRAmodel object 
#' @param data data Object of MRAmodel object
#' @param model_structure Structure object of MRAmodel object
#' @param sample_range Numeric vector containing all starting values for new_link
#' @param verbose Whether the function should be verbose or not
#' @author Bertram Klinger \email{bertram.klinger@@charite.de}
addLink <-  function(new_link,adj,rank,init_residual,model,initial_response,expdes,data,model_structure,sample_range,variable_links=c(),verbose=F){
  adj[new_link] = 1
  model_structure$setAdjacencyMatrix( adj )
  model$setModel ( expdes, model_structure, model$use_log )
  if (class(model) == "Rcpp_ModelSet"){
    model$setNbModels(length(initial_response))
    if (length(variable_links)>0) {  
      for (ii in 1:length(initial_response)){ initial_response[[ii]]$local_response[new_link]=1} # set new link to one to preserve constant verus variable links
      paramstmp = unlist(lapply(1:length(initial_response), function(x) model$getParameterFromLocalResponse(initial_response[[x]]$local_response, initial_response[[x]]$inhibitors)))
      leftVar = which( apply(matrix(paramstmp,ncol=length(initial_response),byrow=F), 1, not_duplicated) )
      model$setVariableParameters(leftVar)
    }else{
        leftVar=numeric()
        }
    }
  best_res = Inf
  for (jj in sample_range){
    if (class(model)== "Rcpp_ModelSet"){
      paramstmp = c()
      for (ii in 1:length(initial_response)){
      initial_response[[ii]]$local_response[new_link]=jj
      paramstmp = c(paramstmp, model$getParameterFromLocalResponse(initial_response[[ii]]$local_response, initial_response[[ii]]$inhibitors)) 
      }
      tmp_result = model$fitmodelset( data, paramstmp, "levmar")  
    }else{
    initial_response$local_response[new_link]=jj
    paramstmp = model$getParameterFromLocalResponse(initial_response$local_response, initial_response$inhibitors)  
    tmp_result = model$fitmodel( data, paramstmp, "levmar" )
    }
    if ( verbose == T )
    message( paste0( "for ", jj ," : ",tmp_result$residuals ) )
    if ( tmp_result$residuals < best_res ){
      best_res = tmp_result$residuals
      result = tmp_result
    }
  }
  response_matrix = model$getLocalResponseFromParameter( result$parameter )
  if (class(model) == "Rcpp_ModelSet"){
    new_rank = model$modelRank()/length(initial_response)+length(leftVar)*(length(initial_response)-1)
  }else{
  new_rank = model$modelRank()
  }
  dr = new_rank-rank
  deltares = init_residual-result$residuals
  extension_mat = matrix(c(new_link,
                           model_structure$names[(new_link-1) %/% dim(adj)[1]+1],
                           model_structure$names[(new_link-1) %% dim(adj)[1]+1],
                           response_matrix$local_response[new_link],
                           result$residuals,
                           new_rank,
                           deltares,
                           dr,
                           ifelse(dr > 0, 1-pchisq(deltares, df=dr), 1) ),nrow=1 )
  colnames(extension_mat) <- c("adj_idx","from","to","value","residual","df","Res_delta","df_delta","pval")
  adj[new_link] = 0
  model_structure$setAdjacencyMatrix( adj )
  model$setModel ( expdes, model_structure, model$use_log )
  if (class(model)== "Rcpp_ModelSet"){
    model$setNbModels(length(initial_response))
    if (length(variable_links)>0) { model$setVariableParameters(match(variable_links,model$getParametersNames()$names))}
    for (ii in 1:length(initial_response)){
      initial_response[[ii]]$local_response[new_link]=0
    }
  }else{
    initial_response$local_response[new_link]=0
  }
    if ( verbose == T ){
    message(paste("[",extension_mat[1], "]" ,
                     ", new : ", new_rank,
                     extension_mat[2],"->",
                     extension_mat[3],
                     ": Delta residual = ",
                     trim_num(extension_mat[7]),
                     "; Delta rank = ",
                     extension_mat[8],
                     ", p-value = ",
                     trim_num(extension_mat[9]) ))  
  }
  return(extension_mat)
}

#' Computes the fitting scores for a new parameter set
#'
#' Test the model with the provided parameter set and returns the fit, the scores and the (possibly updated) parameter set
#' @inheritParams computeFitScore
#' @param new_parameters A vector of parameters to use for the new fit
#' @return An objet of class MRAmodel
#' @export
# TODO specialisation of update ??
testModel <- function(mra_model, new_parameters, refit_model=FALSE) {
    if (length(new_parameters) != length(mra_model$parameters)) {
        stop("The number of parameters is incorrect")
    }
    tmp_model = mra_model
    tmp_model$parameters = new_parameters
    tmp_model = computeFitScore(tmp_model, refit_model)
    tmp_model$bestfit = getModelError(tmp_model)$residual

    return(tmp_model)
}

#' Refit the model
#'
#' Refit the model with a specified parameter set while keeping parameters constant
#' @param mra_model A MRAmodel object
#' @param parameter_set A vector of values used as parameters for the model. There must be as many values as there are parameters, or one that will be used for all parameters.
#' @param vary_param A vector of index or name of the parameters to refit, the others will be kept constant. Repetitions or redundant information (index and name designating the same parameter) are removed.
#' @param inits Number of random initialisations for the variable parameters
#' @param nb_cores Number of processes to use for the refitting. 0 to use all cores of the machines but one.
#' @param method Method to use for the sample generation for the random initialisations
#' @param fit_name Name of the refit for the title of the plots
#' @return The refitted model as an MRAmodel object.
#' @seealso printParametersNames
#' @name refitModel
#' @export
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr} 
refitModel <- function(mra_model, parameter_set=c(), vary_param=c(), inits=100, nb_cores=1, method="randomlhs", fit_name="") {
    if (length(parameter_set) == 0) {
        stop("No 'parameter_set' provided in 'refitModel'")
    } else if (length(parameter_set) == 1) {
        parameter_set = rep(parameter_set, length(mra_model$parameters))
    } else if (length(parameter_set) != length(mra_model$parameters)) {
        stop("Incompatible 'parameter_set', wrong number of parameters in 'refitModel'")
    }
    if (length(vary_param) == 0) {
        stop("No parameters to vary ('vary_param==c()') in 'refitModel', use 'computeFitScore' if you just want to recompute the model scores")
    }
    if (nb_cores == 0) {
        nb_cores = detectCores()-1
    }

    pnames = getParametersNames(mra_model)
    keep_constant = 1:length(mra_model$parameters)
    for (ii in 1:length(vary_param)) {
        if (!is.numeric(suppressWarnings(as.numeric(vary_param[ii])))) {
            if (!vary_param[ii] %in% pnames) {
                stop(paste0(vary_param[ii], " is not a valid parameter name"))
            }
            vary_param[ii] = which(pnames==vary_param[ii])
        } else if (!vary_param[ii] %in% keep_constant) {
            stop(paste0(vary_param[ii], " is not a valid index for the model"))
        }
    }
    vary_param = unique(as.numeric(vary_param))
    samples = getSamples(length(vary_param), inits-1, method, nb_cores)
    init_pset = matrix(rep(parameter_set, inits), byrow=TRUE, nrow=inits)
    for (ii in 1:length(vary_param)) {
        init_pset[2:inits,vary_param[ii]] = samples[,ii]
        keep_constant = keep_constant[-which( keep_constant==vary_param[ii] )]
    }
    results = parallel_initialisation(mra_model$model, mra_model$data, init_pset, nb_cores, keep_constant)
    order_id = order(results$residuals)
    residuals_plot(results$residuals, fit_name)

    new_model = cloneModel(mra_model) 
    new_model$parameters = results$params[order_id[1],]
    new_model$bestfit = results$residuals[order_id[1]]
    new_model$infos = c(new_model$infos, paste0("Refitted with variable parameters c(", pastecoma(vary_param), ")") )
    if (fit_name != "") { new_model$name = fit_name }
    return( computeFitScore(new_model, FALSE) )
}

#' Fit a model using the parameter set from another model
#' @rdname refit
#' @export
fitFromModel <- function(mra_model, parameters_model, vary_param=c(), inits=100, nb_cores=1, method="randomlhs", fit_name="") {
# Add controls
    plp = parameters_model$model$getLocalResponseFromParameter(parameters_model$parameters)
    new_pset = mra_model$model$getParameterFromLocalResponse(plp$local_response, plp$inhibitors)

    return( refitModel(mra_model, new_pset, vary_param, inits, nb_cores, method, fit_name) )
}
