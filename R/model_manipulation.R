########################### model_manipulation.R ###########################
# Functions to change and visualise the model

#' Print the value of each path from the model, and add the profile likelihood infos if they are provided
#' @param model_description An MRAmodel object
#' @return Nothing
#' @export
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
printParameters <- function(model_description, precision=2) {
  model = model_description$model
  parameters = model_description$parameters
  
  print("Parameters :")
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
      print(text)
    }
  } else {
    for (i in 1:length(paths)) {
      print (paste( simplify_path_name(paths[i]), ":", signif(parameters[i], precision) ))
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
  mismatch = (stim_data - simulation) / error
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
  print("Treatments : ")
  print(treatments)
  colnames(mismatch) = colnames(stim_data) = colnames(simulation) = nodes[design$measured_nodes + 1]
  rownames(mismatch) = rownames(stim_data) = rownames(simulation) = treatments
  
  # Comparison of the data and the stimulation in term of error fold change and log fold change
  plot_heatmap(mismatch,"(data - simulation) / error")
  plot_heatmap(stim_data-simulation,"log2(data/simulation)")
  # Log fold changes for the data and the stimulation with comparable color code
  lim=min(10,abs(range(quantile(stim_data,0.05, na.rm=T),
                       quantile(simulation,0.05, na.rm=T),
                       quantile(stim_data,0.95, na.rm=T),
                       quantile(simulation,0.95, na.rm=T))))
  plot_heatmap(stim_data, "Log-fold change Experimental data",lim,T)
  plot_heatmap(simulation, "Log-fold change Simulated data",lim,T)
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

#' Tries to add one link each and returns a list of ordered chi-squared differences, highlighting the significant ones
#' @param model_description MRAmodel object describing thze model and its best fit, containing the data
#' @param nb_cores Number of cores that should be used for the computation
#' @param inits Number of initialisation steps which should be performed (see method for the exact meaning of this value)
#' @param method Method to be used for the initialisation, available methods are :
#'      random : Perform a Latin Hypercube Sampling to choose \emph{inits} starting points then perform a gradient descent to find the local optimum for each of those points.
#'      correlation : Deduce some parameters from the correlation between the measurements for the target node and all of its input nodes, then perform random to find the other parameters. Recommended, very efficient for small datasets.
#'      genetic : Genetic algorithm with mutation only. \emph{inits} is the total number of points sampled.
#'      annealing : Simulated annealing and gradient descent on the best result. \emph{inits} is the maximum number of iteration without any change before the algorithm decides it reached the best value. Use not recommended.
#' @return An MRAmodel object describing the model and its best fit, containing the data
#' @export
#' @seealso createModel, initModel
#' @author Bertram Klinger \email{bertram.klinger@@charite.de}
#' @examples
#' ext_list = createModel(model) # Produces a model for the network described in links.tab using the data in data_MIDES.csv
suggestExtension <- function(model_description,nb_cores=1,inits=1000,method="geneticlhs",print=F){
  
  # Extra fitting informations from the model description
  model = model_description$model
  init_params = model_description$parameters
  initial_response = model$getLocalResponseFromParameter( init_params )
  expdes = model_description$design
  model_structure = model_description$structure
  adj = model_structure$adjacencyMatrix
  data = model_description$data
  if (is.na(model_description$bestfit)) {stop("Data are required to reduce the model")}
  
  writeLines("Performing model extension…")
  initresidual = model_description$bestfit
  rank = model$modelRank()
  links.to.test=which(adj==0 & diag(1,nrow(adj),ncol(adj))==0)
  writeLines(paste0(length(links.to.test)," links will be tested..."))
  tmp_adj=adj
  
  # Each link is added and compared to the previous model
  extension_mat=NULL
  sample=c(-1,0,1)
  for (i in links.to.test) {
    tmp_adj[i]=1
    model_structure$setAdjacencyMatrix( tmp_adj )  
    model$setModel ( expdes, model_structure )
    best_res=Inf
    for (j in sample){
      initial_response$local_response[i]=j
      paramstmp = model$getParameterFromLocalResponse(initial_response$local_response, initial_response$inhibitors)
      tmp_result = model$fitmodel(data,paramstmp)
      writeLines(paste0("for ", j ," : ",tmp_result$residuals))
      if (tmp_result$residuals<best_res){
        best_res=tmp_result$residuals
        result=tmp_result
      }
    }
    response_matrix = model$getLocalResponseFromParameter( result$parameter )
    new_rank = model$modelRank()
    dr = new_rank-rank
    deltares = initresidual-result$residuals
    extension_mat=rbind(extension_mat,c(i,
                                        model_structure$names[(i-1) %/% dim(adj)[1]+1],
                                        model_structure$names[(i-1) %% dim(adj)[1]+1],
                                        response_matrix$local_response[i],
                                        result$residuals,
                                        new_rank,
                                        deltares,
                                        dr,
                                        1-pchisq(deltares, df=dr)))  
    writeLines(paste("[",which(links.to.test==i), "]" ,
                     " old :", rank, 
                     ", new : ", new_rank,
                     extension_mat[nrow(extension_mat),2],"->",
                extension_mat[nrow(extension_mat),3],
                ": Delta residual = ",
                extension_mat[nrow(extension_mat),7],
                "; Delta rank = ",
                extension_mat[nrow(extension_mat),8],
                ", p-value = ",
                extension_mat[nrow(extension_mat),9] ))
    
    tmp_adj[i]=0; ## Slightly accelerates the computation
    initial_response$local_response[i]=0
    model_structure$setAdjacencyMatrix( tmp_adj )
    model$setModel ( expdes, model_structure )
  }
  colnames(extension_mat) <- c("adj_idx","from","to","value","residual","df","Res_delta","df_delta","pval")
  
  extension_mat=extension_mat[order(as.numeric(as.matrix(extension_mat$Res_delta)),decreasing=T),]
  extension_mat=data.frame(extension_mat,"adj_pval"=p.adjust(extension_mat$pval,method="BH"))
  
  # convert numeric parts into numeric
  
  writeLines("Extension trial completed!")
  if (any(extension_mat$adj_pval<=0.05)){
  print("Significant link extensions:")
  print(extension_mat[extension_mat$adj_pval<=0.05,])
  }else{
  print("no significant extension could be found!")  
  }
  if(print)
    write.table(extension_mat,"Additional_link_suggestion.txt",quote = F,row.names = F)
  
  return(extension_mat)
  
  #TODO So far only extension is started from previously best parameterisation and -1 0 and 1 as starting value for the new parameter-> local extension   
  # only returns suggestions for the improvment if the model does not give any insights into whatever
}

