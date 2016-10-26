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
#' @return Nothing
#' @export
#' @seealso createModel, importModel
#' @family Model plots
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
plotModelAccuracy <- function(model_description) {
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

  message("Treatments : ")
  message(paste(treatments, collapse=" "))
  colnames(mismatch) = colnames(stim_data) = colnames(simulation) = nodes[design$measured_nodes + 1]
  rownames(mismatch) = rownames(stim_data) = rownames(simulation) = treatments

# Comparison of the data and the stimulation in term of error fold change and log fold change
  plotHeatmap(mismatch,"(data - simulation) / error")
  plotHeatmap(stim_data-simulation,"log2(data/simulation)")
# Log fold changes for the data and the stimulation with comparable color code
  lim=min(10, max(abs( range(quantile(stim_data,0.05, na.rm=T),
                       quantile(simulation,0.05, na.rm=T),
                       quantile(stim_data,0.95, na.rm=T),
                       quantile(simulation,0.95, na.rm=T)) )))
  plotHeatmap(stim_data, "Log-fold change Experimental data",lim,T)
  plotHeatmap(simulation, "Log-fold change Simulated data",lim,T)

  invisible(list(mismatch=mismatch, stim_data=stim_data, simulation=simulation))
}

#' Compute the error of the model
#'
#' @param mra_model An MRAmodel object
#' @return A list with the simulation, the mismatch between the simulation and the data, and the residual of the fit
getModelError <- function(mra_model) {
    simulation = mra_model$model$simulate(mra_model$data, mra_model$parameters)$prediction
    mismatch = (mra_model$data$stim_data - simulation) / mra_model$data$error
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
    bb=barplot(mra_model$Rscores, xaxt="n", ...)
    lablist = names(mra_model$Rscores)
    text(bb, par("usr")[3] - 0.05, labels=lablist, srt=45, pos=1, xpd=TRUE)

    invisible(bb)
}

#' @rdname selectMinimalModel
reduceModel <- function(model_description, accuracy=0.95) {
    selectMinimalModel(model_description, accuracy)
}

#' Selection of a minimal model by the removal of non significant links with an F-test
#' @param model_description An MRAmodel object, as the one produced by createModel or importModel
#' @param accuracy Probability threshold, the type I error for each link will be 1-accuracy. Multiple testing is not taken into account.
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
  real_data = model_description$data$stim_data
  data_count = sum(!is.na(real_data) & !is.nan(real_data))
  
  message("Performing model reduction...")
  init_residual = model_description$bestfit
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
      result = model$fitmodel(data, paramstmp)
      response.matrix = model$getLocalResponseFromParameter( result$parameter )
      residuals = c(residuals, result$residuals);   
      params = cbind(params, c(response.matrix))
      new_rank = model$modelRank()
      ranks = c(ranks, new_rank)
      
      if (verbose) {
        dr = rank - new_rank
        message(paste("old :", rank, ", new : ", new_rank))
        new_residual = residuals[length(residuals)]
        dfreedom = data_count - rank
        deltares = new_residual - init_residual
        f_score = ((new_residual - init_residual) / dr) / (init_residual/dfreedom)
        message(paste(model_structure$names[(i-1) %/% dim(adj)[1]+1], "->", model_structure$names[(i-1) %% dim(adj)[1]+1], ": Delta residual = ", deltares, "; Delta rank = ", dr, ", p-value = ", pf(f_score, dr, dfreedom) ))
      }
      
      newadj[i]=1; ## Slightly accelerate the computation
    }
    
    order.res=order(residuals)
    # The loss of degree of freedom is equal to the difference in the ranks of the matrices
    new_rank = ranks[order.res[1]]
    new_residual = residuals[order.res[1]]
    dr = rank - new_rank
    deltares = new_residual - init_residual
    dfreedom = data_count - rank
    f_score = ((new_residual - init_residual) / dr) / (init_residual/dfreedom)
    # Some boundary cases might give low improvement of the fit
    if (deltares < 0) { warning(paste("Negative delta residual :", deltares)) ; deltares = -deltares  }
    if (f_score < qf(accuracy, df1=dr, df2=dfreedom)) {
      adj[links.to.test[order.res[1]]]=0
      rank = new_rank
      initial_response=params[,order.res[1]]
      init_residual = residuals[order.res[1]]
      #message(initial_response)
      message(paste0("Remove ",
                   model_structure$names[((links.to.test[order.res[1]]-1) %/% (dim(adj)[1])) +1], "->", # Line
                   model_structure$names[((links.to.test[order.res[1]]-1) %% (dim(adj)[1])) +1])); # Column (+1 because of the modulo and the R matrices starting by 1 instead of 0)
      
      message(paste( "New residual = ", residuals[order.res[1]], ", Delta residual = ", deltares, ",  p-value = ", pf(f_score, df1=dr, df2=dfreedom) ))
      message("------------------------------------------------------------------------------------------------------")
      
      model_description$bestfit = sort(residuals)[1]
    } else {
      reduce=FALSE
    }
  }
  message("Reduction complete")
  # We recover the final model
  ## Basal activity and data do not change
  model_description$structure$setAdjacencyMatrix(adj)
  model_description$model$setModel(expdes, model_description$structure)
  model_description$parameters = model_description$model$getParameterFromLocalResponse(initial_response$local_response, initial_response$inhibitors)
  model_description$infos = c(model_description$infos, "Reduced model")
  model_description$param_range = list()
  model_description$lower_values = c()
  model_description$upper_values = c()
# TODO either create another model or update the statistics
  
  return(computeFitScore(model_description))
}

#' Tries to add one link each and returns a list of links ordered by their chi-squared differences to the original model
#' This list can then be based and compared to literature knowledge and if considered suitable manually added to the starting network and rerun with the createModel function.
#' @param model_description MRAmodel object describing thze model and its best fit, containing the data
#' @param parallel Boolean number indicating whether addition is executed in a parallel fashion
#' @param mc Number of cores that should be used for the computation
#' @param sample_range Numeric vector containing all starting values for the new link (DEFAULT: c(10^(2:-1),0,-10^(-1:2)))
#' @param print Boolean indicating whether the result should be printed in a text file "Additional_link_suggestion.txt"
#' @export
#' @seealso selectMinimalModel, createModel
#' @author Bertram Klinger \email{bertram.klinger@@charite.de}
#' @examples \dontrun{
#' ext_list = suggestExtension(mramodel)
#' }
suggestExtension <- function(model_description,parallel = F, mc = 1, sample_range=c(10^(2:-1),0,-10^(-1:2)), print = F){
  # Extra fitting informations from the model description
  model = model_description$model
  initial_response = model$getLocalResponseFromParameter( model_description$parameters )
  expdes = model_description$design
  model_structure = model_description$structure
  adj = model_structure$adjacencyMatrix
  data = model_description$data
  if (is.na(model_description$bestfit)) {
    stop("A prior fit is required to reduce the model")
  }
  if (length(sample_range)==0){
    stop("'sample_range' should have at least one numeric value") 
  }
  
  message("Performing model extension...")
  init_residual = model_description$bestfit
  rank = model$modelRank()
  #determine the links that should be added, exclude self links and links acting on a stimulus (if not measured)
  exclude=diag(1,nrow(adj),ncol(adj))
  if (length(setdiff(expdes$stim_nodes,expdes$measured_nodes))>0){
    exclude[setdiff(expdes$stim_nodes,expdes$measured_nodes)+1,]=1
  }
  links_to_test=which( adj==0 & exclude==0)
  message(paste0(length(links_to_test)," links will be tested..."))
  
  # Each link is added and compared to the previous model
  if (parallel == T){
    extension_mat=mclapply(links_to_test,addLink,adj,rank,init_residual,model,initial_response,expdes,data,model_structure,sample_range,mc.cores=mc)  
    extension_mat=as.data.frame(do.call("rbind",extension_mat))
  }else{
    cnames=c("adj_idx","from","to","value","residual","df","Res_delta","df_delta","pval")
    extension_mat=data.frame(matrix(NA,nrow=length(links_to_test),ncol=length(cnames),byrow=T))
    for (ii in 1:length(links_to_test))
      extension_mat[ii,]=addLink(links_to_test[ii],adj,rank,init_residual,model,initial_response,expdes,data,model_structure,sample_range,verbose=T)
  colnames(extension_mat) <- cnames  
  } 
  extension_mat=extension_mat[order(as.numeric(as.matrix(extension_mat$Res_delta)),decreasing=T),]
  extension_mat=data.frame(extension_mat,"adj_pval"=p.adjust(as.numeric(as.matrix(extension_mat$pval)),method="BH"))
  
  message("Extension tests completed!")
  sig_res = sum(as.numeric(as.matrix(extension_mat$adj_pval))<=0.05)
  if (sig_res > 0){
    select=match(c("from","to","value","Res_delta","adj_pval"),colnames(extension_mat))
    message(paste0(sig_res ," significant link extensions found"))
    sig_res=min(sig_res,10)
    message(paste0("printing the first ",sig_res," :\n"))
    message(paste(colnames(extension_mat)[select],collapse="\t"))
    for (ii in 1:sig_res){
    message(paste(as.matrix(extension_mat[as.numeric(as.matrix(extension_mat$adj_pval))<=0.05,][ii,select]),collapse="\t"))
    }
  }
  if(print)
    write.table(extension_mat,"Additional_link_suggestion.txt",quote = F,row.names = F,sep="\t")
  
  return(extension_mat)
  
  #TODO so far only locally explores extension by assuming the starting values of all previously fitted parameters and a range of different starting values for the candidate link
}

#' add Link routine
#'
#' @param new_link Integer link whose addition is to be tested
#' @param adj integer Matrix original adjacency matrix excluding the new_link 
#' @param rank Rank of the input model
#' @param init_residual Numeric sum-squared error of original network
#' @param model MRAmodel object of original network
#' @param initial_response List containing the local_response matrix and inhibitor strength of original network
#' @param expdes Design object of MRAmodel object 
#' @param data data Object of MRAmodel object
#' @param model_structure Structure object of MRAmodel object
#' @param sample_range Numeric vector containing all starting values for new_link
#' @param verbose Whether the function should be verbose or not
addLink <-  function(new_link,adj,rank,init_residual,model,initial_response,expdes,data,model_structure,sample_range,verbose=F){
  adj[new_link] = 1
  model_structure$setAdjacencyMatrix( adj )
  model$setModel ( expdes, model_structure )
  best_res = Inf
  for (jj in sample_range){
    initial_response$local_response[new_link]=jj
    paramstmp = model$getParameterFromLocalResponse(initial_response$local_response, initial_response$inhibitors)
    tmp_result = model$fitmodel( data,paramstmp )
    if ( verbose == T )
    message( paste0( "for ", jj ," : ",tmp_result$residuals ) )
    if ( tmp_result$residuals < best_res ){
      best_res = tmp_result$residuals
      result = tmp_result
    }
  }
  response_matrix = model$getLocalResponseFromParameter( result$parameter )
  new_rank = model$modelRank()
  dr = new_rank-rank
  deltares = init_residual-result$residuals
  data_count = sum(!is.na(data$stim_data) & !is.nan(data$stim_data))
  dfreedom = data_count - length(result$parameters)
  f_score = (deltares/dr) / (result$residuals/dfreedom)
  extension_mat = matrix(c(new_link,
                           model_structure$names[(new_link-1) %/% dim(adj)[1]+1],
                           model_structure$names[(new_link-1) %% dim(adj)[1]+1],
                           response_matrix$local_response[new_link],
                           result$residuals,
                           new_rank,
                           deltares,
                           dr,
                           1-pf(f_score, dr, dfreedom)),nrow=1)  
  colnames(extension_mat) <- c("adj_idx","from","to","value","residual","df","Res_delta","df_delta","pval")
  adj[new_link] = 0
  model_structure$setAdjacencyMatrix( adj )
  model$setModel ( expdes, model_structure )
    if ( verbose == T ){
    message(paste("[",extension_mat[1], "]" ,
                     ", new : ", new_rank,
                     extension_mat[2],"->",
                     extension_mat[3],
                     ": Delta residual = ",
                     ifelse(as.numeric(extension_mat[7])>1,round(as.numeric(extension_mat[7]),2),signif(as.numeric(extension_mat[7]),2)),
                     "; Delta rank = ",
                     extension_mat[8],
                     ", p-value = ",
                     signif(as.numeric(extension_mat[9]),2) ))  
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

