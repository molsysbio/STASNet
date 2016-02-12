########################### model_manipulation.R ###########################
# Functions to change and visualise the model

#' Print the value of each path from the model, with the profile likelihood infos if they are provided
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

#' Plot heatmaps of the model simulation against the data weighted by the error, as well as the log fold change for the data and the prediction
#' @param model_description An MRAmodel object
#' @return Nothing
#' @export
#' @seealso plotModelSimulation, createModel, importModel
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

print("Treatments : ")
print(treatments)
colnames(mismatch) = colnames(stim_data) = colnames(simulation) = nodes[design$measured_nodes + 1]
rownames(mismatch) = rownames(stim_data) = rownames(simulation) = treatments

# Comparison of the data and the stimulation in term of error fold change and log fold change
plot_heatmap(mismatch,"(data - simulation) / error")
plot_heatmap(stim_data-simulation,"log2(data/simulation)")
# Log fold changes for the data and the stimulation with comparable color code
lim=min(10, max( range(quantile(stim_data,0.05, na.rm=T),
                     quantile(simulation,0.05, na.rm=T),
                     quantile(stim_data,0.95, na.rm=T),
                     quantile(simulation,0.95, na.rm=T)) ))
plot_heatmap(stim_data, "Log-fold change Experimental data",lim,T)
plot_heatmap(simulation, "Log-fold change Simulated data",lim,T)
}

#' Plot the scores of each antibody
#'
#' Plot the scores of the fit for each antibody, which is how much
#' of the variation in the data is explained by the model
#' @param mra_model An MRAmodel object
#' @export
#' @author Mathurin Dorel \email{dorel@@horus.ens.fr}
plotModelScores.MRAmodel <- function(mra_model, ...) {
    bb=barplot(mra_model$Rscores, xaxt="n", ...)
    lablist = names(mra_model$Rscores)
    text(bb, par("usr")[3] - 0.05, labels=lablist, srt=45, pos=1, xpd=TRUE)

    invisible(bb)
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
      result = model$fitmodel( data,paramstmp)
      response.matrix = model$getLocalResponseFromParameter( result$parameter )
      residuals = c(residuals, result$residuals);   
      params = cbind(params, c(response.matrix))
      new_rank = model$modelRank()
      ranks = c(ranks, new_rank)
      
      if (verbose) {
        dr = rank - new_rank
        print(paste("old :", rank, ", new : ", new_rank))
        deltares = residuals[length(residuals)]-init_residual
        print(paste(model_structure$names[(i-1) %/% dim(adj)[1]+1], "->", model_structure$names[(i-1) %% dim(adj)[1]+1], ": Delta residual = ", deltares, "; Delta rank = ", dr, ", p-value = ", pchisq(deltares, df=dr) ))
      }
      
      newadj[i]=1; ## Slightly accelerate the computation
    }
    
    order.res=order(residuals)
    # The loss of degree of freedom is equal to the difference in the ranks of the matrices
    new_rank = ranks[order.res[1]]
    dr = rank - new_rank
    deltares = residuals[order.res[1]]-init_residual
    # Some boundary cases might give low improvement of the fit
    if (deltares < 0) { warning(paste("Negative delta residual :", deltares)) ; deltares = -deltares  }
    if (deltares < qchisq(accuracy, df=dr)) {
      adj[links.to.test[order.res[1]]]=0
      rank = new_rank
      initial_response=params[,order.res[1]]
      init_residual = residuals[order.res[1]]
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

#' Tries to add one link each and returns a list of links ordered by their chi-squared differences to the original model
#' This list can then be based and compared to literature knowledge and if considered suitable manually added to the starting network and rerun with the createModel function.
#' @param model_description MRAmodel object describing thze model and its best fit, containing the data
#' @param parallel Boolean number indicating whether addition is executed in a parallel fashion
#' @param mc Number of cores that should be used for the computation
#' @param print Boolean indicating whether the result should be printed in a text file "Additional_link_suggestion.txt"
#' @param inits DEFUNCT Number of initialisation steps which should be performed (see method for the exact meaning of this value)
#' @param method DEFUNCT Method to be used for the initialisation, available methods are :
#'      random : Perform a Latin Hypercube Sampling to choose \emph{inits} starting points then perform a gradient descent to find the local optimum for each of those points.
#'      correlation : Deduce some parameters from the correlation between the measurements for the target node and all of its input nodes, then perform random to find the other parameters. Recommended, very efficient for small datasets.
#'      genetic : Genetic algorithm with mutation only. \emph{inits} is the total number of points sampled.
#'      annealing : Simulated annealing and gradient descent on the best result. \emph{inits} is the maximum number of iteration without any change before the algorithm decides it reached the best value. Use not recommended.
#' @return a data.frame sorting the links from most significant to least significant
#' @export
#' @seealso selectMinimalModel, createModel
#' @author Bertram Klinger \email{bertram.klinger@@charite.de}
#' @examples
#' ext_list = suggestExtension(MRAmodel) 
suggestExtension <- function(model_description,parallel = F,mc = 1,print = F,inits = 1000,method = "geneticlhs"){
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
  init_residual = model_description$bestfit
  rank = model$modelRank()
  links_to_test=which( adj==0 & diag(1,nrow(adj),ncol(adj))==0)
  writeLines(paste0(length(links_to_test)," links will be tested..."))
  sample = c(10^c(2:-6),0,-10^c(2:-6)) 
  
  # Each link is added and compared to the previous model
  if (parallel == T){
    extension_mat=mclapply(links_to_test,addLink,adj,rank,init_residual,model,initial_response,expdes,data,model_structure,sample,mc.cores=mc)  
    extension_mat=as.data.frame(do.call("rbind",extension_mat))
  }else{
    tmp_adj=adj
    extension_mat=NULL
    for (i in links_to_test) {
      tmp_adj[i] = 1
      model_structure$setAdjacencyMatrix( tmp_adj )
      model$setModel ( expdes, model_structure )
      best_res = Inf
      for (j in sample){
        initial_response$local_response[i]=j
        paramstmp = model$getParameterFromLocalResponse(initial_response$local_response, initial_response$inhibitors)
        tmp_result = model$fitmodel(data,paramstmp)
        writeLines( paste0( "for ", j ," : ",tmp_result$residuals ) )
        if ( tmp_result$residuals < best_res ){
          best_res = tmp_result$residuals
          result = tmp_result
        }
      }
      response_matrix = model$getLocalResponseFromParameter( result$parameter )
      new_rank = model$modelRank()
      dr = new_rank-rank
      deltares = init_residual-result$residuals
      extension_mat = rbind(extension_mat,c(i,
                                            model_structure$names[(i-1) %/% dim(adj)[1]+1],
                                            model_structure$names[(i-1) %% dim(adj)[1]+1],
                                            response_matrix$local_response[i],
                                            result$residuals,
                                            new_rank,
                                            deltares,
                                            dr,
                                            1-pchisq(deltares, df=dr)))  
      writeLines(paste("[",which(links_to_test == i), "]" ,
                       ", new : ", new_rank,
                       extension_mat[nrow(extension_mat),2],"->",
                       extension_mat[nrow(extension_mat),3],
                       ": Delta residual = ",
                       extension_mat[nrow(extension_mat),7],
                       "; Delta rank = ",
                       extension_mat[nrow(extension_mat),8],
                       ", p-value = ",
                       extension_mat[nrow(extension_mat),9] ))
      
      tmp_adj[i] = 0
      model_structure$setAdjacencyMatrix( tmp_adj )
      model$setModel ( expdes, model_structure )
      
    }
    colnames(extension_mat) <- c("adj_idx","from","to","value","residual","df","Res_delta","df_delta","pval")
    extension_mat=data.frame(extension_mat)
  }
  extension_mat=extension_mat[order(as.numeric(as.matrix(extension_mat$Res_delta)),decreasing=T),]
  extension_mat=data.frame(extension_mat,"adj_pval"=p.adjust(as.numeric(as.matrix(extension_mat$pval)),method="BH"))
  
  writeLines("Extension trial completed!")
  if (length(as.numeric(as.matrix(extension_mat$adj_pval))<=0.05)>0){
    print("Significant link extensions:")
    print(extension_mat[as.numeric(as.matrix(extension_mat$adj_pval))<=0.05,])
  }
  if(print)
    write.table(extension_mat,"Additional_link_suggestion.txt",quote = F,row.names = F)
  
  return(extension_mat)
  
  #TODO so far only locally explores extension by assuming the starting values of all previously fitted parameters and a starting value of the new parameter of either -1,0, or 1
}

#' add Link routine
#'
#' @param new_link integer link whose addition is to be tested
#' @param adj integer matrix original adjacency matrix excluding the new_link 
#' @param init_residual numeric sum-squared error of original network
#' @param model MRAmodel object of original network
#' @param initial_response list containing the local_response matrix and inhibitor strength of original network
#' @param expdes design object of MRAmodel object 
#' @param data data object of MRAmodel object
#' @param model_structure structure object of MRAmodel object
#' @param sample numeric vector containing all starting values for new_link
addLink <-  function(new_link,adj,rank,init_residual,model,initial_response,expdes,data,model_structure,sample){
  adj[new_link] = 1
  model_structure$setAdjacencyMatrix( adj )
  model$setModel ( expdes, model_structure )
  best_res = Inf
  for (jj in sample){
    initial_response$local_response[new_link]=jj
    paramstmp = model$getParameterFromLocalResponse(initial_response$local_response, initial_response$inhibitors)
    tmp_result = model$fitmodel( data,paramstmp )
    if ( tmp_result$residuals < best_res ){
      best_res = tmp_result$residuals
      result = tmp_result
    }
  }
  response_matrix = model$getLocalResponseFromParameter( result$parameter )
  new_rank = model$modelRank()
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
                           1-pchisq(deltares, df=dr)),nrow=1)  
  colnames(extension_mat) <- c("adj_idx","from","to","value","residual","df","Res_delta","df_delta","pval")
  return(extension_mat)
}
