getModelStructure <- function(links, names=c(), struct_name="") {
  modelStructure =new(STASNet:::ModelStructure, as.character(links[,1]), as.character(links[,2]), as.character(names), struct_name)
  return(modelStructure);
}

getExperimentalDesign <- function(model.structure, stim.nodes, inhib.nodes, measured.nodes, stimuli, inhibitor, basal.activity = c()) {
    # Consistency checks
  if ( sum(is.na(match(as.character(stim.nodes),model.structure$names))) > 0 ) {
    stop(paste0("problem matching stim.nodes names, could not match ", stim.nodes[!stim.nodes%in%model.structure$names], collapse="\n") );
  }
  if ( sum(is.na(match(as.character(inhib.nodes),model.structure$names))) > 0 ) {
    stop(paste0("problem matching inhib.nodes names, could not match ", inhib.nodes[!inhib.nodes%in%model.structure$names], collapse="\n") );
  }
  if ( sum(is.na(match(as.character(measured.nodes),model.structure$names))) > 0 ) {
    stop(paste0("problem matching measured.nodes names, could not match ", measured.nodes[!measured.nodes%in%model.structure$names], collapse="\n") );
  }
  if ( sum(is.na(match(as.character(basal.activity),model.structure$names))) > 0 ) {
    message("Unmatched basal names :")
    message(paste0(basal.activity[is.na(match(as.character(basal.activity),model.structure$names))], collapse=" "))
  }
  if (length(basal.activity) != 0 && sum(basal.activity %in% model.structure$names) == 0 ) {
      warning("No basal names were matched!")
  }
  basal.activity = as.character(basal.activity[as.character(basal.activity) %in% model.structure$names])
  if ( !is.null(dim(stimuli)) && !is.null(dim(inhibitor)) &&  dim(stimuli)[1]!=dim(inhibitor)[1] ) {
    stop("number of experiments does not match in stimuli / inhibitor");
  } else if ( !is.null(dim(stimuli)) && dim(stimuli)[2]!=length(c(stim.nodes))) {
    stop("number stimuli does not match in stimuli / stim.nodes");
  } else if ( !is.null(dim(inhibitor)) && dim(inhibitor)[2]!=length(c(inhib.nodes)) ) {
    stop("number inhibitors does not match in inhib.nodes / inhibitor");
  }

  ExperimentalDesign <- STASNet:::ExperimentalDesign;
  expdes=new(ExperimentalDesign);

  expdes$stim_nodes=match(as.character(stim.nodes),model.structure$names)-1;
  expdes$inhib_nodes=match(as.character(inhib.nodes),model.structure$names)-1;
  expdes$measured_nodes=match(as.character(measured.nodes),model.structure$names)-1;
  expdes$basal_activity=model.structure$names %in% as.character(basal.activity);
  
  expdes$set_inhibitor(inhibitor);
  expdes$set_stimuli(stimuli);

  return(expdes)

}

#' Clone an MRAmodel or MRAmodelSet object 
#'
#' Copy a MRAmodel or MRAmodelSet object into a new independent variable
#' @param old_model A MRAmodel or MRAmodelSet object.
#' @return An MRAmodel/MRAmodelSet object with the same properties but separated from the old model
#' @export
#' @author Bertram Klinger \email{bertram.klinger@charite.de}
#' #' @examples \dontrun{
#' clonedModel = cloneModel(mramodel)
#' }
cloneModel <- function(old_model){
  
  type = class(old_model)
  
  if ("MRAmodelSet" %in% type ){
    model = new(STASNet:::ModelSet)
    data = new(STASNet:::DataSet)
    for (ii in 1:old_model$nb_models){
      data$addData(old_model$data$datas_list[[ii]], FALSE)
    }
  }else if ("MRAmodel" %in% type ){
    model = new(STASNet:::Model)
    data=new(STASNet:::Data)
  }else{
    stop(paste0("Wrong input class '",type,",' must be of class 'MRAmodel' or 'MRAmodelSet'!")) 
  }
  
  idx = which(old_model$structure$adjacencyMatrix==1)
  names = old_model$structure$names
  from = names[1+((idx-1) %/% length(names))]
  to = names[ifelse(idx %% length(names)==0, length(names), idx %% length(names))]
  links_list = cbind(from,to)
  structure = STASNet:::getModelStructure(links = links_list, struct_name = old_model$structure$title, names)
  design = STASNet:::getExperimentalDesign(model.structure = structure,
                                           stim.nodes = structure$names[old_model$design$stim_nodes+1],
                                           inhib.nodes = structure$names[old_model$design$inhib_nodes+1],
                                           measured.nodes = structure$names[old_model$design$measured_nodes+1],
                                           stimuli = old_model$design$stimuli,
                                           inhibitor = old_model$design$inhibitor,
                                           basal.activity = old_model$basal)

    model$setModel(design = design, structure = structure, old_model$use_log)
    
    data$set_unstim_data ( old_model$data$unstim_data )
    data$set_scale( old_model$data$scale )
    data$set_stim_data( old_model$data$stim_data )
    data$set_error( old_model$data$error)
    if (old_model$use_log) {
        data$use_log()
    }
    
  if ("MRAmodelSet" %in% type ){
    model$setNbModels(old_model$nb_models)
    new_model = STASNet:::MRAmodelSet(nb_models = old_model$nb_models,
                                      model = model,
                                      design = design,
                                      structure = structure,
                                      basal = old_model$basal,
                                      data = data,
                                      cv = old_model$cv,
                                      parameters = old_model$parameters,
                                      bestfit = old_model$bestfit,
                                      name = old_model$names,
                                      infos = old_model$infos,
                                      param_range = old_model$param_range,
                                      lower_values = old_model$lower_values,
                                      upper_values = old_model$upper_values,
                                      unused_perturbations = old_model$unused_perturbations,
                                      unused_readouts = old_model$unused_readouts,
                                      min_cv = old_model$min_cv,
                                      default_cv = old_model$default_cv,
                                      use_log = old_model$use_log)
      
    if ( length(old_model$variable_parameters) > 0 ){ 
      new_model = setVariableParameters(new_model, old_model$variable_parameters) 
    }
  }else{


    
      new_model = STASNet:::MRAmodel(model = model,
                                   design = design,
                                   structure = structure,
                                   basal = old_model$basal,
                                   data = data,
                                   cv = old_model$cv,
                                   parameters = old_model$parameters,
                                   bestfit = old_model$bestfit,
                                   name = old_model$name,
                                   infos = old_model$infos,
                                   param_range = old_model$param_range,
                                   lower_values = old_model$lower_values,
                                   upper_values = old_model$upper_values,
                                   unused_perturbations = old_model$unused_perturbations,
                                   unused_readouts = old_model$unused_readouts,
                                   min_cv = old_model$min_cv,
                                   default_cv = old_model$default_cv,
                                   use_log = old_model$use_log)
  }
  return(new_model)
}
